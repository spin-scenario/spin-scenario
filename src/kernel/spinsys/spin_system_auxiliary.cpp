/* Copyright 2019 The Spin-Scenario Authors. All Rights Reserved.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
    http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
==============================================================================*/

#include "spin_system_auxiliary.h"
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>
#include <functional>

namespace ssl {
namespace spinsys {

cx_vec acquire(const sp_cx_vec &rho0, const sp_cx_vec &coil, const sp_cx_mat &L, int points, double sw) {
  double dt = 1 / sw;
  cx_vec fid(points);
  sp_cx_vec rho = rho0;
  for (int i = 0; i < points; i++) {
    sp_cx_vec rho1 = step(rho, L, dt);
    fid[i] = projection(rho1, coil);
    rho = rho1;
  }
  return fid;

}

/* according to the motion equation in Liouville space, d_rho/d_t=(-iL0+R) rho = -i(L0+iR) rho => rho=expm[-i*(L+iR)*t]*rho.
Denote (L0+iR) as L, this function can be used in the following cases:
1. free evolution in the presence of relaxation with a duration time t (without pulse hamiltonian items).
2. small time step evolution with excitation pulse (with pulse hamiltonian items), which is frequently used in the optimal pulse design.
3. ideal evolution for a given flip angle \beta. \beta=2*_pi*u*t, thus here H should just be the pure pulse hamiltonian item, it is convenient
to manipulate the spin sys by _pi/2,_pi or any other flip angle.
*/
sp_cx_vec step(const sp_cx_vec &rho0, const sp_cx_mat &L, double dt) {
  sp_cx_vec rho = rho0;
  //////////////////////////////
  //rho=expm(-1i*L*dt)*rho;

  // Get the norm of the std::vector
  double norm_rho = norm(rho0);

  // Get the norm of the i*L*dt matrix
  double norm_mat = norm(L) * fabs(dt);
  // Scale the std::vector
  rho /= norm_rho;

  // Determine the number of time steps
  int nsteps = ceil(norm_mat / 5);

  sp_cx_vec next_term;
  // Run the Krylov procedure
  for (int i = 0; i < nsteps; i++) {
    next_term = rho;
    int k = 0;
    rho.setZero();

    while (max_abs(next_term) > std::numeric_limits<double>::epsilon()) {
      rho += next_term;
      k++;
      next_term = (-ci * dt / (double) (k * nsteps)) * (L * next_term);
    }

  }
  //// Scale the std::vector back
  rho *= norm_rho;
  return rho;
}

cx_vec step(const cx_vec &rho0, const cx_mat &L, double dt) {
  return (-ci * dt * L).exp() * rho0;
}

//sp_cx_vec step0(const sp_cx_vec &rho0, const sp_cx_mat &L, double dt)
//{
//  return (-ci * dt * L).exp() * rho0;
//}
// define function to be applied coefficient-wise
double ramp(double x) {
  if (x != 0)
    return x;
  else
    return g_inf;
}

sp_cx_mat expmv_tspan(const sp_cx_vec &b0,
                      const sp_cx_mat &A0,
                      double t0,
                      double tmax,
                      int q,
                      const mat &M0,
                      bool shift,
                      bool bal) {

  sp_cx_mat A = A0;
  sp_cx_vec b = b0;
  mat M = M0;
  double tol = 2 ^(-53); // double

  int n = A.rows();

//  cx_mat D;
//
//  if (bal) {
//    Eigen::ComplexEigenSolver<cx_mat> ces;
//    ces.compute(A);
//
//    cx_mat B = ces.eigenvalues().asDiagonal();
//
//    if (norm1(B) < norm1(A)) {
//      A = B;
//      D = ces.eigenvectors();
//      b = D.inverse() * b;
//
//    } else
//      bal = false;
//  }

  bool force_estm = 0;
  double temp = (tmax - t0) * norm1(A);

  if (temp > 63.152)
    force_estm = 1;

  //std::cout<<force_estm<<"force_estm\n";

  if (M.size() == 1)
    M = select_taylor_degree(A, shift, false, force_estm);

  sp_cx_mat X(n, q + 1);

  int s;
  int m;
  degree_selector(tmax - t0, M, m, s);
  double h = (tmax - t0) / (double) q;

  X.col(0) = expmv(b, A, 0, M, shift, false); //rho0

  // Easy case.
  if (q <= s) {
    std::cout << "Easy case\n";
    for (int k = 1; k <= q; k++) {
      X.col(k) = expmv(X.col(k - 1), A, h, M, shift, false);
    }

//    if (bal) {
//      X = D * X;
//    }
    return X;
  }

  cd mu = 0;
  if (shift) {
    mu = trace(A) / (double) n;
    A = A - mu * Identity(n);
  }

  int d = floor((double) q / (double) s);
  int j = floor((double) q / (double) d);
  int r = q - d * j;

  sp_cx_vec z = X.col(0);

  int m_opt;
  degree_selector((double) d, M, m_opt, s);

  int dr = d;

  for (int i = 1; i <= j + 1; i++) {
    if (i > j)
      dr = r;

    cx_mat K = cx_mat::Zero(n, m_opt + 1);
    K.col(0) = z;
    int m = 0;

    for (int k = 1; k <= dr; k++) {
      sp_cx_vec f = z;
      double c1 = norm_inf(z);

      int p;
      for (p = 1; p <= m_opt; p++) {
        if (p > m)
          K.col(p) = (h / (double) p) * (A * K.col(p - 1));

        cx_vec temp = pow(k, p) * K.col(p);
        f = f + temp;
        double c2 = norm_inf(temp);

        if ((c1 + c2) <= tol * norm_inf(f))
          break;
        c1 = c2;
      }
      m = std::max(m, p);
      X.col(k + (i - 1) * d) = exp((double) k * h * mu) * f;
    }

    if (i <= j)
      z = X.col(i * d);

  }

//  if (bal) {
//    X = D * X;
//  }
  return X;

}
cx_mat expmv_tspan(const cx_vec &b0,
                   const cx_mat &A0,
                   double t0,
                   double tmax,
                   int q,
                   const mat &M0,
                   bool shift,
                   bool bal) {
  cx_mat A = A0;
  cx_vec b = b0;
  mat M = M0;
  double tol = 2 ^(-53); // double

  int n = A.rows();

  cx_mat D;

  if (bal) {
    Eigen::ComplexEigenSolver<cx_mat> ces;
    ces.compute(A);

    cx_mat B = ces.eigenvalues().asDiagonal();

    if (norm1(B) < norm1(A)) {
      A = B;
      D = ces.eigenvectors();
      b = D.inverse() * b;

    } else
      bal = false;
  }

  bool force_estm = 0;
  double temp = (tmax - t0) * norm1(A);

  if (temp > 63.152)
    force_estm = 1;

  //std::cout<<force_estm<<"force_estm\n";

  if (M.size() == 1)
    M = select_taylor_degree(A, shift, false, force_estm);

  cx_mat X = cx_mat::Zero(n, q + 1);

  int s;
  int m;
  degree_selector(tmax - t0, M, m, s);
  double h = (tmax - t0) / (double) q;

  X.col(0) = expmv(b, A, 0, M, shift, false); //rho0

  // Easy case.
  if (q <= s) {
    std::cout << "Easy case\n";
    for (int k = 1; k <= q; k++) {
      X.col(k) = expmv(X.col(k - 1), A, h, M, shift, false);
    }

    if (bal) {
      X = D * X;
    }
    return X;
  }

  cd mu = 0;
  if (shift) {
    mu = A.trace() / (double) n;
    A = A - mu * cx_mat::Identity(n, n);
  }

  int d = floor((double) q / (double) s);
  int j = floor((double) q / (double) d);
  int r = q - d * j;

  cx_vec z = X.col(0);

  int m_opt;
  degree_selector((double) d, M, m_opt, s);

  int dr = d;

  for (int i = 1; i <= j + 1; i++) {
    if (i > j)
      dr = r;

    cx_mat K = cx_mat::Zero(n, m_opt + 1);
    K.col(0) = z;
    int m = 0;

    for (int k = 1; k <= dr; k++) {
      cx_vec f = z;
      double c1 = norm_inf(z);

      int p;
      for (p = 1; p <= m_opt; p++) {
        if (p > m)
          K.col(p) = (h / (double) p) * (A * K.col(p - 1));

        cx_vec temp = pow(k, p) * K.col(p);
        f = f + temp;
        double c2 = norm_inf(temp);

        if ((c1 + c2) <= tol * norm_inf(f))
          break;
        c1 = c2;
      }
      m = std::max(m, p);
      X.col(k + (i - 1) * d) = exp((double) k * h * mu) * f;
    }

    if (i <= j)
      z = X.col(i * d);

  }

  if (bal) {
    X = D * X;
  }
  return X;
}

sp_cx_mat Identity(int n) {
  sp_cx_mat sm(n, n);
  for (int i = 0; i < n; i++)
    sm.coeffRef(i, i) = 1;
  return sm;
}

sp_cx_vec expmv(const sp_cx_vec &b0, const sp_cx_mat &A0, double t, const mat &M0, bool shift, bool bal) {
  sp_cx_mat A = A0;
  sp_cx_vec b = b0;
  mat M = M0;
  double tol = 2 ^(-53);

//  cx_mat D;
//
//  if (bal) {
//    Eigen::ComplexEigenSolver<cx_mat> ces;
//    ces.compute(A);
//
//    cx_mat B = ces.eigenvalues().asDiagonal();
//
//    if (norm1(B) < norm1(A)) {
//      A = B;
//      D = ces.eigenvectors();
//      b = D.inverse() * b;
//
//    } else
//      bal = false;
//  }

  int n = A.rows();

  cd mu;
  if (shift) {
    mu = trace(A) / (double) n;
    A = A - mu * Identity(n);
  }

  double tt;
  if (M.size() == 1) {
    sp_cx_mat C = t * A;
    M = select_taylor_degree(C, false, false, false);
    //std::cout<<M<<"\n\n";
    tt = 1;
  } else
    tt = t;

  int s = 1;
  int m;
  if (t == 0) // The case tA = 0
    m = 0;
  else
    degree_selector(tt, M, m, s);

  cd eta = 1;
  if (shift) {
    eta = exp(t * mu / (double) s);
  }

  sp_cx_vec f = b;

  for (int i = 1; i <= s; i++) {
    double c1 = norm_inf(b);

    for (int k = 1; k <= m; k++) {

      b = (t / (s * k)) * (A * b);
      f = f + b;
      double c2 = norm_inf(b);

      if ((c1 + c2) <= tol * norm_inf(f))
        break;
      c1 = c2;
    }

    f = eta * f;
    b = f;
  }

  if (bal) {
    //f = D * f;
  }

  return f;
}
cx_vec expmv(const cx_vec &b0, const cx_mat &A0, double t, const mat &M0, bool shift, bool bal) {
  cx_mat A = A0;
  cx_vec b = b0;
  mat M = M0;
  double tol = 2 ^(-53);

  cx_mat D;

  if (bal) {
    Eigen::ComplexEigenSolver<cx_mat> ces;
    ces.compute(A);

    cx_mat B = ces.eigenvalues().asDiagonal();

    if (norm1(B) < norm1(A)) {
      A = B;
      D = ces.eigenvectors();
      b = D.inverse() * b;

    } else
      bal = false;
  }

  int n = A.rows();

  cd mu;
  if (shift) {
    mu = A.trace() / (double) n;
    A = A - mu * cx_mat::Identity(n, n);
  }

  double tt;
  if (M.size() == 1) {
    M = select_taylor_degree(t * A, false, false, false);
    //std::cout<<M<<"\n\n";
    tt = 1;
  } else
    tt = t;

  int s = 1;
  int m;
  if (t == 0) // The case tA = 0
    m = 0;
  else
    degree_selector(tt, M, m, s);

  cd eta = 1;
  if (shift) {
    eta = exp(t * mu / (double) s);
  }

  cx_vec f = b;

  for (int i = 1; i <= s; i++) {
    double c1 = norm_inf(b);

    for (int k = 1; k <= m; k++) {

      b = (t / (s * k)) * (A * b);
      f = f + b;
      double c2 = norm_inf(b);

      if ((c1 + c2) <= tol * norm_inf(f))
        break;
      c1 = c2;
    }

    f = eta * f;
    b = f;
  }

  if (bal) {
    f = D * f;
  }

  return f;
}

double norm1(const cx_mat &m) {
  return m.cwiseAbs().colwise().sum().maxCoeff();
}
double norm_inf(const cx_mat &m) {
  return m.cwiseAbs().rowwise().sum().maxCoeff();
}

double norm1(const sp_cx_mat &m) {
  std::vector<double> sum_cols;
  int ncols = m.cols();
  sp_mat m_abs = m.cwiseAbs();

  for (int i = 0; i < ncols; i++)
    sum_cols.push_back(m_abs.col(i).sum());
  return *max_element(sum_cols.begin(), sum_cols.end());
}

double norm_inf(const sp_cx_vec &m) {
  std::vector<double> sum_rows;
  int nrows = m.rows();
  sp_vec m_abs = m.cwiseAbs();

  for (int i = 0; i < nrows; i++)
    sum_rows.push_back(m_abs.coeff(i));
  return *max_element(sum_rows.begin(), sum_rows.end());
}
void degree_selector(double t, const mat &M, int &m, int &s) {
  int m_max = M.rows();
  int p = M.cols();

  vec v = vec::LinSpaced(m_max, 1, m_max);
  mat C = (abs(t) * M).array().ceil().matrix().transpose() * v.asDiagonal();
  C = C.unaryExpr([](double x) {return ramp(x);});
  double cost;
  if (p > 1) {
    cost = (C.colwise().minCoeff()).minCoeff(&m);
  } else {
    //when C is one column. Happens if p_max = 2.
    cost = C.col(0).minCoeff(&m);
  }

  if (cost == g_inf)
    cost = 0;
  s = std::max<int>((int) (cost / (double) m), 1);
}

mat select_taylor_degree(const cx_mat &A0, bool shift, bool bal, bool force_estm) {
  cx_mat A = A0;

  int p_max = 8;
  int m_max = 55;

  int n = A.rows();

  if (bal) {
    Eigen::ComplexEigenSolver<cx_mat> ces;
    ces.compute(A);

    cx_mat B = ces.eigenvalues().asDiagonal();
    if (norm1(B) < norm1(A)) {
      A = B;
    }
  }

  if (shift) {
    cd mu = A.trace() / (double) n;
    A = A0 - mu * cx_mat::Identity(n, n);
  }

  vec alpha;

  if (!force_estm) {
    double normA = norm1(A);
    alpha = vec::Ones(p_max - 1);
    //Base choice of m on normA, not the alpha_p.
    if ((normA <= 4.0 * g_expmv_theta(m_max - 1) * (double) (p_max * (p_max + 3)) / (double) (m_max)))
      alpha *= normA;
  } else {
    vec eta = vec::Zero(p_max);
    alpha = vec::Zero(p_max - 1);
    for (int p = 1; p <= p_max; p++) {
      cx_mat tmp = A.pow(1 + p);
      double c = norm1(tmp);
      eta(p - 1) = pow(c, 1.0 / (double) (p + 1));
    }

    for (int p = 0; p < p_max - 1; p++)
      alpha(p) = std::max<double>(eta(p), eta(p + 1));
  }

  mat M = mat::Zero(m_max, p_max - 1);

  for (int p = 2; p <= p_max; p++)
    for (int m = p * (p - 1) - 1; m <= m_max; m++)
      M(m - 1, p - 2) = alpha(p - 2) / g_expmv_theta(m);
  return M;
}

mat select_taylor_degree(const sp_cx_mat &A0, bool shift, bool bal, bool force_estm) {
  sp_cx_mat A = A0;

  int p_max = 8;
  int m_max = 55;

  int n = A.rows();

//  if (bal) {
//    Eigen::ComplexEigenSolver<cx_mat> ces;
//    ces.compute(A);
//
//    cx_mat B = ces.eigenvalues().asDiagonal();
//    if (norm1(B) < norm1(A)) {
//      A = B;
//    }
//  }

  if (shift) {
    cd mu = trace(A) / (double) n;
    A = A0 - mu * Identity(n);
  }

  vec alpha;
  if (!force_estm) {
    double normA = norm1(A);
    alpha = vec::Ones(p_max - 1);
    //Base choice of m on normA, not the alpha_p.
    if ((normA <= 4.0 * g_expmv_theta(m_max - 1) * (double) (p_max * (p_max + 3)) / (double) (m_max)))
      alpha *= normA;
  } else {
    vec eta = vec::Zero(p_max);
    alpha = vec::Zero(p_max - 1);
    for (int p = 1; p <= p_max; p++) {
      sp_cx_mat tmp = spow(A, 1 + p);
      double c = norm1(tmp);
      eta(p - 1) = pow(c, 1.0 / (double) (p + 1));
    }

    for (int p = 0; p < p_max - 1; p++)
      alpha(p) = std::max<double>(eta(p), eta(p + 1));
  }

  mat M = mat::Zero(m_max, p_max - 1);

  for (int p = 2; p <= p_max; p++)
    for (int m = p * (p - 1) - 1; m <= m_max; m++)
      M(m - 1, p - 2) = alpha(p - 2) / g_expmv_theta(m);
  return M;
}

cd projection(const sp_cx_vec &rho, const sp_cx_vec &coil) {
  return (coil.adjoint() * rho).toDense().trace();
}

cd projection(const cx_vec &rho, const cx_vec &coil) {
  return (coil.adjoint() * rho).trace();
}

double transfer_fidelity(const sp_cx_vec &rho, const sp_cx_vec &targ) {
  return traced(targ.adjoint() * rho).real() - 1;
  //cd val = traced(rho.adjoint() * targ)-cd(1,0);
  //return fabs(val);
}
double transfer_fidelity(const sp_cx_mat &rho, const sp_cx_mat &targ) {
  return traced(targ.adjoint() * rho).real();
  //cd val= traced(rho.adjoint() * targ);
  //return 0.5*abs(val) * abs(val);

  // cd val = traced(rho.adjoint() * targ)-cd(1,0);
  // return fabs(val);
}
sp_cx_vec norm_state(const sp_cx_vec &rho) {
  sp_cx_vec rho_new = rho;
  double norm = rho.norm();
  rho_new = (1 / norm) * rho_new;
  // enable Levante-Ernst thermalization.
  rho_new.coeffRef(0) = 1;
  return rho_new;
}

void levante_ernst_correction(sp_cx_vec &rho) {
  rho.coeffRef(0) = 1;
}

void levante_ernst_correction(cx_vec &rho) {
  rho[0] = 1;
}

double max_abs(const sp_cx_vec &m) {
  double max_val = std::numeric_limits<double>::min();
  double v;
  for (sp_cx_vec::InnerIterator it(m); it; ++it) {
    v = sqrt(std::norm(it.value()));
    if (v > max_val)
      max_val = v;
  }
  return max_val;
}
sp_cx_mat commutator(const sp_cx_mat &op1, const sp_cx_mat &op2) {
  return op1 * op2 - op2 * op1;
}
sp_cx_mat propagator(const sp_cx_mat &L, double dt) {
  sp_cx_mat A = -ci * dt * L;
  double mat_norm = norm(A);
  // Warn the user if the norm is too big
  if (mat_norm > 1024) {
  }
  // Determine scaling and squaring parameters
  double n_squarings = std::max(0.0, ceil(log2(mat_norm)));
  double scaling_factor = std::pow(2, n_squarings);

  // Scale the matrices
  if (scaling_factor > 1)
    A /= scaling_factor;
  //  Get the propagator by 'taylor'
  // Run the Taylor series procedure for the propagator
  int size = L.rows();
  sp_cx_mat P(size, size);
  P.setIdentity();
  sp_cx_mat next_term(size, size);
  next_term.setIdentity();
  double next_term_norm = std::numeric_limits<double>::max();
  double n = 1;
  while (next_term_norm > 1e-9) {
    // Compute the next term
    next_term = cd(1 / n, 0) * next_term * A;

    // Eliminate small elements
    cleanup(next_term, 1e-10);

    // Compute the residual norm
    sp_mat t = next_term.cwiseAbs();
    next_term_norm = maxm(t);

    //Add to the total and increment the counter
    P += next_term;
    n += 1;
  }
  // Run the squaring stage
  for (int i = 1; i <= n_squarings; i++) {
    //Square the propagator
    P = P * P;
    cleanup(P, 1e-10);
  }
  return P;
}

composition::composition() : B0_(0) {
}

composition::~composition() {
}

double composition::get_proton_freq() const {
  isotope proton("1H");
  double MHz = B0_ * proton.gamma() / (2 * _pi * 1e6);
  return MHz;
}
void composition::init() {
  basis_table_ = basis_table();

  // fundamental zeeman frequencies
  base_freq_.resize(nspins());
  for (size_t i = 0; i < spins_.size(); i++)
    base_freq_[i] = B0_ * spins_[i].gamma();  // Note magnet should be set before.
  // The sign of freq to be confirmed.
}
std::vector<size_t> composition::parse_labeled_spins(const std::string symbol) const {
  std::vector<size_t> idx;
  for (size_t i = 0; i < spins_.size(); i++)
    if (symbol == spins_[i].symbol())
      idx.push_back(i);
  return idx;
}
imat composition::basis_table() const {
  size_t n = spins_.size();
  size_t dim = hilbert_space_dim();
  imat table = imat::Zero(dim * dim, n);
  for (size_t i = 0; i < n; i++) {
    size_t dim_cur = hilbert_space_dim(i);
    size_t dim_front = hilbert_space_dim_forward(i);
    size_t dim_back = hilbert_space_dim_backward(i);

    ivec front = ivec::Ones(dim_front * dim_front);
    ivec cur = ivec::Ones(dim_cur * dim_cur);
    for (size_t j = 0; j < dim_cur * dim_cur; j++)
      cur(j) = j;

    ivec back = ivec::Ones(dim_back * dim_back);
    table.col(i) = kroneckerProduct(kroneckerProduct(front, cur).eval(),
                                    back).eval();
  }
  return table;
  // Todo: spin_system.bas.basis=unique(basis_spec,'rows');
}
cx_mat composition::basis_state() const {
  cx_mat ss(basis_table_.rows(), basis_table_.cols());
  for (int i = 0; i < basis_table_.rows(); i++)
    for (int j = 0; j < basis_table_.cols(); j++) {
      Tlm tmp(basis_table_(i, j));
      ss(i, j) = cd(tmp.L, tmp.M);
    }
  return ss;
}
mat composition::freq_basis() const {
  mat freqs = mat::Ones(basis_table_.rows(), basis_table_.cols());
  for (size_t i = 0; i < nspins(); i++)
    freqs.col(i) *= base_freq_[i];
  return freqs;
}

double composition::get_freq(size_t id) const {
  double MHz = base_freq_[id] / (2 * _pi * 1e6);
  return MHz;
}

size_t composition::hilbert_space_dim() const {
  size_t m = 1;
  for (size_t i = 0; i < spins_.size(); i++)
    m *= spins_[i].hs();
  return m;
}
size_t composition::hilbert_space_dim_forward(size_t id) const {
  size_t m = 1;
  for (size_t i = 0; i < id; i++)
    m *= spins_[i].hs();
  return m;
}
size_t composition::hilbert_space_dim_backward(size_t id) const {
  size_t m = 1;
  for (size_t i = id + 1; i < spins_.size(); i++)
    m *= spins_[i].hs();
  return m;
}
sp_cx_mat composition::mprealloc() const {
  size_t dim = liouville_space_dim();
  return sp_cx_mat(dim, dim);
}
sp_cx_vec composition::state(const std::string list) const {
  return op(list, kLeft).cast<cd>() * unit_state();
}
sp_mat composition::op(const std::string list, op_side type) const {
  std::string list_copy = list;
  boost::trim(list_copy);
  std::vector<std::string> par_vec;
  boost::split(par_vec, list_copy, boost::is_any_of("\t, "), boost::token_compress_on);
  std::string op_type = par_vec.back();
  if (op_type == "comm") {
    type = kComm;
    par_vec.pop_back();
  } else if (op_type == "lcomm") {
    type = kLeftComm;
    par_vec.pop_back();
  } else if (op_type == "rcomm") {
    type = kRightComm;
    par_vec.pop_back();
  } else if (op_type == "left") {
    type = kLeft;
    par_vec.pop_back();
  } else if (op_type == "right") {
    type = kRight;
    par_vec.pop_back();
  } else if (op_type == "acomm") {
    type = kAcomm;
    par_vec.pop_back();
  }

  bool is_product_superop = false;
  try {
    boost::lexical_cast<size_t>(par_vec[0]);
    is_product_superop = true;
  } catch (boost::bad_lexical_cast &) {
    // if it throws, it's not a number.
  }
  std::vector<int> flags = std::vector<int>(nspins(), 0);
  std::vector<op_coeff> all;
  double coeff = 1;
  if (is_product_superop) { // op1 case.
    for (size_t i = 0; i < par_vec.size(); i++) {
      size_t id = boost::lexical_cast<size_t>(par_vec[i]);
      std::string op_label = par_vec[++i];

      check_valid_spin_id(id);
      check_valid_op_label(op_label);

      op_coeff tmp;
      tmp.id = id - 1;
      boost::to_lower(op_label);
      if (op_label == "ie") {
        tmp.linear_index = kIe;
        tmp.coeff = 1;
      }
      if (op_label == "i+" || op_label == "ip") {
        tmp.linear_index = kIp;
        tmp.coeff = -sqrt(2);
      }
      if (op_label == "iz") {
        tmp.linear_index = kIz;
        tmp.coeff = 1;
      }
      if (op_label == "i-" || op_label == "im") {
        tmp.linear_index = kIm;
        tmp.coeff = sqrt(2);
      }
      if (op_label == "t") {  // Tl,m type.
        Tlm tensor;
        tensor.L = boost::lexical_cast<int>(par_vec[++i]);
        tensor.M = boost::lexical_cast<int>(par_vec[++i]);
        tmp.linear_index = tensor.linear_indexing();
        tmp.coeff = 1;
      }
      flags[tmp.id] = tmp.linear_index;
      coeff *= tmp.coeff;
      all.push_back(tmp);
    }
  } else {
    std::vector<size_t> spins_id = parse_labeled_spins(par_vec[0]);
    if (spins_id.size() == 0) {
      std::string s = "invalid spin symbol ** " + par_vec[0] + " **";
      throw std::runtime_error(s.c_str());
    }
    std::string op_label = par_vec.back();
    check_valid_op_label(op_label);
    boost::to_lower(op_label);
    for (size_t i = 0; i < spins_id.size(); i++) {
      op_coeff tmp;
      tmp.id = spins_id[i];

      if (op_label == "ie") {
        tmp.linear_index = kIe;
        tmp.coeff = 1;
      }
      if (op_label == "i+" || op_label == "ip") {
        tmp.linear_index = kIp;
        tmp.coeff = -sqrt(2);
      }
      if (op_label == "iz") {
        tmp.linear_index = kIz;
        tmp.coeff = 1;
      }
      if (op_label == "i-" || op_label == "im") {
        tmp.linear_index = kIm;
        tmp.coeff = sqrt(2);
      }
      if (op_label == "t") {  // Tl,m type.
        Tlm tensor;
        tensor.L = boost::lexical_cast<int>(par_vec[1]);
        tensor.M = boost::lexical_cast<int>(par_vec[2]);
        tmp.linear_index = tensor.linear_indexing();
        tmp.coeff = 1;
      }
      flags[tmp.id] = tmp.linear_index;
      all.push_back(tmp);
    }

  }
  //BasisType formalism = basis_.type;
  //switch (formalism) {
  //case sphten_liouv: {
  sp_mat m = empty_op();
  if (!is_product_superop)
    for (size_t k = 0; k < all.size(); k++) {
      std::vector<op_coeff> tmp;
      tmp.push_back(all[k]);
      m += all[k].coeff * p_superop(tmp, flags, type);
    }
  else
    m = coeff * p_superop(all, flags, type);
  return m;
  //}
  //break;
  //case zeeman_hilb: {
  //	sp_mat m;
  //	for (size_t i = 0; i < flags.size(); i++) {
  //		LM lm = linear_index(flags[i]);
  //		std::vector<sp_mat> T = irr_sph_ten2(getHilbertSpaceDim(i), lm.L);
  //		if (i == 0)
  //			m = T[lm.L - lm.M];
  //		else {
  //			sp_mat tmp = kroneckerProduct(m, T[lm.L - lm.M]);
  //			m = tmp;   // to be improved!
  //		}
  //	}
  //	return coeff * m;
  //}
  //				  break;
  //default:
  //        break;
  //}
}
void composition::check_valid_spin_id(const size_t id) const {
  size_t n = nspins();
  if (id > n || id == 0) {
    std::string s = "invalid spin index ** " + boost::lexical_cast<std::string>
        (id) + " ** ( 1 ~ " + boost::lexical_cast<std::string>(n) + " )";
    throw std::runtime_error(s.c_str());
  }
}
void composition::check_valid_op_label(const std::string linear_index) const {
  if (linear_index != "Ie" && linear_index != "Iz" && linear_index != "I+"
      && linear_index != "Ip" && linear_index != "I-" && linear_index != "Im"
      && linear_index != "T") {
    std::string s = "invalid operator notation ** " + linear_index + " **";
    throw std::runtime_error(s.c_str());
  }
}

sp_mat composition::empty_op() const {
  size_t dim = liouville_space_dim();
  return sp_mat(dim, dim);
}
sp_cx_vec composition::unit_state() const {
  //BasisType formalism = sys.basis_type();
  //  switch (formalism) {
  //  case sphten_liouv:
  size_t dim = liouville_space_dim();
  sp_cx_vec unit(dim);
  unit.coeffRef(0) = 1.0;
  return unit;
  //  break;
  //case zeeman_hilb:
  //  break;
  //default:
  //  break;
  //}
}
sp_mat composition::p_superop(const std::vector<op_coeff> &active_spins,
                              const std::vector<int> &flags,
                              const op_side type) const {

  if (type != kComm && type != kLeftComm && type != kRightComm && type != kLeft
      && type != kRight) {
    //throw ErrorInfo("unknown commutation relation for operator setup!");
  }

  if (type == kComm)
    return p_superop(active_spins, flags, kLeftComm)
        - p_superop(active_spins, flags, kRightComm);
  size_t active_num = active_spins.size();
  std::vector<ivec> source = std::vector<ivec>(active_num);  // source state index
  std::vector<ivec> destin = std::vector<ivec>
      (active_num);  // destination state index
  std::vector<vec> structure = std::vector<vec>
      (active_num);  // structure coefficients table
  for (size_t i = 0; i < active_num; i++) {
    size_t mult = hilbert_space_dim(active_spins[i].id);
    size_t mult2 = mult * mult;
    mat *pt_left = new mat[mult2];
    mat *pt_right = new mat[mult2];
    ist_product_table(mult, pt_left, pt_right);

    mat pt = mat::Zero(mult2, mult2);
    if (type == kLeftComm || type == kLeft) {
      for (size_t k = 0; k < mult2; k++)
        pt.col(k) = pt_left[k].row(flags[active_spins[i].id]);
    } else if (type == kRightComm || type == kRight) {
      for (size_t k = 0; k < mult2; k++)
        pt.col(k) = pt_right[k].row(flags[active_spins[i].id]);
    }

    std::vector<size_t> rows;
    std::vector<size_t> cols;
    std::vector<double> vals;
    for (size_t col = 0; col < mult2; col++)
      for (size_t row = 0; row < mult2; row++) {
        double val = pt(row, col);
        if (val) {
          rows.push_back(row);
          cols.push_back(col);
          vals.push_back(val);
        }
      }
    destin[i] = ivec(rows.size());
    source[i] = ivec(rows.size());
    structure[i] = vec(rows.size());

    for (size_t p = 0; p < rows.size(); p++) {
      destin[i](p) = rows[p];
      source[i](p) = cols[p];
      structure[i](p) = vals[p];
    }
  }
  // Compute the structure coefficients for the relevant sub-algebra
  imat from = source[0];
  imat to = destin[0];
  mat coeff = structure[0];
  for (size_t i = 1; i < active_num; i++) {
    ivec ones1;
    ones1.setOnes(source[i].rows());
    ivec ones2;
    ones2.setOnes(from.rows());
    imat from1 = kroneckerProduct(from, ones1);
    imat from2 = kroneckerProduct(ones2, source[i]);
    from.resize(from1.rows(), from1.cols() + from2.cols());
    from.block(0, 0, from1.rows(), from1.cols()) = from1;
    from.block(0, from1.cols(), from2.rows(), from2.cols()) = from2;

    ones1.setOnes(destin[i].rows());
    ones2.setOnes(to.rows());
    imat to1 = kroneckerProduct(to, ones1);
    imat to2 = kroneckerProduct(ones2, destin[i]);
    to.resize(to1.rows(), to1.cols() + to2.cols());
    to.block(0, 0, to1.rows(), to1.cols()) = to1;
    to.block(0, to1.cols(), to2.rows(), to2.cols()) = to2;

    coeff = kroneckerProduct(coeff, structure[i]).eval();
  }
  //std::cout << "from:\n" << from << "\n\n" << "to:\n" << to << "\n\n" << "coeff:\n"
  //    << coeff << "\n";
  size_t dim = liouville_space_dim();
  sp_mat mat(dim, dim);

  // Lift the basis columns corresponding to the relevant spins
  std::vector<size_t> active_ids;
  for (size_t k = 0; k < active_spins.size(); k++)
    active_ids.push_back(active_spins[k].id);
  imat basis_col = basis_cols(active_ids);

  std::vector<int> used_rows;
  // For commutation superoperators remove commuting paths
  if (type == kLeftComm || type == kRightComm) {
    for (int i = 0; i < from.rows(); i++)
      if (from.row(i).sum() != 0 && to.row(i).sum() != 0)
        used_rows.push_back(i);
  } else {
    for (int i = 0; i < from.rows(); i++)
      used_rows.push_back(i);
  }

  // Loop over source states
  for (size_t i = 0; i < used_rows.size(); i++) {
    // Retrieve the source subspace
    ivec source_subsp_idx;
    source_subsp_idx.setOnes(basis_col.rows());
    ivec comp(basis_col.rows());
    for (int j = 0; j < from.cols(); j++) {
      comp.setZero();
      for (int k = 0; k < basis_col.rows(); k++)
        if (basis_col(k, j) == from(used_rows[i], j))
          comp(k) = 1;
      source_subsp_idx = source_subsp_idx.cwiseProduct(comp);
    }

    std::vector<int> source_row_select;
    for (int j = 0; j < source_subsp_idx.rows(); j++)
      if (source_subsp_idx(j))
        source_row_select.push_back(j);
    imat source_subsp = basis_cols(active_ids, source_row_select);

    // Retrieve the destination subspace
    ivec destin_subsp_idx;
    destin_subsp_idx.setOnes(basis_col.rows());
    for (int j = 0; j < to.cols(); j++) {
      comp.setZero();
      for (int k = 0; k < basis_col.rows(); k++)
        if (basis_col(k, j) == to(used_rows[i], j))
          comp(k) = 1;
      destin_subsp_idx = destin_subsp_idx.cwiseProduct(comp);
    }

    std::vector<int> destin_row_select;
    for (int j = 0; j < destin_subsp_idx.rows(); j++)
      if (destin_subsp_idx(j))
        destin_row_select.push_back(j);
    imat destin_subsp = basis_cols(active_ids, destin_row_select);

    // Fill the operator  isequal(source_subsp, destin_subsp)
    if (source_subsp == destin_subsp) {
      // If the subspaces fully match, use the raw indices
      int subsp_dim = source_subsp.rows();
      for (int k = 0; k < subsp_dim; k++) {
        int a = source_row_select[k];
        int b = destin_row_select[k];
        if (mat.coeff(a, b))
          mat.coeffRef(a, b) += coeff(used_rows[i], 0);
        else
          mat.insert(a, b) = coeff(used_rows[i], 0);
      }
      //mat.makeCompressed();
    } else { // Otherwise, use brute-force state-by-state matching

      Lia tt = ismember(source_subsp, destin_subsp);
      ivec does_it_go_anywhere = tt.lta;
      ivec where_it_goes_if_it_does = tt.locb;

      for (int k = 0; k < does_it_go_anywhere.size(); k++) {
        if (does_it_go_anywhere[k]) {
          int a = source_row_select[k];
          int b = destin_row_select[where_it_goes_if_it_does[k] - 1];
          if (mat.coeff(a, b))
            mat.coeffRef(a, b) += coeff(used_rows[i], 0);
          else
            mat.insert(a, b) = coeff(used_rows[i], 0);
        }
      }
    }
  }
  return mat;
}
// Structure coefficient tables for the associateive envelopes of su(mult) algebras
void composition::ist_product_table(size_t mult, mat *product_table_left,
                                    mat *product_table_right) const {
  // Get the irreducible spherical tensors
  std::vector<sp_mat> T = irr_sph_ten2(mult);

  // Preallocate the arrays
  size_t mult2 = mult * mult;
  for (size_t i = 0; i < mult2; i++) {
    product_table_left[i] = mat::Zero(mult2, mult2);
    product_table_right[i] = product_table_left[i];
  }

  // Get the structure coefficients
  sp_mat a, b;
  for (size_t m = 0; m < mult2; m++)
    for (size_t k = 0; k < mult2; k++) {
      a = T[k] * T[k].transpose();
      b = T[m] * T[m].transpose();
      double norm = sqrt(trace(a) * trace(b));
      for (size_t n = 0; n < mult2; n++) {
        a = T[n] * T[m] * T[k].transpose();
        b = T[m] * T[n] * T[k].transpose();
        product_table_left[k](n, m) = trace(a) / norm;
        product_table_right[k](n, m) = trace(b) / norm;
      }
    }
}
std::vector<sp_mat> composition::irr_sph_ten2(size_t mult) const {
  std::vector<sp_mat> T;
  for (int k = 0; k < (int) mult; k++) {
    std::vector<sp_mat> Tk = irr_sph_ten2(mult, k);
    for (size_t i = 0; i < Tk.size(); i++)
      T.push_back(Tk[i]);
  }
  return T;
}
std::vector<sp_mat> composition::irr_sph_ten2(size_t mult, int k) const {
  std::vector<sp_mat> Tk;
  pauli_matrices L(mult);
  // Get the top state
  sp_mat T0 = (std::pow(-1, k) * std::pow(2, -k / 2.0)) * spow(L.p, k);
  Tk.push_back(T0);
  for (int j = 0; j <= 2 * k - 1; j++) {
    int q = k - j;
    sp_mat Tj = (1.0 / sqrt(double((k + q) * (k - q + 1))))
        * (L.m * Tk[j] - Tk[j] * L.m);
    Tk.push_back(Tj);
  }
  return Tk;
}
imat composition::basis_cols(std::vector<size_t> active_ids) const {
  size_t cols = active_ids.size();
  imat m(basis_table_.rows(), cols);
  for (size_t i = 0; i < cols; i++)
    m.col(i) = basis_table_.col(active_ids[i]);
  return m;
}
imat composition::basis_cols(std::vector<size_t> active_ids,
                             std::vector<int> row_select) const {
  size_t cols = active_ids.size();
  size_t rows = row_select.size();

  size_t spin_num = nspins();
  imat n(rows, spin_num - cols);
  if (cols == spin_num || rows == 0)
    return n;

  imat m(basis_table_.rows(), spin_num - cols);
  std::vector<int> if_active(spin_num, 0);
  for (size_t i = 0; i < cols; i++)  // set active flags (1)
    if_active[active_ids[i]] = 1;

  size_t j = 0;
  for (size_t i = 0; i < spin_num; i++) {
    if (!if_active[i])
      m.col(j++) = basis_table_.col(i);
  }

  for (size_t i = 0; i < rows; i++)
    n.row(i) = m.row(row_select[i]);

  return n;
}

///////////////////////////////////////////////////////////////////////////////
interaction::interaction(
    const composition &comp) : comp_(comp), assume_("") {
}

interaction::~interaction() {
}

void interaction::set_zeeman(std::string list) {
  boost::trim(list); // trim the spaces on both left and right sides.
  boost::to_lower(list);

  std::vector<std::string> par_vec;
  boost::split(par_vec, list, boost::is_any_of("\t, "), boost::token_compress_on);

  std::vector<std::string>::iterator iter;
  iter = par_vec.begin();
  while (iter != par_vec.end()) {
    if (*(iter + 1) == "scalar") {
      size_t id = boost::lexical_cast<size_t>(*iter);

      std::string val = *(iter + 2);
      std::vector<std::string> val_vec;
      boost::split(val_vec, val, boost::is_any_of(":"), boost::token_compress_on);

      std::string unit = *(iter + 3);

      if (val_vec.size() == 1) // normal case.
        set_zeeman_scalar(id, boost::lexical_cast<double>(val), unit);
      else if (val_vec.size() == 3) { // broadband case.
        vec3 par;
        par[0] = boost::lexical_cast<double>(val_vec[0]);
        par[1] = boost::lexical_cast<double>(val_vec[1]);
        par[2] = boost::lexical_cast<double>(val_vec[2]);
        //set_zeeman_scalar(id, par[0], unit); // this is the default cs value.
        set_zeeman_scalar(id, (par[0]+par[1])/2, unit); // this is the default cs value.
        set_zeeman_scalar_broadband(id, par, unit);
      }
      iter += 4;
    } else if (*(iter + 1) == "eigs-euler") {
      // TODO.
    } else {
      std::string s = "unidentified zeeman parameters: " + list;
      throw std::runtime_error(s.c_str());
    }
  }
#ifdef SSL_OUTPUT_ENABLE
  std::string s = str(boost::format("%s %s.\n") % "zeeman effect pars set to be" % list);
  ssl_color_text("info", s);
#endif
}
void interaction::set_Jcoupling(std::string list) {
  boost::trim(list); // trim the spaces on both left and right sides.
  boost::to_lower(list);

  std::vector<std::string> par_vec;
  boost::split(par_vec, list, boost::is_any_of("\t, "), boost::token_compress_on);

  std::vector<std::string>::iterator iter;
  iter = par_vec.begin();
  while (iter != par_vec.end()) {
    if (*(iter + 2) == "scalar") {
      size_t id1 = boost::lexical_cast<size_t>(*iter);
      size_t id2 = boost::lexical_cast<size_t>(*(iter + 1));

      std::string val = *(iter + 3);
      std::vector<std::string> val_vec;
      boost::split(val_vec, val, boost::is_any_of(":"), boost::token_compress_on);

      std::string unit = *(iter + 4);

      if (val_vec.size() == 1) // normal case.
        set_Jcoupling_scalar(id1, id2, boost::lexical_cast<double>(val), unit);
      else if (val_vec.size() == 3) { // broadband case.
        vec3 par;
        par[0] = boost::lexical_cast<double>(val_vec[0]);
        par[1] = boost::lexical_cast<double>(val_vec[1]);
        par[2] = boost::lexical_cast<double>(val_vec[2]);
        //set_Jcoupling_scalar(id1, id2, par[0], unit); // center j-coup.
        set_Jcoupling_scalar(id1, id2, (par[0]+par[1])/2, unit); // center j-coup.
        set_Jcoupling_scalar_broadband(id1, id2, par, unit);
      }
      iter += 5;
    } else if (*(iter + 1) == "eigs-euler") {
      // TODO.
    } else {
      std::string s = "unidentified zeeman parameters: " + list;
      throw std::runtime_error(s.c_str());
    }
  }
#ifdef SSL_OUTPUT_ENABLE
  std::string s = str(boost::format("%s %s.\n") % "J coupling pars set to be" % list);
  ssl_color_text("info", s);
#endif
}
void interaction::set_zeeman(std::vector<std::string> list) {
  if (!list.size())
    return;
  std::string par_zeeman = boost::algorithm::join(list, " ");
  set_zeeman(par_zeeman);
}
void interaction::set_zeeman_scalar(size_t id, double val, std::string unit) {
  if (id > comp_.nspins()) {
    std::string s = "zeeman scalar id out of range: " + boost::lexical_cast<std::string>(id);
    throw std::runtime_error(s.c_str());
  }
  if (unit == "ppm") { // do nothing.
  } else if (unit == "hz")
    val /= comp_.get_freq(id - 1); // Hz into ppm.
  else if (unit == "khz") {
    val *= 1e3;
    val /= comp_.get_freq(id - 1); // kHz into ppm.
  } else {
    std::string s = "unidentified zeeman scalar unit: " + unit;
    throw std::runtime_error(s.c_str());
  }
  zeeman_.scalars[id - 1] = val;
}

void interaction::set_zeeman_scalar_broadband(size_t id, vec3 val, std::string unit) {
  if (id > comp_.nspins()) {
    std::string s = "zeeman scalar id out of range: " + boost::lexical_cast<std::string>(id);
    throw std::runtime_error(s.c_str());
  }

  int num = int(val[2]);
  if (num < 2) {
    std::string s = "broadband zeeman scalar set failed, spin No. " + boost::lexical_cast<std::string>(id);
    throw std::runtime_error(s.c_str());
  }

  broadband_cs bb;
  //bb.nominal_offset = vec::LinSpaced(num, val[0] - val[1] / 2, val[0] + val[1] / 2);
  if (unit == "khz") {
    val[0] *= 1e3;
    val[1] *= 1e3;
  }

  bb.nominal_offset = vec::LinSpaced(num, val[0], val[1]);

  if (unit == "ppm") {
  } else if (unit == "hz") {
    val[0] /= comp_.get_freq(id - 1); // Hz into ppm.
    val[1] /= comp_.get_freq(id - 1); // Hz into ppm.
  } else if (unit == "khz") {
    //val[0] *= 1e3;
    //val[1] *= 1e3;
    val[0] /= comp_.get_freq(id - 1); // kHz into ppm.
    val[1] /= comp_.get_freq(id - 1); // kHz into ppm.
  } else {
    std::string s = "unidentified zeeman scalar unit: " + unit;
    throw std::runtime_error(s.c_str());
  }


  bb.spin_id = id - 1;
  //bb.offset = vec::LinSpaced(num, -val[1] / 2, val[1] / 2);
  bb.offset = vec::LinSpaced(num, -(val[1]-val[0]) / 2, (val[1]-val[0]) / 2);

  zeeman_.bb_scalars.push_back(bb);
}

void interaction::set_Jcoupling(std::vector<std::string> list) {
  if (!list.size())
    return;
  std::string par_coupling = boost::algorithm::join(list, " "); // "1 2 scalar 14 Hz"
  set_Jcoupling(par_coupling);
}
void interaction::set_Jcoupling_coords(const sol::table &t) {
  for (auto &kv : t) {
    sol::object val = kv.second;
    const sol::table &nt = val.as<sol::table>();
    vec3 pos;
    int i = 0;
    for (auto &kvv : nt) {
      pos[i++] = kvv.second.as<double>();
    }
    coupling_.coordinates.push_back(pos);
    //std::cout<<pos<<"\n";
    //std::cout<<pos.norm()<<"\n";
  }
}

void interaction::set_relaxation(std::vector<std::string> list) {
  if (!list.size())
    return;
  std::string par_relax = boost::algorithm::join(list, " ");
  set_relaxation(par_relax);
}
void interaction::set_relaxation(std::string list) {
  boost::trim(list); // trim the spaces on both left and right sides.
  boost::to_lower(list);

  std::vector<std::string> par_vec;
  boost::split(par_vec, list, boost::is_any_of("\t, "), boost::token_compress_on);

  std::vector<std::string>::iterator t1_pos = find(par_vec.begin(),
                                         par_vec.end(),
                                         "t1");
  std::vector<std::string>::iterator t2_pos = find(par_vec.begin(),
                                         par_vec.end(),
                                         "t2");

  size_t nspins = comp_.nspins();
  if (t1_pos != par_vec.end()) {
    double unit = 1;// s.
    if (*(t1_pos + nspins + 1) == "ms")
      unit *= 1e-3;
    for (size_t i = 1; i <= nspins; i++)
      relax_.r1.push_back(1.0 / (boost::lexical_cast<double>(*(t1_pos + i)) * unit));
  }
  if (t2_pos != par_vec.end()) {
    double unit = 1;// s.
    if (*(t2_pos + nspins + 1) == "ms")
      unit *= 1e-3;
    for (size_t i = 1; i <= nspins; i++)
      relax_.r2.push_back(1.0 / (boost::lexical_cast<double>(*(t2_pos + i)) * unit));
  }

#ifdef SSL_OUTPUT_ENABLE
  std::string s = str(boost::format("%s %s.\n") % "relaxation pars set to be" % list);
  ssl_color_text("info", s);
#endif
}

void interaction::set_Jcoupling_scalar(size_t id1, size_t id2, double val,
                                       std::string unit) {
  if (unit == "hz") {
    coupling_.scalar(id1 - 1, id2 - 1) = val;
    if (id1 != id2)
      coupling_.scalar(id2 - 1, id1 - 1) = val;
  } else {
    std::string s = "unidentified Jcoupling scalar unit: " + unit;
    throw std::runtime_error(s.c_str());
  }
}

void interaction::set_Jcoupling_scalar_broadband(size_t id1, size_t id2, vec3 val, std::string unit) {
  if (unit == "hz") {
    int num = int(val[2]);
    if (num < 2) {
      std::string s = "broadband Jcoupling scalar set failed, spin No. " + boost::lexical_cast<std::string>(id1) + " and "
          + boost::lexical_cast<std::string>(id2);
      throw std::runtime_error(s.c_str());
    }

    broadband_jcoup bb;
    //bb.nominal_offset = vec::LinSpaced(num, val[0] - val[1] / 2, val[0] + val[1] / 2);
    bb.nominal_offset = vec::LinSpaced(num, val[0], val[1]);
    bb.spin_id1 = id1 - 1;
    bb.spin_id2 = id2 - 1;
    bb.offset = vec::LinSpaced(num, -(val[1]-val[0]) / 2, (val[1]-val[0]) / 2);
    coupling_.bb_scalars.push_back(bb);
  } else {
    std::string s = "unidentified Jcoupling scalar unit: " + unit;
    throw std::runtime_error(s.c_str());
  }
}
std::vector<std::vector<cs_par>> interaction::parsing_zeeman_broadband() {
  // case 1: broadband spin system. ONLY for the specific spin.
  std::vector<std::vector<cs_par>> result0;
  if (zeeman_.bb_scalars.size() == 0)
      return result0;
  for (int i = 0; i < zeeman_.bb_scalars[0].offset.size(); i++) {
      cs_par tmp;
      tmp.id = zeeman_.bb_scalars[0].spin_id;
      tmp.val = zeeman_.bb_scalars[0].offset[i];
      std::vector<cs_par> b;
      b.push_back(tmp);
      result0.push_back(b);
  }
  std::vector<std::vector<cs_par>> result = result0;
  for (int i = 1; i < zeeman_.bb_scalars.size(); i++) {
      result = parsing_zeeman_broadband(result, zeeman_.bb_scalars[i]);
  }
  return result;


  // case 2: broadband magnet field. NOTE: chemical shifts of all spins are shifted rigidly.
  //std::vector<std::vector<cs_par>> result;
  //if (zeeman_.bb_scalars.size() == 0)
  //  return result;

  //vec offset = zeeman_.bb_scalars[0].offset;
  //for (int i = 0; i < offset.size(); i++) {
  //  std::vector<cs_par> b;
  //  for (size_t j = 0; j < comp_.nspins(); j++) {
  //    cs_par tmp;
  //    tmp.id = j;
  //    tmp.val = offset[i];
  //    b.push_back(tmp);
  //  }
  //  result.push_back(b);
  //}
  //return result;
}
std::vector<std::vector<jcoup_par>> interaction::parsing_Jcoupling_broadband() {
  std::vector<std::vector<jcoup_par>> result;
  if (coupling_.bb_scalars.size() == 0)
    return result;

  std::vector<std::vector<jcoup_par>> result0;
  for (int i = 0; i < coupling_.bb_scalars[0].offset.size(); i++) {
    jcoup_par tmp;
    tmp.id1 = coupling_.bb_scalars[0].spin_id1;
    tmp.id2 = coupling_.bb_scalars[0].spin_id2;
    tmp.val = coupling_.bb_scalars[0].offset[i];
    std::vector<jcoup_par> b;
    b.push_back(tmp);
    result0.push_back(b);
  }
  result = result0;
  for (size_t i = 1; i < coupling_.bb_scalars.size(); i++) {
    result = parsing_Jcoupling_broadband(result, coupling_.bb_scalars[i]);
  }
  return result;
}

std::vector<std::vector<cs_par>> interaction::parsing_zeeman_broadband(const std::vector<std::vector<cs_par>> &old,
                                                             const broadband_cs &bb) {
  std::vector<std::vector<cs_par>> result;
  cs_par cur;
  cur.id = bb.spin_id;
  for (int i = 0; i < bb.offset.size(); i++) {
    cur.val = bb.offset[i];
    for (size_t j = 0; j < old.size(); j++) {
      std::vector<cs_par> tmp = old[j];
      tmp.push_back(cur);
      result.push_back(tmp);
    }
  }
  return result;
}

std::vector<std::vector<jcoup_par>> interaction::parsing_Jcoupling_broadband(const std::vector<std::vector<jcoup_par>> &old,
                                                                   const broadband_jcoup &bb) {
  std::vector<std::vector<jcoup_par>> result;
  jcoup_par cur;
  cur.id1 = bb.spin_id1;
  cur.id2 = bb.spin_id2;
  for (int i = 0; i < bb.offset.size(); i++) {
    cur.val = bb.offset[i];
    for (size_t j = 0; j < old.size(); j++) {
      std::vector<jcoup_par> tmp = old[j];
      tmp.push_back(cur);
      result.push_back(tmp);
    }
  }
  return result;
}

void interaction::alloc() {
  size_t nspins = comp_.nspins();
  zeeman_ = zeeman(nspins);
  coupling_ = coupling(nspins);
}
void interaction::set_assumption(std::string assume) {
  if (assume_ == assume)
    return;
  assume_ = assume;
  if (assume_ != "labframe" && assume_ != "nmr" && assume_ != "esr") {
    std::string s = "unknown interaction assumption ** " + assume + " **!";
    throw std::runtime_error(s.c_str());
  }
  size_t nspins = comp_.nspins();
  if (assume == "labframe") {
    // Process zeeman interactions
    for (size_t i = 0; i < nspins; i++)
      zeeman_.strengths[i] = full;

    // Process couplings
    for (size_t i = 0; i < nspins; i++)
      for (size_t j = i; j < nspins; j++)
        coupling_.strengths[i * nspins + j] = strong;  // full coupling tensors should be used for all spins
  }
  if (assume == "nmr" || assume == "esr") {
    // secular shielding terms for nuclear zeeman interactions.
    // secular coupling terms for spins belonging to the same species.
    // weak coupling terms for spins belonging to different species.
    // Process zeeman interactions
    for (size_t i = 0; i < nspins; i++)
      zeeman_.strengths[i] = secular;

    // Process couplings
    for (size_t i = 0; i < nspins; i++)
      for (size_t j = i; j < nspins; j++) {
        if (comp_.spins()[i] == comp_.spins()[j])
          coupling_.strengths[i * nspins + j] = secular;
        else
          coupling_.strengths[i * nspins + j] = weak;
      }
  }
}
void interaction::init() {
  size_t nspins = comp_.nspins();
  // -----------------Preallocate zeeman tensor array-----------------
  mat33 unit = mat33::Identity();
  for (size_t i = 0; i < nspins; i++) {
    // SCALARS case
    if (zeeman_.scalars[i])
      zeeman_.matrices[i] += zeeman_.scalars[i] * unit;
    // EIGEN_EULERS case
    if (zeeman_.eigs[i] != vec3::Zero()) {
      mat33 m = unit;
      if (zeeman_.eulers[i] != vec3::Zero())
        m = euler2dcm(zeeman_.eulers[i]);
      zeeman_.matrices[i] += m.transpose() * zeeman_.eigs[i].asDiagonal() * m;  // TO BE CAREFULL!
    }
  }
  // convert to angular frequencies
  for (size_t i = 0; i < nspins; i++)
    // For electrons, assume that the g-factor is given and compute the offset from the free electron g-factor
    if (comp_.spins()[i].symbol() == "E") {
      double freeg = 2.00231930436220;
      zeeman_.matrices[i] -= unit * freeg;
      zeeman_.matrices[i] *= (comp_.base_freq_[i] / freeg);
    } else
      // For nuclei, assume that the chemical shift is given and compute the corresponding offset
      zeeman_.matrices[i] *= -1e-6 * comp_.base_freq_[i];

  // -----------------Preallocate the coupling tensor array-----------------
  // Process the coordinates
  if (coupling_.coordinates.size()) {
    vec3 orient;
    mat33 mat;
    for (size_t i = 0; i < nspins; i++)
      for (size_t j = i + 1; j < nspins;
           j++) {  // single counting (1,2) (1 3) ....
        orient = coupling_.coordinates[i] - coupling_.coordinates[j];
        double r = orient.norm();
        //if (r < tol_.prox_cutoff)  // Dipolar interaction for a pair of spins at a separation
        // (in Angstron) greater than the threshold is set to 0.
        {
          //std::cout << r << "###########\n";
          //std::cout << orient << "\n";
          double x = orient(_cx) / r;
          double y = orient(_cy) / r;
          double z = orient(_cz) / r;
          mat(0, 0) = 1 - 3 * x * x;
          mat(0, 1) = -3 * x * y;
          mat(0, 2) = -3 * x * z;

          mat(1, 0) = -3 * y * x;
          mat(1, 1) = 1 - 3 * y * y;
          mat(1, 2) = -3 * y * z;

          mat(2, 0) = -3 * z * x;
          mat(2, 1) = -3 * z * y;
          mat(2, 2) = 1 - 3 * z * z;

          //std::cout << mat << "##\n";

          r *= 1e-10;  // into A
          double coeff = kMu * kHbar * comp_.gamma(i) * comp_.gamma(j) / (r * r * r); // constant_dipolar.

          coupling_.matrices[i * nspins + j] = coeff * mat;
          //std::cout << coupling_.matrices[i * nspins + j] << "\n";
        }
      }
  }
  // Absorb the user-specified couplings
  //if (coupling_.scalar.norm() > tol_inter_cutoff) {
  //mat33 unit = mat33::Identity();
  for (size_t i = 0; i < nspins; i++)
    for (size_t j = i; j < nspins; j++) {
      // SCALARS case
      if (coupling_.scalar(i, j)) {
        coupling_.matrices[i * nspins + j] += 2 * _pi * coupling_.scalar(i, j) * unit;
        //std::cout << coupling_.matrices[i * nspins + j] << "\n";
      }
      // EIGEN_EULERS case
      if (coupling_.eigs[i] != vec3::Zero()) {
        mat33 m = unit;
        if (coupling_.eulers[i] != vec3::Zero())
          m = euler2dcm(coupling_.eulers[i * nspins + j]);
        coupling_.matrices[i * nspins + j] +=
            2 * _pi * m.transpose() * coupling_.eigs[i * nspins + j].asDiagonal() * m;  // TO BE CAREFULL!
      }
      if (i != j)
        coupling_.matrices[j * nspins + i] = coupling_.matrices[i * nspins + j];
    }
  //}
  /* for (size_t k = 0; k < coupling_.cm.size(); k++) {
           size_t i = coupling_.cm[k].i;
           size_t j = coupling_.cm[k].j;
           coupling_.matrix[i * nspins + j] += 2 * _pi* coupling_.cm[k].m;
   }*/
  cacl_equilibrium_state(); //MUST BE CALLED PRIOR TO 'cacl_relaxation'.
  cacl_relaxation();
}

void interaction::init_broadband() {
  // firstly, reset interaction matrix.
  for (size_t i = 0; i < zeeman_.matrices.size(); i++)
    zeeman_.matrices[i].setZero();

  for (size_t i = 0; i < coupling_.matrices.size(); i++)
    coupling_.matrices[i].setZero();

  size_t nspins = comp_.nspins();
  // -----------------Preallocate zeeman tensor array-----------------
  mat33 unit = mat33::Identity();
  for (size_t i = 0; i < nspins; i++) {
    // SCALARS case
    if (zeeman_.scalars[i])
      zeeman_.matrices[i] += zeeman_.scalars[i] * unit;
    // EIGEN_EULERS case
    //if (zeeman_.eigs[i] != vec3::Zero()) {
    //	mat33 m = unit;
    //	if (zeeman_.eulers[i] != vec3::Zero())
    //		m = euler2dcm(zeeman_.eulers[i]);
    //	zeeman_.matrices[i] += m.transpose() * zeeman_.eigs[i].asDiagonal() * m;  // TO BE CAREFULL!
    //}
  }
  // convert to angular frequencies
  for (size_t i = 0; i < nspins; i++)
    // For electrons, assume that the g-factor is given and compute the offset from the free electron g-factor
    if (comp_.spins()[i].symbol() == "E") {
      double freeg = 2.00231930436220;
      zeeman_.matrices[i] -= unit * freeg;
      zeeman_.matrices[i] *= (comp_.base_freq_[i] / freeg);
    } else
      // For nuclei, assume that the chemical shift is given and compute the corresponding offset
      zeeman_.matrices[i] *= -1e-6 * comp_.base_freq_[i];

  for (size_t i = 0; i < nspins; i++)
    for (size_t j = i; j < nspins; j++) {
      // SCALARS case
      if (coupling_.scalar(i, j))
        coupling_.matrices[i * nspins + j] += 2 * _pi * coupling_.scalar(i, j) * unit;
      // EIGEN_EULERS case
      //if (coupling_.eigs[i] != vec3::Zero()) {
      //	mat33 m = unit;
      //	if (coupling_.eulers[i] != vec3::Zero())
      //		m = euler2dcm(coupling_.eulers[i * nspins + j]);
      //	coupling_.matrices[i * nspins + j] += 2 * _pi * m.transpose()* coupling_.eigs[i * nspins + j].asDiagonal() *m;  // TO BE CAREFULL!
      //}
      if (i != j)
        coupling_.matrices[j * nspins + i] = coupling_.matrices[i * nspins + j];
    }
}

mat33 euler2dcm(const vec3 &euler) {
  double alpha = euler(0);
  double beta = euler(1);
  double gamma = euler(2);
  mat33 m = mat33::Zero();
  m = Eigen::AngleAxisd(alpha, vec3::UnitZ()) * Eigen::AngleAxisd(beta,
                                                                  vec3::UnitY())
      * Eigen::AngleAxisd(gamma, vec3::UnitZ());
  return m;
}
ham_op interaction::hamiltonian(const op_side type, bool build_aniso) const {
  /* if (assume_ == kNoAssumption) {
           std::cout
                           << format("%1% %2%.\n") % "S-S-L warning: "
                           % "simulation assumption should be specified heretofore";
           exit(0);
   }*/
  size_t nspins = comp_.nspins();
  std::vector<op_item> iso_ham;
  std::vector<op_item> *ani_ham = NULL;
  if (build_aniso)
    ani_ham = new std::vector<op_item>[25];
  // Process zeeman interactions
  for (size_t i = 0; i < nspins; i++) {
    // Write the isotropic part.
    double zeeman_iso;
    switch (zeeman_.strengths[i]) {
      case full: {
        //Add the carrier frequency
        zeeman_iso = zeeman_.matrices[i].trace() / 3
            + comp_.base_freq_[i];
        if (fabs(zeeman_iso) > tol_inter_cutoff) {
#ifdef SSL_ENABLE_DEBUG_SPIN_SYS
          std::cout
                  << format("%30T %=10s %-20s %-20d \n") % "---->"
                  % "complete isotropic zeeman interaction for spin " % i;
          std::cout
                  << format("%31T %40s %20.4e Hz\n") % "(Lz) x  "
                  % (zeeman_iso / (2 * k_pi));
#endif
          iso_ham.push_back(op_item(zeeman_iso, "Iz", "", i, -1));
        }
      }
        break;
      case secular: {
        //Skip the carrier frequency
        zeeman_iso = zeeman_.matrices[i].trace() / 3;
        if (fabs(zeeman_iso) > tol_inter_cutoff) {
#ifdef SSL_ENABLE_DEBUG_SPIN_SYS
          std::cout
                  << format("%30T %=10s %-20s %-20d \n") % "---->"
                  % "offset isotropic zeeman interaction for spin " % i;
          std::cout
                  << format("%31T %40s %20.4e Hz\n") % "(Lz) x  "
                  % (zeeman_iso / (2 * k_pi));
#endif
          iso_ham.push_back(op_item(zeeman_iso, "Iz", "", i, -1));
        }
      }
        break;
      default:break;
    }

    // Process anisotropic part if required.
    if (build_aniso) {
      rank T = mat2sphten(zeeman_.matrices[i]);
      if (T.r2.lpNorm<1>() > tol_inter_cutoff) {
        // Write irreducible spherical tensors.
        op_item *tt = new op_item[5];
        switch (zeeman_.strengths[i]) {
          case full: {
#ifdef SSL_ENABLE_DEBUG_SPIN_SYS
            std::cout
                    << format("%30T %=10s %-20s %-20d \n") % "---->"
                    % "complete anisotropic zeeman interaction for spin " % i;
            std::cout
                    << format("%31T %40s %20s Hz\n") % "-0.5*(Lp) x "
                    % (T.r2(1) / (2 * k_pi));
            std::cout
                    << format("%31T %40s %20s Hz\n") % "sqrt(2/3)*(Lz) x "
                    % (T.r2(2) / (2 * k_pi));
            std::cout
                    << format("%31T %40s %20s Hz\n") % "0.5*(Lm) x "
                    % (T.r2(3) / (2 * k_pi));
#endif
            tt[0] = op_item();
            tt[1] = op_item(-0.5, "I+", "", i, -1);
            tt[2] = op_item(sqrt(2.0 / 3), "Iz", "", i, -1);
            tt[3] = op_item(0.5, "I-", "", i, -1);
            tt[4] = op_item();
          }
            break;
          case secular: {
#ifdef SSL_ENABLE_DEBUG_SPIN_SYS
            std::cout
                    << format("%30T %=10s %-20s %-20d \n") % "---->"
                    % "Z part of the anisotropic zeeman interaction for spin "
                    % i;
            std::cout
                    << format("%31T %40s %20s Hz\n") % "sqrt(2/3)*(Lz) x "
                    % (T.r2(2) / (2 * k_pi));
#endif
            tt[0] = op_item();
            tt[1] = op_item();
            tt[2] = op_item(sqrt(2.0 / 3), "Iz", "", i, -1);
            tt[3] = op_item();
            tt[4] = op_item();
          }
            break;
          default:break;
        }
        for (size_t k = 0; k < 5; k++)
          if (tt[k].used)
            for (size_t m = 0; m < 5; m++) {
              op_item tmp = tt[k];
              tmp.coeff *= T.r2(m);
              std::vector<op_item> Qmn;
              Qmn.push_back(tmp);
              copy(Qmn.begin(), Qmn.end(), back_inserter(ani_ham[k * 5 + m]));
            }
      }

      size_t index = i * nspins + i;// (2 * comp_.nspins + 1 - i) * i / 2;
      T = mat2sphten(coupling_.matrices[index]);
      if (T.r2.cwiseAbs().sum() > tol_inter_cutoff) {
        std::vector<op_item> *tt = new std::vector<op_item>[5];
        switch (coupling_.strengths[index]) {
          case strong: {
#ifdef SSL_ENABLE_DEBUG_SPIN_SYS
            std::cout
                    << format("%30T %=10s %-20s %=3d \n") % "---->"
                    % "complete quadratic coupling for spins " % i;
            std::cout
                    << format("%31T %40s %20s Hz\n") % "T(2,+2) x "
                    % (T.r2(0) / (2 * k_pi));
            std::cout
                    << format("%31T %40s %20s Hz\n") % "T(2,+1) x "
                    % (T.r2(1) / (2 * k_pi));
            std::cout
                    << format("%31T %40s %20s Hz\n") % "T(2,0) x "
                    % (T.r2(2) / (2 * k_pi));
            std::cout
                    << format("%31T %40s %20s Hz\n") % "T(2,-1) x "
                    % (T.r2(3) / (2 * k_pi));
            std::cout
                    << format("%31T %40s %20s Hz\n") % "T(2,-2) x "
                    % (T.r2(4) / (2 * k_pi));
#endif

            std::vector<op_item> tmp;
            tmp.push_back(op_item(1.0, "T(2,2)", "", i, -1));
            tt[0] = tmp;
            tmp.clear();

            tmp.push_back(op_item(1.0, "T(2,1)", "", i, -1));
            tt[1] = tmp;
            tmp.clear();

            tmp.push_back(op_item(1.0, "T(2,0)", "", i, -1));
            tt[2] = tmp;
            tmp.clear();

            tmp.push_back(op_item(1.0, "T(2,-1)", "", i, -1));
            tt[3] = tmp;
            tmp.clear();

            tmp.push_back(op_item(1.0, "T(2,-2)", "", i, -1));
            tt[4] = tmp;
            tmp.clear();
          }
            break;
          case secular: {
#ifdef SSL_ENABLE_DEBUG_SPIN_SYS
            std::cout
                    << format("%30T %=10s %-20s %=3d \n") % "---->"
                    % "secular part of the quadratic coupling_ for spin " % i;
            std::cout
                    << format("%31T %40s %20s Hz\n") % "T(2,0) x "
                    % (T.r2(2) / (2 * k_pi));
#endif
            std::vector<op_item> tmp;
            tt[0] = tmp;
            tt[1] = tmp;

            tmp.push_back(op_item(1.0, "T(2,0)", "", i, -1));
            tt[2] = tmp;
            tmp.clear();

            tt[3] = tmp;
            tt[4] = tmp;
          }
            break;
          default:break;
        }
        for (size_t k = 0; k < 5; k++)
          for (size_t m = 0; m < 5; m++) {
            std::vector<op_item> Qmn;
            for (size_t p = 0; p < tt[k].size(); p++) {
              op_item tmp = tt[k][p];
              tmp.coeff *= T.r2(m);
              Qmn.push_back(tmp);
            }
            copy(Qmn.begin(), Qmn.end(), back_inserter(ani_ham[k * 5 + m]));
          }
      }

    }

  }

  // Process coupling tensors
  for (size_t i = 0; i < nspins; i++)
    for (size_t j = i + 1; j < nspins; j++) {
      size_t index = i * nspins + j;// (2 * comp_.nspins + 1 - i) * i / 2 + j - i;
      // Get the isotropic coupling constant
      double coupling_iso = coupling_.matrices[index].trace() / 3;
      if (fabs(coupling_iso) > tol_inter_cutoff)
        switch (coupling_.strengths[index]) {
          case strong:
          case secular: {
#ifdef SSL_ENABLE_DEBUG_SPIN_SYS
            std::cout
                    << format("%30T %=10s %-20s %=3d  %=3d\n") % "---->"
                    % "complete isotropic coupling for spins " % i % j;
            std::cout
                    << format("%31T %40s %20s Hz\n") % "(LxSx+LySy+LzSz) x "
                    % (coupling_iso / (2 * k_pi));
#endif
            iso_ham.push_back(op_item(coupling_iso, "Iz", "Iz", i, j));
            iso_ham.push_back(op_item(0.5 * coupling_iso, "I+", "I-", i, j));
            iso_ham.push_back(op_item(0.5 * coupling_iso, "I-", "I+", i, j));
          }
            break;
          case weak: {
#ifdef SSL_ENABLE_DEBUG_SPIN_SYS
            std::cout
                    << format("%30T %=10s %-20s %=3d  %=3d\n") % "---->"
                    % "(z,z) part of the isotropic coupling_ for spins " % i
                    % j;
            std::cout
                    << format("%31T %40s %20s Hz\n") % "(LzSz) x "
                    % (coupling_iso / (2 * k_pi));
#endif
            iso_ham.push_back(op_item(coupling_iso, "Iz", "Iz", i, j));
          }
            break;
          default:break;
        }

      //Process anisotropic part if required
      if (build_aniso) {
        rank T = mat2sphten(coupling_.matrices[index]);
        if (T.r2.cwiseAbs().sum() > tol_inter_cutoff) {
          std::vector<op_item> *tt = new std::vector<op_item>[5];
          switch (coupling_.strengths[index]) {
            case strong: {
#ifdef SSL_ENABLE_DEBUG_SPIN_SYS
              std::cout
                      << format("%30T %=10s %-20s %=3d  %=3d\n") % "---->"
                      % "complete anisotropic coupling for spins " % i % j;

              std::cout
                      << format("%31T %40s %20s Hz\n") % "0.5*(LpSp) x "
                      % (T.r2(0) / (2 * k_pi));
              std::cout
                      << format("%31T %40s %20s Hz\n") % "-0.5*(LzSp+LpSz) x "
                      % (T.r2(1) / (2 * k_pi));
              std::cout
                      << format("%31T %40s %20s Hz\n")
                      % "sqrt(2/3)*(LzSz-0.25*(LpSm+LmSp)) x "
                      % (T.r2(2) / (2 * k_pi));
              std::cout
                      << format("%31T %40s %20s Hz\n") % "0.5*(LzSm+LmSz) x "
                      % (T.r2(3) / (2 * k_pi));
              std::cout
                      << format("%31T %40s %20s Hz\n") % "0.5*(LmSm) x "
                      % (T.r2(4) / (2 * k_pi));
#endif
              // Eq(3) in Kuprov 2011
              std::vector<op_item> tmp;
              tmp.push_back(op_item(0.5, "I+", "I+", i, j));
              tt[0] = tmp;
              tmp.clear();

              tmp.push_back(op_item(-0.5, "Iz", "I+", i, j));
              tmp.push_back(op_item(-0.5, "I+", "Iz", i, j));
              tt[1] = tmp;
              tmp.clear();

              tmp.push_back(op_item(1.00 * sqrt(2.0 / 3), "Iz", "Iz", i, j));
              tmp.push_back(op_item(-0.25 * sqrt(2.0 / 3), "I+", "I-", i, j));
              tmp.push_back(op_item(-0.25 * sqrt(2.0 / 3), "I-", "I+", i, j));
              tt[2] = tmp;
              tmp.clear();

              tmp.push_back(op_item(0.5, "Iz", "I-", i, j));
              tmp.push_back(op_item(0.5, "I-", "Iz", i, j));
              tt[3] = tmp;
              tmp.clear();

              tmp.push_back(op_item(0.5, "I-", "I-", i, j));
              tt[4] = tmp;
              tmp.clear();
            }
              break;
            case secular: {
#ifdef SSL_ENABLE_DEBUG_SPIN_SYS
              std::cout
                      << format("%30T %=10s %-20s %=3d  %=3d\n") % "---->"
                      % "secular part of the anisotropic coupling for spins "
                      % i % j;

              std::cout
                      << format("%31T %40s %20s Hz\n")
                      % "sqrt(2/3)*(LzSz-0.25*(LpSm+LmSp)) x "
                      % (T.r2(2) / (2 * k_pi));
#endif
              // Eq(3) in Kuprov 2011
              std::vector<op_item> tmp;
              tt[0] = tmp;
              tt[1] = tmp;

              tmp.push_back(op_item(1.00 * sqrt(2.0 / 3), "Iz", "Iz", i, j));
              tmp.push_back(op_item(-0.25 * sqrt(2.0 / 3), "I+", "I-", i, j));
              tmp.push_back(op_item(-0.25 * sqrt(2.0 / 3), "I-", "I+", i, j));
              tt[2] = tmp;
              tmp.clear();

              tt[3] = tmp;
              tt[4] = tmp;
            }
              break;
            case weak: {
#ifdef SSL_ENABLE_DEBUG_SPIN_SYS
              std::cout
                      << format("%30T %=10s %-20s %=3d  %=3d\n") % "---->"
                      % "(z,z) part of the anisotropic coupling for spins " % i
                      % j;

              std::cout
                      << format("%31T %40s %20s Hz\n") % "sqrt(2/3)*LzSz x "
                      % (T.r2(2) / (2 * k_pi));
#endif
              // Eq(3) in Kuprov 2011
              std::vector<op_item> tmp;
              tt[0] = tmp;
              tt[1] = tmp;

              tmp.push_back(op_item(1.00 * sqrt(2.0 / 3), "Iz", "Iz", i, j));
              tt[2] = tmp;
              tmp.clear();

              tt[3] = tmp;
              tt[4] = tmp;
            }
              break;
            default:break;
          }
          for (size_t k = 0; k < 5; k++)
            for (size_t m = 0; m < 5; m++) {
              std::vector<op_item> Qmn;
              for (size_t p = 0; p < tt[k].size(); p++) {
                op_item tmp = tt[k][p];
                tmp.coeff *= T.r2(m);
                Qmn.push_back(tmp);
              }
              copy(Qmn.begin(), Qmn.end(), back_inserter(ani_ham[k * 5 + m]));
            }
        }
      }
    }
#ifdef SSL_ENABLE_DEBUG_SPIN_SYS
  std::cout
          << format("\n%=20s %=10s %-20s %3s items.\n") % "Hamiltonian" % "---->"
          % "Build the isotropic Hamiltonian: " % iso_ham.size();
#endif

  // H: isotropic part
  ham_op ham;

  ham.isotropic = aux_hamiltonian(iso_ham, type);
  // Q: twenty-five matrices giving the irreducible components of the anisotropic part
#ifdef SSL_ENABLE_DEBUG_SPIN_SYS
  std::cout
          << format("\n%=20s %=10s %-20s \n") % "Hamiltonian" % "---->"
          % "Build the rotational basis.";
  std::cout
          << format("%30T %=10s %-20s\n") % "---->"
          % "building irreducible spherical component.";
#endif
  if (build_aniso) {
    ham.anisotropic = new sp_cx_mat[25];
    for (int k = 0; k < 5; k++)
      for (int m = 0; m < 5; m++) {
#ifdef SSL_ENABLE_DEBUG_SPIN_SYS
        std::cout
                << format("%41T k= %2s %3T  m= %2s %3T %2s items\n") % (2 - k)
                % (2 - m) % ani_ham[k * 5 + m].size();
#endif
        ham.anisotropic[k * 5 + m] = aux_hamiltonian(ani_ham[k * 5 + m], type);
      }
  }
  return ham;
}
sp_cx_mat interaction::aux_hamiltonian(const std::vector<op_item> &descr,
                                       const op_side type) const {
  size_t dim = comp_.liouville_space_dim();
  sp_cx_mat H = sp_cx_mat(dim, dim);

  for (size_t i = 0; i < descr.size(); i++) {
    std::string par = "";
    switch (descr[i].used) {  // tensor type
      case 0:  // not used
        break;
      case 1: {
        par += boost::lexical_cast<std::string>(descr[i].id1 + 1);
        par += " " + descr[i].op1;
      }
        break;
      case 2: {
        par += boost::lexical_cast<std::string>(descr[i].id1 + 1);
        par += " " + descr[i].op1;
        par += " " + boost::lexical_cast<std::string>(descr[i].id2 + 1);
        par += " " + descr[i].op2;
      }
        break;
      default:break;
    }
    if (par != "" && descr[i].coeff != cd(0, 0)) {
      sp_mat real = comp_.op(par, type);
      //sp_mat imag = sp_mat(real.rows, real.cols);
      //H += sp_cx_mat(real, imag) * descr[i].coeff;
      H += real * descr[i].coeff;
    }
  }
  return H;
}
rank mat2sphten(const mat33 &m) {
  rank rank;
  // rank 0 component
  rank.r0 = -(1 / sqrt(3)) * m.trace();

  // rank 1 components
  cd a, b, c, d, e;
  a = -0.5 * (m(2, 0) - m(0, 2) - ci * (m(2, 1) - m(1, 2)));
  b = -(ci / sqrt(2)) * (m(0, 1) - m(1, 0));
  c = -0.5 * (m(2, 0) - m(0, 2) + ci * (m(2, 1) - m(1, 2)));

  rank.r1 << a, b, c;

  // rank 2 components
  a = +(0.5) * (m(0, 0) - m(1, 1) - ci * (m(0, 1) + m(1, 0)));
  b = -(0.5) * (m(0, 2) + m(2, 0) - ci * (m(1, 2) + m(2, 1)));
  c = (1 / sqrt(6)) * (2 * m(2, 2) - m(0, 0) - m(1, 1));
  d = +(0.5) * (m(0, 2) + m(2, 0) + ci * (m(1, 2) + m(2, 1)));
  e = +(0.5) * (m(0, 0) - m(1, 1) + ci * (m(0, 1) + m(1, 0)));

  rank.r2 << a, b, c, d, e;
  return rank;
}
void interaction::cacl_relaxation() {
  sp_cx_mat mat_R = comp_.mprealloc();
  switch (relax_.theory) {
    case redfield: {
      std::cout
          << boost::format("%=20s %=10s %-20s \n") % "Relaxation" % "---->"
              % "statrt to calculate redfield superoperator.";
      // If correlation times are zero, just return a zero matrix
      if (relax_.tau.size() == 0) {
        relax_.R = mat_R;
        return;
      }


      // Determine the largest and the smallest correlation time
      double tau_max = *max_element(relax_.tau.begin(), relax_.tau.end());
      double tau_min = *min_element(relax_.tau.begin(), relax_.tau.end());

      // Get the rotational basis, including the non-secular terms
      set_assumption("labframe");

      sp_cx_mat L0;
      sp_cx_mat *Q;
      ham_op ham = hamiltonian(kComm, true);
      L0 = ham.isotropic;
      Q = ham.anisotropic;

      // Determine the norm of L0
      double L0_norm = norm(L0);
      double xi = tol_rlx_integration;
      // Set the upper integration limit according to the accuracy goal
      double t_max = tau_max * log(1 / xi);
      // Set the number of integration steps according to the accuracy goal
      double nsteps = ceil(
          1 / (10 * std::pow(xi, 1.0 / 6))
              * std::pow(t_max / std::min(tau_max, 1 / L0_norm), 7.0 / 6));
      // Get the numerical integration step
      double timestep = t_max / nsteps;
      // Kill the terms in the static Liouvillian that are irrelevant on the time scale of the integration
      cleanup(L0, 1 / t_max);
      std::cout
          << boost::format("%=20s %=10s %-47s %-10s \n") % "Relaxation" % "---->"
              % "norm of the static Hamiltonian superoperator: " % L0_norm;
      std::cout
          << boost::format("%=32T %-47s %-10s \n")
              % "correlation function integration time step: " % timestep;
      std::cout << boost::format("%=32T %-47s %-10s \n") % "integration steps: " % nsteps;
      // Get the static Liouvillian propagator
      sp_cx_mat P = propagator(L0, timestep / 4);
      std::cout
          << boost::format("%=20s %=10s %-20s \n") % "Relaxation" % "---->"
              % "integrating spherical component.";
      // Loop over the multi-index
      int k, m, p, q;
      k = m = p = q = -1;
      for (int kmpq = 0; kmpq < 625; kmpq++) {
        std::vector<int> indices = dec2base(kmpq, 5);
        switch (indices.size()) {
          case 4: {
            k = indices[3];
            m = indices[2];
            p = indices[1];
            q = indices[0];
          }
            break;
          case 3: {
            k = 0;
            m = indices[2];
            p = indices[1];
            q = indices[0];
          }
            break;
          case 2: {
            k = 0;
            m = 0;
            p = indices[1];
            q = indices[0];
          }
            break;
          case 1: {
            k = 0;
            m = 0;
            p = 0;
            q = indices[0];
          }
            break;
          default:break;
        }
        // Use the isotropic rotational diffusion model
        double D = 1 / (6 * relax_.tau[0]);
        double liouv_zero = tol_liouv_zero;
        //Compute the term in the rotational expansion sum
        double eps = 2e-16;
        if (norm(Q[k * 5 + m]) > eps && norm(Q[p * 5 + q]) > eps
            && G(k, m, p, q, tau_min, D) > eps) {
          sp_cx_mat A = comp_.mprealloc();
          // Set the initial value for the the operator between the exponentials in the BRW integral
          sp_cx_mat B = Q[p * 5 + q].adjoint();
          std::cout
              << boost::format("%33T k=%2s %5T  m=%2s %5T p=%2s %5T q=%2s\n") % (2 - k)
                  % (2 - m) % (2 - p) % (2 - q);
          // Compute the BRW integral using O(h^7) Boole's quadrature
          for (int i = 0; i < nsteps; i++) {
            A += (timestep / 90) * 7 * G(k, m, p, q, (i + 0.00) * timestep, D)
                * B;

            B = P.adjoint() * B * P;
            cleanup(B, liouv_zero);
            A += (timestep / 90) * 32 * G(k, m, p, q, (i + 0.25) * timestep, D)
                * B;

            B = P.adjoint() * B * P;
            cleanup(B, liouv_zero);
            A += (timestep / 90) * 12 * G(k, m, p, q, (i + 0.50) * timestep, D)
                * B;

            B = P.adjoint() * B * P;
            cleanup(B, liouv_zero);
            A += (timestep / 90) * 32 * G(k, m, p, q, (i + 0.75) * timestep, D)
                * B;

            B = P.adjoint() * B * P;
            cleanup(B, liouv_zero);
            A += (timestep / 90) * 7 * G(k, m, p, q, (i + 1.00) * timestep, D)
                * B;
          }
          sp_cx_mat C = Q[k * 5 + m] * A;
          cleanup(C, liouv_zero);
          mat_R -= C;
        }
      }
      // Decide the fate of dynamic frequency shifts
      if (!relax_.dfs) {
        cx_mat tmp = (mat_R.toDense().real()).cast<cd>();
        mat_R = tmp.sparseView();
        std::cout
            << boost::format("%=20s %=10s %-20s \n") % "Relaxation" % "---->"
                % "dynamic frequency shifts have been ignored.";
      } else
        std::cout
            << boost::format("%=20s %=10s %-20s \n") % "Relaxation" % "---->"
                % "dynamic frequency shifts have been kept.";

      switch (relax_.keep) {
        case diagonal: {
          cx_mat tmp = mat_R.diagonal().asDiagonal();
          mat_R = tmp.sparseView();
          std::cout
              << boost::format("%=20s %=10s %-20s \n") % "Relaxation" % "---->"
                  % "diagonal relaxation superoperator";
        }
          break;
        case kite: {
          /* if (basis_type() != sphten_liouv)
                   fprintf(
                           stderr,
                           "Error:kite option is only available for sphten-liouv formalism!");*/
          mat_R = R2kite(mat_R);
          std::cout
              << boost::format("%=20s %=10s %-20s \n") % "Relaxation" % "---->"
                  % "kite relaxation superoperator";
        }
          break;
        case secular_: {
          /* if (basis_type() != sphten_liouv)
                   fprintf(
                           stderr,
                           "Error:secular option is only available for sphten-liouv formalism!");*/
          mat_R = R2secular(mat_R);
          std::cout
              << boost::format("%=20s %=10s %-20s \n") % "Relaxation" % "---->"
                  % "secular relaxation superoperator";
        }
          break;
        case labframe:
          std::cout
              << boost::format("%=20s %=10s %-20s \n") % "Relaxation" % "---->"
                  % "full relaxation superoperator";
          break;
        default:break;
      }
      if (relax_.equ == levitt) {
        /*        if (sys_.relax_temp() != 0) {
        double beta = kHbar / (sys_.relax_temp() * kbol);
        sys_.set_assumption("labframe");
        sp_cx_mat H;
        sp_cx_mat Q[25];
        H = hamiltonian(Q, sys_, "right");
        R = R * propagator(H, Ci * beta);
        }*/
      }
    }
      break;
    case t1_t2: {
      if (relax_.r1.size() == 0 || relax_.r2.size() == 0) {
        relax_.R = mat_R;
        return;
      }
      //  Preallocate the relaxation superoperator
      cx_mat lm_mat = comp_.basis_state();
      int dim = lm_mat.rows();
      int spins = lm_mat.cols();
      cd val;
      sp_cx_mat R(dim, dim);
      std::vector<double> r1 = relax_.r1;
      std::vector<double> r2 = relax_.r2;
      double current_rlx_rate;

      for (int i = 0; i < dim; i++) {
        current_rlx_rate = 0;
        //  Loop over the constituent spins
        for (int j = 0; j < spins; j++) {
          val = lm_mat(i, j);
          if (val.real() == 0) {
            // L=0
            // Spins in identity states do not contribute
          } else {
            // Longitudinal states contribute R1, transverse states contribute R2
            // M=0
            if (val.imag() == 0)
              current_rlx_rate += r1[j];
            else
              current_rlx_rate += r2[j];
          }
        }
        if (current_rlx_rate != 0)
          R.insert(i, i) = -current_rlx_rate;
      }
      /* std::cout
               << format("%=20s %=10s %-20s \n") % "Relaxation" % "---->"
               % "T1-T2 relaxation superoperator";*/
      mat_R = R;
    }
      break;
    default:break;
  }
  // Levante-Ernst thermalization
  if (relax_.equ == levitt)
    mat_R.col(0) = -mat_R * equilibrium_state();
  relax_.R = mat_R;
}
void interaction::cacl_equilibrium_state() {
  sp_cx_vec rho;
  //switch (basis_type()) {
  //case sphten_liouv: {
  sp_cx_vec unit = comp_.unit_state();
  set_assumption("labframe");
  ham_op ham = hamiltonian(kLeft, 0);
  sp_cx_mat H = ham.isotropic;

  if (relax_.temperature == 0) { // If zero temperature is supplied, use the high temperature approximation
    rho = H * unit;
    rho /= norm(rho);
  } else {
    //double beta = kHbar / (sys.relax_temp() * kbol);
    //rho = expv(-beta * H, unit);
    //cd factor = (rho.cwiseProduct(unit)).norm();  // dot(unit,rho)??
    //rho /= factor;
  }
  //}
  //break;
  //case zeeman_hilb:
  //        break;
  //default:
  //        break;
  //}
  equ_state_ = rho;
}
double interaction::krondelta(int a, int b) const {
  return (a == b) ? 1 : 0;
}
double interaction::G(int k, int m, int p, int q, double tau, double D) const {
  return (1.0 / 5) * krondelta(k, p) * krondelta(m, q) * exp(-6.0 * D * tau);
}
sp_cx_mat interaction::R2kite(const sp_cx_mat &R) const {
  imat M = (comp_.basis_state().imag()).cast<int>();
  M = M.cwiseAbs();
  ivec long_states(M.rows());

  std::vector<int> long_states_idx;
  for (int i = 0; i < M.rows(); i++) {
    long_states[i] = M.row(i).sum();
    if (long_states[i] == 0)
      long_states_idx.push_back(i);
  }

  ivec long_states_idx0(long_states_idx.size());
  for (size_t i = 0; i < long_states_idx.size(); i++)
    long_states_idx0(i) = long_states_idx[i];
  //Index the relaxation superoperator
  ivec rows;
  ivec cols;
  cx_vec vals;
  find(R, rows, cols, vals);
  int num = rows.rows();
  Lia a = ismember(rows, long_states_idx0);
  Lia b = ismember(cols, long_states_idx0);
  ivec ab = a.lta.cwiseProduct(b.lta);
  ivec row_col = rows;
  row_col.setZero();
  for (int i = 0; i < num; i++)
    if (rows(i) == cols(i))
      row_col(i) = 1;

  ivec choice = ab + row_col;
  cx_vec choice_coeff(choice.rows());
  choice_coeff.setZero();
  for (int i = 0; i < choice.rows(); i++)
    if (choice(i))
      choice_coeff(i) = 1;  // reset to 1
  vals = vals.cwiseProduct(choice_coeff);

  sp_cx_mat R_kite(R.rows(), R.rows());
  for (int i = 0; i < num; i++) {
    if (vals(i) != cd(0, 0))
      R_kite.insert((int) rows(i), (int) cols(i)) = vals(i);
  }
  return R_kite;
}
sp_cx_mat interaction::R2secular(const sp_cx_mat &R) const {
  // Compute state carrier frequencies
  mat M = comp_.basis_state().imag();
  mat freqs_mat = comp_.freq_basis();
  freqs_mat = freqs_mat.cwiseProduct(M);

  vec freqs(freqs_mat.rows());
  for (int i = 0; i < freqs_mat.rows(); i++)
    freqs(i) = freqs_mat.row(i).sum();

  std::vector<double> freqs_unique = eigen2stl(freqs);
  sort(freqs_unique.begin(), freqs_unique.end());
  std::vector<double>::iterator it = unique(freqs_unique.begin(),
                                       freqs_unique.end());
  freqs_unique.erase(it, freqs_unique.end());

  //Index the relaxation superoperator
  ivec rows;
  ivec cols;
  cx_vec vals;
  find(R, rows, cols, vals);

  // Set the initial keep mask
  ivec keep_mask(vals.rows());
  keep_mask.setZero();

  for (size_t i = 0; i < freqs_unique.size(); i++) {
    std::vector<int> current_frq_group;
    for (int j = 0; j < freqs.rows(); j++)
      if (freqs(j) == freqs_unique[i])
        current_frq_group.push_back(j);

    ivec current_frqs = stl2eigen(current_frq_group);

    Lia a = ismember(rows, current_frqs);
    Lia b = ismember(cols, current_frqs);
    ivec ab = a.lta.cwiseProduct(b.lta);
    keep_mask += ab;
  }

  cx_vec choice_coeff(keep_mask.rows());
  choice_coeff.setZero();
  for (int i = 0; i < keep_mask.rows(); i++)
    if (keep_mask(i))
      choice_coeff(i) = 1;  // reset to 1
  vals = vals.cwiseProduct(choice_coeff);

  sp_cx_mat R_sec(R.rows(), R.rows());
  for (int i = 0; i < vals.rows(); i++) {
    if (vals(i) != cd(0, 0))
      R_sec.insert((int) rows(i), (int) cols(i)) = vals(i);
  }
  return R_sec;
}
void find(sp_cx_mat m, ivec &rows, ivec &cols, cx_vec &vals) {
  int nnz = m.nonZeros();
  rows.setZero(nnz);
  cols.setZero(nnz);
  vals.setZero(nnz);
  int i = 0;
  for (int k = 0; k < m.outerSize(); ++k)
    for (Eigen::SparseMatrix<cd>::InnerIterator it(m, k); it; ++it) {
      rows(i) = it.row();
      cols(i) = it.col();
      vals(i) = it.value();
      i++;
    }
}
double norm(const sp_cx_mat &m) { // 1-norm
  std::vector<double> cols;
  std::vector<double> sum_cols;
  for (int k = 0; k < m.outerSize(); ++k)
    for (Eigen::SparseMatrix<cd>::InnerIterator it(m, k); it; ++it) {
      int col = it.col();
      double val = sqrt(norm(it.value()));
      typedef std::vector<double>::iterator iter;
      iter i = find(cols.begin(), cols.end(), col);
      if (i != cols.end()) { // already have it
        int pos = -1;
        for (int j = 0; j < (int) cols.size(); j++)
          if (col == cols[j]) {
            pos = j;
            break;
          }
        if (pos != -1)
          sum_cols[pos] += val;
      } else {
        cols.push_back(col);
        sum_cols.push_back(val);
      }
    }
  if (sum_cols.size() == 0)
    return 0;
  else
    return *max_element(sum_cols.begin(), sum_cols.end());
}
ivec stl2eigen(std::vector<int> m) {
  ivec result(m.size());
  for (size_t i = 0; i < m.size(); i++)
    result(i) = m[i];
  return result;
}
std::vector<double> eigen2stl(vec m) {
  std::vector<double> result;
  for (int i = 0; i < m.rows(); i++)
    result.push_back((double) m(i));
  return result;
}
std::vector<int> dec2base(int num, int base) {
  std::vector<int> result;
  while (num >= base) {
    result.push_back(num % base);
    num /= base;
  }
  result.push_back(num);
  return result;
}
void cleanup(sp_cx_mat &m, double nonzero_tol) {
  cd val;
  int i, j;
  for (int k = 0; k < m.outerSize(); ++k)
    for (Eigen::SparseMatrix<cd>::InnerIterator it(m, k); it; ++it) {
      i = it.row();
      j = it.col();
      val = it.value();
      val.real(round(val.real() / nonzero_tol));
      val.imag(round(val.imag() / nonzero_tol));
      m.coeffRef(i, j) = val * nonzero_tol;
    }
}
double maxm(sp_mat m) {
  double max_val = std::numeric_limits<double>::min();
  for (int k = 0; k < m.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(m, k); it; ++it) {

      if (it.value() > max_val)
        max_val = it.value();
    }
  return max_val;
}
}
}
