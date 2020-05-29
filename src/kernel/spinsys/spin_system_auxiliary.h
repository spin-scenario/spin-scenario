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

#pragma once

#include "isotope.h"
#include<kernel/utilities/ssl_config.h>
using namespace ssl::utility;

namespace ssl {
namespace spinsys {

const double kHbar = 1.05457266e-34; // Planck's  constant (Js)
const double kMu = 1.e-7;          // mu_0/(4*pi) in unit J-sec^2/C^2m

struct rf_ham {
  size_t channels;
  sp_cx_mat *Lx;  // channel 1,2,...
  sp_cx_mat *Ly;  // channel 1,2,...
  cx_mat *Lx_dense;  // channel 1,2,...
  cx_mat *Ly_dense;  // channel 1,2,...
  std::vector<std::string> chs;
  void init(std::vector<std::string> ch_list) {
    channels = ch_list.size();
    Lx = new sp_cx_mat[channels];
    Ly = new sp_cx_mat[channels];
    Lx_dense = new cx_mat[channels];
    Ly_dense = new cx_mat[channels];
    chs = ch_list;
  }
  int channel_index(std::string symbol) {
    std::vector<std::string>::iterator pos = find(chs.begin(), chs.end(), symbol);
    if (pos != chs.end())
      return distance(chs.begin(), pos);
    std::string s = "no ** " + symbol + " ** channel for 'rf_ham'";
    throw std::runtime_error(s.c_str());
    return -1;
  }
};
cx_vec acquire(const sp_cx_vec &rho0, const sp_cx_vec &coil, const sp_cx_mat &L, int points, double sw);
sp_cx_vec step(const sp_cx_vec &rho0, const sp_cx_mat &L, double dt);
cx_vec step(const cx_vec &rho0, const cx_mat &L, double dt);

cx_vec expmv(const cx_vec &b0, const cx_mat &A0, double t, const mat &M0, bool shift, bool bal);
cx_mat expmv_tspan(const cx_vec &b0,
                   const cx_mat &A0,
                   double t0,
                   double tmax,
                   int q,
                   const mat &M0,
                   bool shift,
                   bool bal);

mat select_taylor_degree(const cx_mat &A0, bool shift, bool bal, bool force_estm);
void degree_selector(double t, const mat &M, int &m, int &s);

double norm1(const cx_mat &m);
double norm_inf(const cx_mat &m);

double norm1(const sp_cx_mat &m);
double norm_inf(const sp_cx_vec &m);

sp_cx_vec expmv(const sp_cx_vec &b0, const sp_cx_mat &A0, double t, const mat &M0, bool shift, bool bal);
sp_cx_mat expmv_tspan(const sp_cx_vec &b0,
                      const sp_cx_mat &A0,
                      double t0,
                      double tmax,
                      int q,
                      const mat &M0,
                      bool shift,
                      bool bal);
mat select_taylor_degree(const sp_cx_mat &A0, bool shift, bool bal, bool force_estm);

//sp_cx_vec step0(const sp_cx_vec &rho0, const sp_cx_mat &L, double dt);

cd projection(const sp_cx_vec &rho, const sp_cx_vec &coil);
cd projection(const cx_vec &rho, const cx_vec &coil);
double transfer_fidelity(const sp_cx_vec &rho, const sp_cx_vec &targ);
double transfer_fidelity(const sp_cx_mat &rho, const sp_cx_mat &targ);
double max_abs(const sp_cx_vec &m);
sp_cx_mat commutator(const sp_cx_mat &op1, const sp_cx_mat &op2);
sp_cx_mat propagator(const sp_cx_mat &L, double dt);
sp_cx_vec norm_state(const sp_cx_vec &rho);

/////////////////////TODO/////////////////////////////////////////
void levante_ernst_correction(sp_cx_vec &rho);
void levante_ernst_correction(cx_vec &rho);

// spherical tensor.
struct Tlm {
  int L;
  int M;
  Tlm() {}
  Tlm(int l, int m) : L(l), M(m) {
  }
  // converts linear indexing state specification to L, M indexing. In the linear indexing convention,
  // the states are listed in the order of increasing L rank, and, within ranks, in the order of decreasing M projection.
  // I = 0 -> (L = 0, M =  0) => Ie
  // I = 1 -> (L = 1, M =  1) => I+
  // I = 2 -> (L = 1, M =  0) => Iz
  // I = 3 -> (L = 1, M = -1) => I-
  Tlm(int I) {
    L = int(sqrt(I));
    M = L * L + L - I;
  }
  // converts L, M spin state specification to linear indexing specification.
  int linear_indexing() {
    return L * L + L - M;
  }
};
enum op_label {
  kIe = 0,
  kIp = 1,
  kIz = 2,
  kIm = 3
};
// coeff for product superoperator
struct op_coeff {
  size_t id;
  double coeff;
  int linear_index;  // linear indexing specification of spin state.
};
enum op_side {
  kComm, // produces commutation superoperator.
  kLeftComm,
  kRightComm,
  kLeft, // produces left side product superoperator.
  kRight, // produces right side product superoperator.
  kAcomm // produces anticommutation superoperator.
};
// basic pauli matirx construction
struct pauli_matrices {
  sp_mat p;
  sp_mat m;
  sp_mat z;
  pauli_matrices(size_t mult) {
    double qn = double(mult - 1) / 2.0;  // quantum number
    vec prjs = vec::LinSpaced(mult, (double) (mult - 1.0) - qn, 0 - qn);
    p = sp_mat(mult, mult);  //  I+
    m = sp_mat(mult, mult);  //  I-
    z = sp_mat(mult, mult);  //  Iz
    for (int i = 0; i < (int) mult; i++) {
      // I+
      if (i + 1 < (int) mult) {
        double val = qn * (qn + 1) - prjs(i + 1) * (prjs(i + 1) + 1);
        p.insert(i, i + 1) = sqrt(val);
      }
      // I-
      if (i - 1 >= 0) {
        double val = qn * (qn + 1) - prjs(i - 1) * (prjs(i - 1) - 1);
        m.insert(i, i - 1) = sqrt(val);
      }
      // Iz
      z.insert(i, i) = prjs(i);
    }
  }
};

// Matlab function '[Lia,Locb] = ismember(A,B,'rows') '
struct Lia {
  ivec lta;  // contains 1 where the data in A is found in B. Elsewhere, it returns 0.
  ivec locb;  // contains the lowest index in B for each row in A that is also a row in B.
};
template<typename T>
Lia ismember(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &A,
             const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &B) {
  Lia result;
  int ca = A.cols();
  int cb = B.cols();
  // A and B must have the same number of columns
  if (ca != cb)
    throw std::runtime_error("ismember function of matrix A and B should have the same columns!");

  int ra = A.rows();
  int rb = B.rows();
  ivec lta(ra);
  ivec locb(ra);
  for (int i = 0; i < ra; i++) {
    int j;
    for (j = 0; j < rb; j++) {
      if (A.row(i) == B.row(j)) {
        lta(i) = 1;
        locb(i) = j + 1;  // Note pos index like 1,2,3 ...
        break;
      } else
        continue;
    }
    // current row of A is not in B
    if (j == rb) {
      lta(i) = 0;
      locb(i) = 0;
    }
  }
  result.lta = lta;
  result.locb = locb;
  return result;
  /*
  *  Debug code for 'ismember'
  imat a(4,3);
  a<<1,2,3,2,3,3,3,4,8,1,2,3;
  imat b(2,3);
  b<<1 ,2 ,3,3, 4 ,8;
  std::cout<<(a.row(0)==b.row(0))<<"\n";
  Lia result=ismember(a,b);
  std::cout<<a<<"\n\n"<<b<<"\n";
  std::cout<<"\n"<<result.lta<<"\n\n"<<result.locb;
  */
}
template<typename T>
Lia ismember(const Eigen::Matrix<T, Eigen::Dynamic, 1> &A,
             const Eigen::Matrix<T, Eigen::Dynamic, 1> &B) {
  Lia result;
  int ca = A.cols();
  int cb = B.cols();
  // A and B must have the same number of columns
  if (ca != cb)
    throw std::runtime_error("ismember function of matrix A and B should have the same columns!");

  int ra = A.rows();
  int rb = B.rows();
  ivec lta(ra);
  ivec locb(ra);
  for (int i = 0; i < ra; i++) {
    int j;
    for (j = 0; j < rb; j++) {
      if (A.row(i) == B.row(j)) {
        lta(i) = 1;
        locb(i) = j + 1;  // Note pos index like 1,2,3 ...
        break;
      } else
        continue;
    }
    // current row of A is not in B
    if (j == rb) {
      lta(i) = 0;
      locb(i) = 0;
    }
  }
  result.lta = lta;
  result.locb = locb;
  return result;
}
// sparse matirx pow function m^k.
template<typename T>
Eigen::SparseMatrix<T> spow(const Eigen::SparseMatrix<T> &m, int k) {
  if (k == 1)
    return m;
  Eigen::SparseMatrix<T> mat = m;
  if (k == 0) {
    mat.setIdentity();
    return mat;
  }
  return mat * spow(mat, k - 1);
}
// temporary 'trace' method for SparseMatrix
template<typename T>
T trace(const Eigen::SparseMatrix<T> &m) {
  T sum = 0;
  for (int k = 0; k < m.outerSize(); ++k)
    sum += m.coeff(k, k);
  return sum;
}

sp_cx_mat Identity(int n);

class composition {
 public:
  composition();
  ~composition();

  void init();
  inline size_t nspins() const {
    return spins_.size();
  }
  std::string isotopes() const {
    std::string s;
    for (auto v : spins_)
      s += v.symbol() + " ";
    return s;
  }
  std::vector<std::string> channels() const {
    std::vector<std::string> s;
    for (auto v : spins_)
      s.push_back(v.symbol());
    //sort(s.begin(), s.end());
    s.erase(unique(s.begin(), s.end()), s.end());
    return s;
  }

  inline const std::vector<isotope> &spins() const {
    const std::vector<isotope> &ref = spins_;
    return ref;
  }
  inline double get_field() const {
    return B0_;
  }
  double get_proton_freq() const; // unit in MHz.
  inline void add_spin(const isotope &iso) {
    spins_.push_back(iso);
  }
  std::vector<size_t> parse_labeled_spins(const std::string symbol) const;
  // basis table matrix of all spins, each column represents one spin
  imat basis_table() const;
  cx_mat basis_state() const;
  mat freq_basis() const;
  double get_freq(size_t id) const;
  // hilbert space size of spins before this spin
  size_t hilbert_space_dim_forward(size_t id) const;
  // hilbert space size of spins after this spin
  size_t hilbert_space_dim_backward(size_t id) const;
  // hilbert space size of all spins
  size_t hilbert_space_dim() const;
  // hilbert space size of the given spin.
  inline size_t hilbert_space_dim(size_t id) const {
    return spins_[id].hs();
  }
  inline size_t liouville_space_dim() const {
    size_t dim = hilbert_space_dim();
    return dim * dim;
  }

  double gamma(size_t i) const {
    return spins_[i].gamma();
  }

  // Valid labels for operators in this type of call are 'Ie', 'Iz', 'I+', 'I-'  and 'Tl,m'
  // op("0 I- 1 I+", "comm", sys)
  // op("I-","13C","comm", sys)
  // If an operator is given as a single std::string and spins are named by passing a single std::string
  // the function returns the sum of the corresponding single-spin superoperators on all spins
  // with that name.
  sp_mat op(const std::string list, op_side type = kComm) const;
  sp_cx_vec state(const std::string list) const;
  sp_cx_mat mprealloc() const;
  sp_cx_vec unit_state() const;
 private:
  void check_valid_spin_id(const size_t id) const;
  void check_valid_op_label(const std::string label) const;
  sp_mat p_superop(const std::vector<op_coeff> &active_spins,
                   const std::vector<int> &flags, const op_side type) const;
  void ist_product_table(size_t mult, mat *product_table_left,
                         mat *product_table_right) const;
  std::vector<sp_mat> irr_sph_ten2(size_t mult, int k) const;
  std::vector<sp_mat> irr_sph_ten2(size_t mult) const;
  imat basis_cols(std::vector<size_t> active_ids) const;
  imat basis_cols(std::vector<size_t> active_ids, std::vector<int> row_select) const;
  sp_mat empty_op() const;
 public:
  double B0_; // static magnet field, unit in Tesla
  vec base_freq_;
 private:
  std::vector<isotope> spins_;
  imat basis_table_;
};

///////////////////////////////////////////////////////////////////////////////
enum strength_option {
  strong = 1,
  full = 2,
  secular = 3,
  weak = 4
};

struct broadband_cs {
  size_t spin_id; // 0-based indexing.
  vec offset; // unit in ppm.
  vec nominal_offset;
};

struct zeeman {
  std::vector<double> scalars;  // in ppm
  std::vector<broadband_cs> bb_scalars;
  std::vector<vec3> eigs;  // in ppm
  std::vector<vec3> eulers;  // in radian

  std::vector<mat33> matrices;
  std::vector<strength_option> strengths;
  zeeman() {}
  zeeman(size_t nspins) {
    scalars = std::vector<double>(nspins, 0);
    eigs = std::vector<vec3>(nspins, vec3::Zero());
    eulers = std::vector<vec3>(nspins, vec3::Zero());
    matrices = std::vector<mat33>(nspins, mat33::Zero());
    strengths = std::vector<strength_option>(nspins, full);
  }
};
//struct CoupMat {
//        size_t i;
//        size_t j;
//        mat33 m;
//        CoupMat(size_t a, size_t b, mat33 mat) {
//                i = a;
//                j = b;
//                m = mat;
//        }
//};
// spherical tensor descriptors
// interaction -> coupling

struct broadband_jcoup {
  size_t spin_id1;
  size_t spin_id2;
  vec offset; // unit in Hz.
  vec nominal_offset;
};

struct coupling {
  mat scalar;  // Isotropic couplings (in Hz).
  std::vector<broadband_jcoup> bb_scalars;
  //imat eigs_euler;  // flag for
  //vec3 * eigs;
  //vec3 *euler;
  std::vector<vec3> eigs;  // in ppm
  std::vector<vec3> eulers;  // in radian

  std::vector<mat33> matrices;
  std::vector<strength_option> strengths;
  //std::vector<CoupMat> cm;
  //strength_option* strengths;
  std::vector<vec3> coordinates;
  coupling() {
    //strengths = NULL;
    //eigs = NULL;
    //euler = NULL;
  }
  coupling(size_t nspins) {
    scalar = mat::Zero(nspins, nspins);
    //eigs_euler = imat::Zero(nspins, nspins);
    eigs = std::vector<vec3>(nspins * nspins, vec3::Zero());
    eulers = std::vector<vec3>(nspins * nspins, vec3::Zero());
    matrices = std::vector<mat33>(nspins * nspins, mat33::Zero());
    strengths = std::vector<strength_option>(nspins * nspins, full);
    //coordinates = std::vector<vec3>(nspins, vec3::Zero());
    //strengths = NULL;
    //eigs = NULL;
    //euler = NULL;
  }
};
const double tol_inter_cutoff = 1e-5;
const double tol_liouv_zero = 1e-7;
// rotation matrix
mat33 euler2dcm(const vec3 &euler);

struct ham_op {
  sp_cx_mat isotropic;  // isotropic part.
  sp_cx_mat *anisotropic;  // twenty-five matrices giving irreducible components of the anisotropic part.
  ham_op() {
    anisotropic = NULL;
  }
};

struct op_item {
  cd coeff;
  std::string op1;
  std::string op2;
  int id1;
  int id2;
  size_t used;  // 0 for non-used 1 for single spin 2 for spin-spin
  op_item() {
    id1 = id2 = 0;
    coeff = 0;
    used = 0;
  }
  op_item(cd a, std::string b, std::string c, int d, int e) {
    coeff = a;
    op1 = b;
    op2 = c;
    id1 = d;
    id2 = e;
    used = 0;
    if (op1 == "" || op2 == "" || id1 == -1 || id2 == -1)
      used = 1;
    else
      used = 2;
  }
};
struct rank {
  double r0;
  cx_vec3 r1;
  cx_vec5 r2;
};
rank mat2sphten(const mat33 &m);
enum relax_theory {
  DAMP,
  t1_t2,
  redfield,
  LINDBLAD,
  NOTTINGHAM
};
enum redfield_option {
  labframe,
  kite,
  secular_,
  diagonal
};
enum equilibrium {
  zero,
  levitt
};
struct relax {
  relax_theory theory;  // a switch controlling the selection of relaxation method, default is NONE.
  redfield_option
      keep;  // a switch controlling the terms to be kept in the relaxation superoperator (only for redfield).
  bool dfs;  // dynamic frequency shifts 1-keep 0-ignore (only for redfield).
  std::vector<double> tau;  // only for redfield
  double temperature;
  equilibrium equ;
  sp_cx_mat R;

  std::vector<double> r1;
  std::vector<double> r2;
  relax() {
    theory = t1_t2;
    equ = levitt;
    temperature = 0;
    keep = labframe;
    dfs = 0;
  }
};
const double tol_rlx_integration = 1e-4;
const double tol_prop_norm = 1e-9;
const double tol_prop_chop = 1e-10;
double norm(const sp_cx_mat &m);
void find(sp_cx_mat m, ivec &rows, ivec &cols, cx_vec &vals);
ivec stl2eigen(std::vector<int> m);
std::vector<double> eigen2stl(vec m);
std::vector<int> dec2base(int num, int base);
void cleanup(sp_cx_mat &m, double nonzero_tol);
double maxm(sp_mat m);

struct cs_par {
  int id;
  double val;
};

struct jcoup_par {
  int id1;
  int id2;
  double val;
};

class interaction {
 public:
  interaction(const composition &comp);
  ~interaction();

  // flexible for Lua binding. E.g. zeeman = {1, "scalar 1.2 ppm"} or zeeman = {1, "scalar", "1.2 ppm"}
  void set_zeeman(std::vector<std::string> list);
  // flexible for Lua binding. E.g. Jcoupling = {1, 2, "scalar 14 Hz"} or Jcoupling = {"1, 2", "scalar", "14 Hz"}
  void set_Jcoupling(std::vector<std::string> list);
  void set_Jcoupling_coords(const sol::table &t);
  void set_relaxation(std::vector<std::string> list);
  void set_zeeman(std::string list);
  void set_Jcoupling(std::string list);
  void set_relaxation(std::string list);
  void init();
  // ONLY used for broadband chemical shifts or j couplings.
  void init_broadband();
  void alloc();
  void set_assumption(std::string assume = "labframe"); // need to be set before calculate the hamiltonian.
  inline std::string get_assumption() const {
    return assume_;
  }
  ham_op hamiltonian(const op_side type, bool build_aniso) const;

  inline sp_cx_mat relaxation() const {
    return relax_.R;
  }
  inline sp_cx_vec equilibrium_state() const {
    return equ_state_;
  }
  std::vector<std::vector<cs_par>> parsing_zeeman_broadband();
  std::vector<std::vector<jcoup_par>> parsing_Jcoupling_broadband();
 private:
  void cacl_relaxation();
  void cacl_equilibrium_state();
  void set_zeeman_scalar_broadband(size_t id, vec3 val, std::string unit = "ppm");
  void set_zeeman_scalar(size_t id, double val, std::string unit = "ppm");
  void set_Jcoupling_scalar(size_t id1, size_t id2, double val,
                            std::string unit = "hz");
  void set_Jcoupling_scalar_broadband(size_t id1, size_t id2, vec3 val, std::string unit = "hz");
  sp_cx_mat aux_hamiltonian(const std::vector<op_item> &descr,
                            const op_side type) const;
  double krondelta(int a, int b) const;
  double G(int k, int m, int p, int q, double tau, double D) const;
  sp_cx_mat R2kite(const sp_cx_mat &R) const;
  sp_cx_mat R2secular(const sp_cx_mat &R) const;
  //sp_cx_mat propagator(sp_cx_mat &L, double dt)const;
  std::vector<std::vector<cs_par>> parsing_zeeman_broadband(const std::vector<std::vector<cs_par>> &old, const broadband_cs &bb);
  std::vector<std::vector<jcoup_par>> parsing_Jcoupling_broadband(const std::vector<std::vector<jcoup_par>> &old,
                                                        const broadband_jcoup &bb);
 public:
  zeeman zeeman_;
  coupling coupling_;
  relax relax_;
  sp_cx_vec equ_state_;
  const composition &comp_;
  std::string assume_;
};

}
}
