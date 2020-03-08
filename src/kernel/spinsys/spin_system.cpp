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

#include "spin_system.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

namespace ssl {
namespace spinsys {
spin_system::spin_system(const string &filename) : inter_(comp_) {

  sol::state lua;
  lua.script_file(filename);
  sol::table t = lua["sys"];
  set_sys(t);
}
spin_system::spin_system(const sol::table &t) : inter_(comp_) {
  set_sys(t);
}
spin_system::~spin_system() {
}
void spin_system::set_sys(const sol::table &t) {
  if (!t.valid() || t.empty())
    throw std::runtime_error("invalid 'spin_system' table parameters (nil or empty).");

  //try {
//  std::string str = retrieve_table("B0", t).as<string>();
//  boost::to_lower(str);
//  vector<string> str_vec;
//  boost::split(str_vec, str, boost::is_any_of(", "), boost::token_compress_on);
//  double val_B0 = boost::lexical_cast<double>(str_vec[0]);
//  if (str_vec[1] == "tesla" || str_vec[1] == "t")
//    set_magnet_field(val_B0);
//  if (str_vec[1] == "mhz")
//    set_proton_freq(val_B0);

  if (g_B0_ < 0) {
    throw std::runtime_error("static magnetic field not set yet!");
  }

  comp_.B0_ = g_B0_;

  set_isotopes(retrieve_table("spin", t).as<string>());
  //setBasis();
  inter_.alloc();

  if (is_retrievable("zeeman", t))
    inter_.set_zeeman(retrieve_table("zeeman", t).as<string>());

  if (is_retrievable("jcoupling", t))
    inter_.set_Jcoupling(retrieve_table("jcoupling", t).as<string>());

  if (is_retrievable("coord", t))
    inter_.set_Jcoupling_coords(retrieve_table("coord", t).as<sol::table>());

  if (is_retrievable("relaxation", t))
    inter_.set_relaxation(retrieve_table("relaxation", t).as<string>());

  inter_.init();
  inter_.set_assumption("nmr");
  /*} catch (const std::runtime_error& e) {
      string s = str(boost::format("%s\n") % string(e.what()));
      ssl_color_text("err", s);
  }*/
}

void spin_system::set_isotopes(const string list) {
  vector<string> symbol_vec;
  boost::split(symbol_vec, list, boost::is_any_of("\t, "), boost::token_compress_on);
  for (size_t i = 0; i < symbol_vec.size(); i++)
    comp_.add_spin(isotope(symbol_vec[i]));
#ifdef SSL_OUTPUT_ENABLE
  string s = str(boost::format("%s %s.\n") % "isotopes set to be" % list);
  ssl_color_text("info", s);
#endif
  comp_.init();
}

void spin_system::set_magnet_field(double tesla) {
  comp_.B0_ = tesla;
#ifdef SSL_OUTPUT_ENABLE
  string s = str(boost::format("%s %.3f Tesla.\n") % "magnet field set to be" % tesla);
  ssl_color_text("info", s);
#endif
}

void spin_system::set_proton_freq(double MHz) {
  isotope proton("1H");
  double tesla = MHz * 2 * _pi / proton.gamma() * 1e6;
#ifdef SSL_OUTPUT_ENABLE
  string s = str(boost::format("%s %.3f MHz.\n") % "proton resonance frequency set to be" % MHz);
  ssl_color_text("info", s);
#endif
  set_magnet_field(tesla);
}
double spin_system::get_proton_freq() const {
  return comp_.get_proton_freq();
}
sp_mat spin_system::op(const string &list) const {
  return op(list, kComm);
}
sp_mat spin_system::op(const sol::table &t) const {
  return op(t, kComm);
}
sp_mat spin_system::op(const string &list, op_side type) const {
  return comp_.op(list, type);
}
sp_mat spin_system::op(const sol::table &t, op_side type) const {
  if (!t.valid() || t.empty())
    throw std::runtime_error("invalid 'op' table parameters (nil or empty).");

  vector<string> s;
  for (size_t i = 0; i < t.size(); i++) {
    sol::object val = t[i + 1];
    string item;
    switch (val.get_type()) {
      case sol::type::number:item = boost::lexical_cast<string>(val.as<int>());
        break;
      case sol::type::string:item = val.as<string>();
        break;
      default:break;
    }
    s.push_back(item);
  }
  string list = boost::algorithm::join(s, " ");
  return comp_.op(list, type);
}

sp_cx_vec spin_system::equilibrium_state() const {
  return inter_.equilibrium_state();
  }

sp_cx_mat spin_system::smart_op(const string expr) const {
  state_par s = state_evaluate(expr);
  int num = s.coeff.size();
  sp_cx_mat result = s.coeff[0] * op(s.expr[0]);
  for (int i = 1; i < num; i++)
    result += s.coeff[i] * op(s.expr[i]);
  return result;
}

sp_cx_vec spin_system::smart_state(const string expr) const {
  state_par s = state_evaluate(expr);
  int num = s.coeff.size();
  sp_cx_vec result = s.coeff[0] * state(s.expr[0]);
  for (int i = 1; i < num; i++)
    result += s.coeff[i] * state(s.expr[i]);
  return result;
}

map<string, sp_cx_vec> spin_system::cartesian_basis_states() const {
  string base_expr = "(1+Ix+Iy+Iz)";
  int nspins = comp_.nspins();
  string expr_in;
  boost::regex reg("I(\\w)");
  vector<string> exprs(nspins);
  for (int i = 0; i < nspins; i++) {
    exprs[i] = base_expr;
    string id = boost::lexical_cast<string>(i + 1);
    string rep = "I" + id + "$1";
    string s = boost::regex_replace(exprs[i], reg, rep);
    exprs[i] = s;
    //cout << rep<<" "<<exprs[i] << "\n";
  }

  expr_in = exprs[0];
  for (int i = 1; i < nspins; i++)
    expr_in += "*" + exprs[i];

  expr_in = "Simplify(" + expr_in + ")";
  //cout << expr_in<<"\n";
  string expr_out = yacas_evaluate(expr_in);
  boost::erase_all(expr_out, ";");
  //cout << expr_out << "\n";
  //boost::replace_all(expr_out, "J", "I");
  //vector<string> expr_vec;
  vector<string> expr;
  boost::algorithm::split(expr, expr_out, boost::is_any_of("+"));
  expr.pop_back(); // !!! remove 1.
  map<string, sp_cx_vec> basis;
  for (size_t i = 0; i < expr.size(); i++) {
    boost::erase_all(expr[i], "*");
    //cout << i + 1 << " " << expr[i] << "\n";
    basis.insert(pair<string, sp_cx_vec>(expr[i], smart_state(expr[i])));
  }
  /*map<string, sp_cx_vec>::iterator iter;
  for (iter = basis.begin(); iter != basis.end(); iter++) {
      string label = iter->first;
      sp_cx_vec val = iter->second;
      cout << label << "\n" << val << "\n";
  }*/
  return basis;
}

sp_cx_vec spin_system::state(const string &list) const {
  return comp_.state(list);
}
sp_cx_vec spin_system::state(const sol::table &t) const {
  if (!t.valid() || t.empty())
    throw std::runtime_error("invalid 'op' table parameters (nil or empty).");
  vector<string> s;
  for (size_t i = 0; i < t.size(); i++) {
    sol::object val = t[i + 1];
    string item;
    switch (val.get_type()) {
      case sol::type::number:item = boost::lexical_cast<string>(val.as<int>());
        break;
      case sol::type::string:item = val.as<string>();
        break;
      default:break;
    }
    s.push_back(item);
  }
  string list = boost::algorithm::join(s, " ");
  return comp_.state(list);
}

ham_op spin_system::hamiltonian(const op_side type, bool build_aniso) const {
  return inter_.hamiltonian(type, build_aniso);
}

sp_cx_mat spin_system::free_hamiltonian() const {
  return hamiltonian(kComm, 0).isotropic;
}

vector<sp_cx_mat> spin_system::free_hamiltonians() {
  vector<sp_cx_mat> L0s;
  vector<vector<cs_par>> result_cs = inter_.parsing_zeeman_broadband();
  vector<vector<jcoup_par>> result_jcoup = inter_.parsing_Jcoupling_broadband();
  if (result_cs.size())
    cout << "chemical shifts num: " << result_cs.size() << "\n";
  if (result_jcoup.size())
    cout << "j-coupling num: " << result_jcoup.size() << "\n";
  // zeeman case.
  vector<double> scalars = inter_.zeeman_.scalars; // backup.
  mat scalar = inter_.coupling_.scalar;
  if (result_cs.size()&&result_jcoup.size()==0)
  for (size_t i = 0; i < result_cs.size(); i++) {
    for (size_t j = 0; j < result_cs[i].size(); j++) {
      const cs_par &cur = result_cs[i][j];
      //cout << (cur.id + 1) << " " << inter_.zeeman_.scalars[cur.id] << " ";
      inter_.zeeman_.scalars[cur.id] += cur.val;
      //cout << (cur.id + 1) << " " << inter_.zeeman_.scalars[cur.id] << " ";
    }
    //cout << "\n";
    inter_.init_broadband();
    L0s.push_back(free_hamiltonian());
    inter_.zeeman_.scalars = scalars; // recover to original pars for the next L0 calculation.
  }
  // J-coupling case.
    if (result_cs.size()==0&&result_jcoup.size())
  for (size_t i = 0; i < result_jcoup.size(); i++) {
    for (size_t j = 0; j < result_jcoup[i].size(); j++) {
      const jcoup_par &cur = result_jcoup[i][j];
      int id1 = cur.id1;
      int id2 = cur.id2;
      inter_.coupling_.scalar(id1, id2) += cur.val;
      if (id1 != id2)
        inter_.coupling_.scalar(id2, id1) += cur.val;
      //cout << (cur.id1 + 1) << " " << (cur.id2 + 1) << " " << inter_.coupling_.scalar(id1, id2) << " ";
    }
    //cout << "\n";
    inter_.init_broadband();
    L0s.push_back(free_hamiltonian());
    inter_.coupling_.scalar = scalar;
  }
	//omp_set_num_threads(omp_core_num);
	if(result_cs.size()&&result_jcoup.size()) {
    //#pragma omp parallel for    
    for (size_t i = 0; i < result_cs.size(); i++) {
      if (i % 100 == 0) cout << i << "==>\n";

      for (size_t j = 0; j < result_cs[i].size(); j++) {
        const cs_par &cur = result_cs[i][j];
        inter_.zeeman_.scalars[cur.id] += cur.val;
      }

      for (size_t ii = 0; ii < result_jcoup.size(); ii++) {

        for (size_t jj = 0; jj < result_jcoup[ii].size(); jj++) {
          const jcoup_par &cur = result_jcoup[ii][jj];
          int id1 = cur.id1;
          int id2 = cur.id2;
          inter_.coupling_.scalar(id1, id2) += cur.val;
          if (id1 != id2) inter_.coupling_.scalar(id2, id1) += cur.val;
        }

        inter_.init_broadband();
        L0s.push_back(free_hamiltonian());

        inter_.coupling_.scalar = scalar;
      }
      inter_.zeeman_.scalars = scalars;  // recover to original pars for the next L0 calculation.
    }
	}
  return L0s;
}

sol::object spin_system::free_hamiltonians(sol::this_state s) {
  sol::state_view lua(s);
  sol::table t = lua.create_table();
  vector<sp_cx_mat> L0s = free_hamiltonians();
  for (size_t i = 0; i < L0s.size(); i++)
    t.add(L0s[i]);
  return t;
}
sp_cx_mat spin_system::relaxation() const {
  return inter_.relaxation();
}

sp_cx_mat spin_system::total_hamiltonian() const {
  return free_hamiltonian() + ci * relaxation();
}
rf_ham spin_system::rf_hamiltonian() const {
  rf_ham rf_ctrl;
  rf_ctrl.init(channels());
  for (size_t i = 0; i < rf_ctrl.channels; i++) {
    // Get the control operators
    string list_plus = rf_ctrl.chs[i] + " I+";
    string list_minus = rf_ctrl.chs[i] + " I-";
    sp_mat Lp = op(list_plus);
    sp_mat Lm = op(list_minus);
    rf_ctrl.Lx[i] = 0.5 * (Lp + Lm).cast<cd>();  // cx
    rf_ctrl.Ly[i] = -0.5 * ci * (Lp - Lm).cast<cd>();  // cy
  }
  return rf_ctrl;
}
vec spin_system::nominal_broadband() {
  if (inter_.zeeman_.bb_scalars.size())
    return inter_.zeeman_.bb_scalars[0].nominal_offset;
  if (inter_.coupling_.bb_scalars.size())
    return inter_.coupling_.bb_scalars[0].nominal_offset;
  return vec(1);
}

}
}
