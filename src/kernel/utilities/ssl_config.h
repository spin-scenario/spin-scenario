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
#ifdef WIN32
#pragma warning (disable: 4625 4526 4710 4711 4820)
#endif

#define _USE_MATH_DEFINES 1
#define NOMINMAX 1
#define SOL_CHECK_ARGUMENTS 1 // enable checking intput arguments of lua functions.
#include <math.h>

#ifdef WIN32
#include <windows.h>
#endif

//#define TENSORFLOW_ENABLED 1
#include "ssl_path.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
using namespace std;
#include <sol.hpp>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/MatrixFunctions>

#include <fftw3.h>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include<yacas/yacas.h>

#include "H5Cpp.h"
using namespace H5;

namespace ssl {
namespace utility {

void ssl_version_output();

// colorful console output.
#ifdef WIN32
extern HANDLE hConsole;
extern WORD wOldColorAttrs;
WORD get_old_color_attrs();

#define SSL_INFO_COLOR1 BACKGROUND_GREEN | BACKGROUND_INTENSITY
#define SSL_WARN_COLOR1 BACKGROUND_RED | BACKGROUND_GREEN | BACKGROUND_INTENSITY
#define SSL_ERR_COLOR1 BACKGROUND_RED | BACKGROUND_INTENSITY
#define SSL_PHASE_COLOR1 BACKGROUND_BLUE | BACKGROUND_INTENSITY

#define SSL_INFO_COLOR2 FOREGROUND_GREEN | FOREGROUND_INTENSITY
#define SSL_WARN_COLOR2 FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_INTENSITY
#define SSL_ERR_COLOR2 FOREGROUND_RED | FOREGROUND_INTENSITY
#define SSL_PHASE_COLOR2 FOREGROUND_BLUE | FOREGROUND_INTENSITY

#endif

void ssl_color_text(const string &option, const string &s, ostream &ostr = cout);

#define SSL_OUTPUT_ENABLE 1

// math definition.
#define _pi M_PI
typedef std::complex<double> cd;

extern double g_inf;
const cd ci = cd(0, 1);

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;
typedef Eigen::ArrayXXd array;

typedef Eigen::Tensor<int, 3> icube;
typedef Eigen::Tensor<double, 3> cube;

typedef Eigen::VectorXi ivec;
typedef Eigen::MatrixXi imat;

typedef Eigen::SparseVector<double> sp_vec;
typedef Eigen::SparseMatrix<double> sp_mat;

typedef Eigen::VectorXcd cx_vec;
typedef Eigen::MatrixXcd cx_mat;

typedef Eigen::SparseVector<cd> sp_cx_vec;
typedef Eigen::SparseMatrix<cd> sp_cx_mat;

typedef Eigen::Matrix3d mat33;
typedef Eigen::Vector2d vec2;
typedef Eigen::Vector2i ivec2;
typedef Eigen::Vector3d vec3;
typedef Eigen::Vector3i ivec3;
typedef Eigen::Vector3cd cx_vec3;
typedef Eigen::Matrix<cd, 5, 1> cx_vec5;
typedef Eigen::Matrix<double, 8, 1> vec8;

extern size_t g_phase_cycle_steps;
extern double g_pw90; // us.
extern double g_max_grad; // mT/m.
extern double g_max_slew_rate; // T/m/s.
extern double g_B0_;
extern int g_openmp_core;
void set_openmp_core(int n);
enum axis {
  _cx = 0,
  _cy = 1,
  _cz = 2
};

double get(size_t i, const vec &v);
double max(const mat &m);

struct seq_param {
  vec2 fov = vec2::Zero(); // unit in mm.
  ivec2 matrix = ivec2::Zero();
  double max_grad = 50; // mT/m.
  double max_slew = 120; // T/m/s.
  double gamma = 42.57; // MHz/T.
  double sw = 50; // kHz.
  double tof = 0;
  string acq_channel = "1H";
  double acq_phase = 0;
  string observer = "";
  void write(ostream &ostr = cout) const {
    string s;
    s = str(boost::format("%70T*\n"));
    ssl_color_text("info", s);
    s = str(boost::format("%s %.f*%.f mm.\n") % "fov set to be" % fov[_cx] % fov[_cy]);
    ssl_color_text("info", s);
    s = str(boost::format("%s %.d*%.d.\n") % "maeigen_readrix size set to be" % matrix[_cx] % matrix[_cy]);
    ssl_color_text("info", s);
    s = str(boost::format("%s %.f mT/m.\n") % "maximum gradient amplitude set to be" % max_grad);
    ssl_color_text("info", s);
    s = str(boost::format("%s %.f T/m/s.\n") % "maximum gradient slew rate set to be" % max_slew);
    ssl_color_text("info", s);
    s = str(boost::format("%70T*\n"));
    ssl_color_text("info", s);
  }
};

struct phantom_space {
  int x0 = -1;
  int x1 = -1;
  int y0 = -1;
  int y1 = -1;
  int z0 = -1;
  int z1 = -1;
  int dx = 1;
  int dy = 1;
  int dz = 1;
};
extern phantom_space g_phantom_space;
void reduce_phantom(const sol::table &t);

extern vec g_expmv_theta;
void load_expmv_theta();

void init_global_lua(sol::state &lua);

extern seq_param *g_seq_param;

// auxiliary functions.
std::map<string, double> phase_map();
extern const std::map<string, double> g_phase_map;

void set_phase_cycle_steps(size_t n);
void set_phase_cycle_steps_api(const sol::table &t);
void set_pw90_api(const sol::table &t);
void set_pw90(double val);
void set_max_grad(double val);
void set_max_slew_rate(double val);

void set_max_grad_api(const sol::table &t);
void set_max_slew_rate_api(const sol::table &t);

void set_grad(double max_amp, double max_slew_rate);
void set_B0(string mag);

void set_B0_api(const sol::table &t);

void apodization(bool app, double decay_rate);

herr_t op_func(hid_t loc_id, const char *name, const H5L_info_t *info, void *operator_data);
void h5write(H5File &file, Group *group, string dataset_name, const string s);
void h5write(H5File &file, Group *group, string dataset_name, const mat &m);
void h5write(H5File &file, Group *group, string dataset_name, const vec &v);
void h5write(H5File &file, Group *group, string dataset_name, const ivec &iv);
void h5write(H5File &file, Group *group, string dataset_name, const icube &cube);
void h5write(H5File &file, Group *group, string dataset_name, const cube & m);

mat h5read_mat(const H5File &file, string dataset_name);
imat h5read_imat(const H5File &file, string dataset_name);
cube h5read_cube(const H5File &file, string dataset_name);
icube h5read_icube(const H5File &file, string dataset_name);

mat h5read_mat(string file, string dataset_name);
// seq time definition.
typedef int timeline;
typedef Eigen::Matrix<timeline, Eigen::Dynamic, 1> tlvec;
const double seq_time_scale = 1; // unit in us.
#define timeline2ms(x/*timeline*/) (double)x*seq_time_scale*1e-3
#define ms2timeline(x/*ms*/) (timeline)((x*1e3+0.5)/seq_time_scale)
//#define Time(x) ((float)((int)((x+0.0005)*1000)))/1000

#define deg2rad(x) x/180*_pi
#define rad2deg(x) x/_pi*180

#define zero_round_off(x) (abs(x) > 1e-8 ? x : 0)

extern int omp_core_num;
extern sol::state *g_lua;
extern CYacas *g_yacas;
extern vector<string> g_h5_string;
#define PATH_SEPARATOR   '/'
#define PATH_SEPARATOR_2 "/"
void declare_path(const char *ptr2);
char *sys_time();

void load_yacas();
// C++ and Lua interface for symbolic-computation.
string yacas_evaluate(const string expr);
double yacas_integral(double from, double to, string func, string var = "x", int precision = 5);
string yacas_integral(string func, string var = "x");
double yacas_func(double pos, string func, string var = "x", int precision = 5);

void yacas_global_vars();

cd amp2xy(cd val);
cd xy2amp(cd val);

cd fid_xy2amp(cd val);
cd fid_amp2xy(cd val);

struct state_par {
  vector<cd> coeff;
  vector<string> expr;
};

state_par state_evaluate(const string expr);

double phase_in_rad(cd val);
double phase_in_degree(cd val);

cx_vec fft(const cx_vec &src, int nfft);

cx_vec fft_1d(const cx_vec &src);
cx_vec fft_1d(const sol::table &t);
cx_mat fft_2d(const cx_mat &src);

mat eigen_read(string file_name);

enum win_shape {
  _rect,
  _hamming,
  _gauss,
  _unknown_window
};
std::map<string, win_shape> win_shape_map();
const std::map<string, win_shape> g_win_shape = win_shape_map();

mat unwrap2d(const mat &wrapped_phi);
vec window_function(win_shape wshape, int length);
win_shape window_interpreter(string win_name);
struct stft_out {
  cx_mat specgram;
  mat specgram_re;
  mat specgram_im;
  mat amp;
  mat ampdB;
  mat phase;
  mat unwrap_phase;
  vec time; // time vector, s
  vec freq; // frequency vector, Hz
  double delta_time;
  double delta_freq;
  stft_out() {};
  stft_out(const cx_mat &a, const vec &b, const vec &c, double d, double e)
      : specgram(a), time(b), freq(c), delta_time(d), delta_freq(e) {};
};
    // http://cn.mathworks.com/matlabcentral/fileexchange/45197-short-time-fourier-transformation--stft--with-matlab-implementation
stft_out stft(const cx_vec &signal, win_shape wshape, int win_length, int hop, int nfft, double fs);

enum cartesian {
  cx = 0,
  cy = 1,
  cz = 2
};

// This is specially designed for Seq Control variable.
enum CtrlID {
  amp = 0,  ///< RF magnitude.
  phi = 1,  ///< RF phase.
  cf = 2,  ///< carrier freq of RF pulse.
  gx = 3,  ///< Grad x axis.
  gy = 4,  ///< Grad y axis.
  gz = 5,  ///< Grad z axis.
  sig = 6,  ///< Acq enabled flag.
  idx = 7  ///< index of fid data.
};
cd traced(const sp_cx_mat &m);
// general for type mat, vec, sp_mat, sp_vec, sp_cx_mat, sp_cx_vec
template<typename T>
T add(const T &a, const T &b) {
  return a + b;
}
template<typename T>
T sub(const T &a, const T &b) {
  return a - b;
}
template<typename T>
T mul(const T &a, const T &b) {
  return a * b;
}
template<typename T>
T mul(const T &a, double b) {
  return a * b;
}
template<typename T>
T mul(double a, const T &b) {
  return a * b;
}
template<typename T>
T div(const T &a, double b) {
  return a * (1 / b);
}
sp_cx_mat mul(const sp_mat &a, cd b);
sp_cx_mat mul(cd a, const sp_mat &b);
sp_cx_mat mul(const sp_cx_mat &a, cd b);
sp_cx_mat mul(cd a, const sp_cx_mat &b);
sp_cx_vec mul(const sp_cx_vec &a, cd b);
sp_cx_vec mul(cd a, const sp_cx_vec &b);
sp_cx_mat div(const sp_mat &a, cd b);
sp_cx_mat div(const sp_cx_mat &a, cd b);
sp_cx_vec div(const sp_cx_vec &a, cd b);

mat table2mat(const sol::table& t, int nrows, int ncols);
vec table2vec(const sol::table &t);
sol::table vec2table(const vec &v);

template<typename T>
void print(sol::variadic_args va, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & /*m*/) {
  cout.precision(3);
  for (auto v : va) {
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> val = v;
    string sep = "\n----------------------------------------\n";
    Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
    cout << val.format(OctaveFmt) << sep;
  }

}
template<typename T>
void print(sol::variadic_args va, const Eigen::Matrix<T, Eigen::Dynamic, 1> & /*v*/) {
  cout.precision(3);
  for (auto v : va) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> val = v;
    string sep = "\n----------------------------------------\n";
       Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
    cout << val.format(OctaveFmt) << sep;
  }
}
template<typename T>
void print(sol::variadic_args va, const Eigen::SparseMatrix<T> & /*m*/) {
  cout.precision(3);
  for (auto v : va) {
    Eigen::SparseMatrix<T> val = v;
    string sep = "\n----------------------------------------\n";
       Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
    cout << val.toDense().format(OctaveFmt) << sep;
  }
}
template<typename T>
void print(sol::variadic_args va, const Eigen::SparseVector<T> & /*v*/) {
  cout.precision(3);
  for (auto v : va) {
    Eigen::SparseVector<T> val = v;
    string sep = "\n----------------------------------------\n";
    Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
    cout << val.toDense().format(OctaveFmt) << sep;
  }
}

void write(string file, sol::variadic_args va, const mat & /*m*/);
void write(string file, sol::variadic_args va, const cx_mat & /*m*/);
void write(string file, sol::variadic_args va, const vec & /*v*/);
void write(string file, sol::variadic_args va, const cx_vec & /*v*/);
void write(string file, sol::variadic_args va, const sp_mat & /*m*/);
void write(string file, sol::variadic_args va, const sp_cx_mat & /*m*/);
void write(string file, sol::variadic_args va, const sp_vec & /*v*/);
void write(string file, sol::variadic_args va, const sp_cx_vec & /*v*/);
template<typename T>
Eigen::SparseVector<T> normalized(const Eigen::SparseVector<T> &v) {
  double norm_val = v.norm();
  return (1 / norm_val) * v;
}

sol::object retrieve_table(string key, const sol::table &t, string supp = "");
string retrieve_table_str(string key, const sol::table &t, string supp = "");
int retrieve_table_int(string key, const sol::table &t, string supp = "");
double retrieve_table_double(string key, const sol::table &t, string supp = "");
size_t retrieve_table_size_t(string key, const sol::table &t, string supp = "");

bool is_retrievable(string key, const sol::table &t);

struct colon_sep {
  double a;
  double b;
  int num;
};

// parse  "a:b:c"
bool parse(string exp, colon_sep &val);

vec stl2vec(const vector<double> &data);

enum WindowFunction {
  kWF_exp_1d, ///<divides the first point by 2 and multiplies the FID by a decaying exponential with the decay rate specified by the user.
  kWF_crisp_1d, ///<divides the first point by 2 and multiplies the FID by a matched cos^8 half-bell.
  kWF_exp_2d
};

/// \brief Performs free induction decay apodization.
void apodization(cx_vec &fid, double decay_rate, const WindowFunction wf = kWF_exp_1d);
}
}

