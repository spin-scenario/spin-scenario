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

#include "phantom.h"
#include <boost/random.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <gnuplot-iostream.h>
#include <time.h>
#include <kernel/utilities/ssl_plot.h>

using namespace ssl::utility;
boost::mt19937 gen;
namespace ssl {
namespace sample {

phantom::phantom() {
}
phantom::phantom(std::string filename, std::string supp_filename) {
  load(filename, supp_filename);
  init_ensemble();
}
phantom::~phantom() {
}
void phantom::init_ensemble() {
  int x0 = -1;
  int x1 = -1;
  int y0 = -1;
  int y1 = -1;
  int z0 = -1;
  int z1 = -1;

  if (g_phantom_space.x0 < 0)
    x0 = 0;
  else
    x0 = g_phantom_space.x0;

  if (g_phantom_space.x1 < 0)
    x1 = dim_[cx] - 1;
  else
    x1 = g_phantom_space.x1;

  if (g_phantom_space.y0 < 0)
    y0 = 0;
  else
    y0 = g_phantom_space.y0;

  if (g_phantom_space.y1 < 0)
    y1 = dim_[cy] - 1;
  else
    y1 = g_phantom_space.y1;

  if (g_phantom_space.z0 < 0)
    z0 = 0;
  else
    z0 = g_phantom_space.z0;

  if (g_phantom_space.z1 < 0)
    z1 = dim_[cz] - 1;
  else
    z1 = g_phantom_space.z1;

  omp_set_num_threads(omp_core_num);
  std::vector<isochromat> *omp_isochromats = new std::vector<isochromat>[omp_core_num];
#pragma omp parallel for
      for (int nz = z0; nz <= z1; nz += g_phantom_space.dz) {
        int id = omp_get_thread_num();
        for (int nx = x0; nx <= x1; nx += g_phantom_space.dx)
          for (int ny = y0; ny <= y1; ny += g_phantom_space.dy) {
            {
              int idx = tissue_dist_(ny, nx, nz);
              // do not deal with background and those tissues not assigned.
              if (idx != 0 && tissue_index_(idx) == 1) {
                isochromat iso;
                iso.data[dB0] = dB0_(ny, nx, nz);
                iso.data[r1] = 1.0 / T1_(ny, nx, nz);
                iso.data[r2] = 1.0 / T2_(ny, nx, nz);
                //iso.data[r2s] = 1.0 / tissue_t1t2_(idx, 1);
                iso.data[cx] = (nx - dim_[cx] / 2 + 0.5) * res_[cx] + offset_[cx];
                iso.data[cy] = (ny - dim_[cy] / 2 + 0.5) * res_[cy] + offset_[cy];
                iso.data[cz] = (nz - dim_[cz] / 2 + 0.5) * res_[cz] + offset_[cz];
                omp_isochromats[id].push_back(iso);
              }

            }
          }
      }
  isochromats_.clear();
  for (int id = 0; id < omp_core_num; id++) {
    std::copy(omp_isochromats[id].begin(), omp_isochromats[id].end(),
              std::back_inserter(isochromats_));
    omp_isochromats[id].clear();
  }
  #ifdef SSL_OUTPUT_ENABLE
  std::string model;
  std::string s = str(boost::format("%s %s.\n") % "total isochromats:" % isochromats_.size());
  ssl_color_text("info", s);
#endif

}
void phantom::view(const sol::table &t) const {
  if (!t.valid() || t.empty())
    throw std::runtime_error("invalid phantom 'view' table parameters (nil or empty).");
  std::string axis;
  int slice = -1;
  std::string prop ="T2";
   if (is_retrievable("prop", t)) 
     prop = retrieve_table_str("prop", t, "phantom view");
 

  for (auto &kv : t) {
    axis = kv.first.as<std::string>();
    if (axis == "prop") continue;
    slice = kv.second.as<int>();
    view(axis, slice, prop);
  }
}
void phantom::view(const std::string &axis, int slice, std::string prop) const {
  int dim = -1; //z
  int max_dim = 0;
  std::string X, Y;
  if (axis == "x" || axis == "X") {
    dim = 1;
    max_dim = dim_[_cx];
    X = "Z";
    Y = "Y";
  }
  if (axis == "y" || axis == "Y") {
    dim = 0;
    max_dim = dim_[_cy];
    X = "Z";
    Y = "X";
  }
  if (axis == "z" || axis == "Z") {
    dim = 2;
    max_dim = dim_[_cz];
    X = "X";
    Y = "Y";
  }
  if (dim == -1) {
    std::string s = "unkown axis: " + axis + "!\n";
    throw std::runtime_error(s.c_str());
  }
  if (slice > max_dim) {
    std::string s = "slice num overflow: " + std::to_string(slice) + " max slice: " + std::to_string(max_dim)+"\n";
    ssl_color_text("err", s);
    return;
  }

  Eigen::Tensor<double, 2> sub;

  if(prop=="T2")
  sub = T2_.chip(slice - 1, dim);

  if(prop=="T1")
  sub = T1_.chip(slice - 1, dim);

  Eigen::Map<mat> m(sub.data(), sub.dimension(0), sub.dimension(1));
  utility::map transfer_map(m.matrix().cast<double>()*1e3); // ms.
  (*g_lua)["_map"] = transfer_map;
  std::string AXIS = axis;
  boost::to_upper(AXIS);
  std::string title = prop + " map of "+ X + Y + " plane view (" + AXIS + "=" + std::to_string(slice) + "/" + std::to_string(max_dim) + ")";
  g_lua->script("plot('title<" + title + "> xlabel<" + X + "> ylabel<" + Y + "> color<Spectral> gnuplot<set size ratio -1>', _map)");
}
void phantom::load(std::string filename, std::string supp_filename) {
  //model_ = unidentified_phantom;

  H5File file;
  file.openFile(filename.c_str(), H5F_ACC_RDONLY);
  tissue_dist_ = h5read_icube(file, "/phantom/tissue_dist");

  dim_[_cy] = tissue_dist_.dimension(0);
  dim_[_cx] = tissue_dist_.dimension(1);
  dim_[_cz] = tissue_dist_.dimension(2);

  mat res = h5read_mat(file, "/phantom/resolution");
  res_.setZero();
  res_[_cx] = res(0);
  res_[_cy] = res(1);
  if(res.size()==3) // for 3d case.
	  res_[_cz] = res(2);
	else
	  res_[_cz] = 1; // 2d case.

  offset_.setZero();

  dB0_ = cube(dim_[_cy], dim_[_cx], dim_[_cz]); 
  dB0_.setZero();

  pd_ = cube(dim_[_cy], dim_[_cx], dim_[_cz]);
  pd_.setConstant(1);

  T1_ = cube(dim_[_cy], dim_[_cx], dim_[_cz]); // rows->y cols->x 
  T1_.setZero();

  T2_ = cube(dim_[_cy], dim_[_cx], dim_[_cz]);
  T2_.setZero();

  T2s_ = cube(dim_[_cy], dim_[_cx], dim_[_cz]);
  T2s_.setZero();


      // load T1/T2 parameters.
      std::string path = utility::g_install_dir + "/share/spin-scenario/config/mida_1.5t_relaxation.dat";
      if (!supp_filename.empty()) path = supp_filename;
      // col 1: index No.
	  // col 2: T1 (ms)
      // col 3: T2 (ms)
	  mat par = eigen_read(path); 
	  if (par.cols()!=3)
		  throw std::runtime_error("Unrecognized phantom T1/T2 parameter file!");

	  int max_idx = (int)par.col(0).maxCoeff()+ 1;
      tissue_index_ = ivec::Zero(max_idx); // indexing based 1.
      tissue_t1t2_ = mat::Zero(max_idx, 2); // indexing based 1.

      for (int i = 0; i < par.rows(); i++) {
        int idx = par(i, 0);
        tissue_index_(idx) = 1;  // effective tissue.
        tissue_t1t2_(idx, 0) = par(i, 1);
        tissue_t1t2_(idx, 1) = par(i, 2);
      }
      tissue_t1t2_ *= 1e-3;  // ms into unit s.	   
      omp_set_num_threads(omp_core_num);
      #pragma omp parallel for
      for (int nz = 0; nz < dim_[_cz]; nz++) {
        for (int nx = 0; nx < dim_[_cx]; nx++)
          for (int ny = 0; ny < dim_[_cy]; ny++) {
            {
              int idx = tissue_dist_(ny, nx, nz);
              // do not deal with background and those tissues not assigned.
              if (idx != 0 && tissue_index_(idx) == 1) {
                T1_(ny, nx, nz) = tissue_t1t2_(idx, 0);
                T2_(ny, nx, nz) = tissue_t1t2_(idx, 1);
              }

            }
          }
      }

#ifdef SSL_OUTPUT_ENABLE
      std::string name;
      //name = h5read_string(file, "/phantom/name");
      std::string s = name + " phantom loaded, cube dimension x/y/z(" +
                      std::to_string(dim_[cx]) + " " +std::to_string(dim_[cy]) + " " 
                      +std::to_string(dim_[cz])+ ")\n";
      ssl_color_text("info", s);
#endif
    }
}
}
