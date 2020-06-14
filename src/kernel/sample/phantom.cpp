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
phantom::phantom(const char *filename) {
  load(filename);
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
  switch (model_) {
    case usr_phantom:
    case mni_brain: {
#pragma omp parallel for
      for (int nz = z0; nz <= z1; nz += g_phantom_space.dz) {
        int id = omp_get_thread_num();
        Eigen::Tensor<double, 2> sub_pd = pd_.chip(nz, 2);
        Eigen::Map<mat> m_pd(sub_pd.data(), sub_pd.dimension(0), sub_pd.dimension(1));

        Eigen::Tensor<double, 2> sub_T1 = T1_.chip(nz, 2);
        Eigen::Map<mat> m_T1(sub_T1.data(), sub_T1.dimension(0), sub_T1.dimension(1));

        Eigen::Tensor<double, 2> sub_T2 = T2_.chip(nz, 2);
        Eigen::Map<mat> m_T2(sub_T2.data(), sub_T2.dimension(0), sub_T2.dimension(1));
        //row = "y";
        //col = "x";

        for (int nx = x0; nx <= x1; nx += g_phantom_space.dx)
          for (int ny = y0; ny <= y1; ny += g_phantom_space.dy) {
            {
              if (m_pd(ny, nx)) {  // only if the spin density is non-zero.
                isochromat iso;
                iso.data[dB0] = dB0_(nx, ny, nz);
                iso.data[r1] = 1.0 / m_T1(ny, nx);
                iso.data[r2] = 1.0 / m_T2(ny, nx);
                //iso.data[r2s] = 1.0 / T2s_(nx, ny, nz);
                iso.data[cx] = (nx - dim_[cx] / 2 + 0.5) * res_[cx] + offset_[cx];
                iso.data[cy] = (ny - dim_[cy] / 2 + 0.5) * res_[cy] + offset_[cy];
                iso.data[cz] = (nz - dim_[cz] / 2 + 0.5) * res_[cz] + offset_[cz];
                iso.data[pd] = m_pd(ny, nx);
                omp_isochromats[id].push_back(iso);
                //}

              }
            }

          }
      }
    }
      break;
    case mida_brain: {
      //clock_t  clockBegin, clockEnd;
      //clockBegin = clock();
#pragma omp parallel for
      for (int nz = z0; nz <= z1; nz += g_phantom_space.dz) {
        int id = omp_get_thread_num();
        for (int nx = x0; nx <= x1; nx += g_phantom_space.dx)
          for (int ny = y0; ny <= y1; ny += g_phantom_space.dy) {
            {
              int idx = tissue_dist_(nx, ny, nz);
              // do not deal with background and those tissues not assigned.
              if (idx != 50 && tissue_index_(idx) == 1) {
                isochromat iso;
                iso.data[dB0] = dB0_(nx, ny, nz);
                iso.data[r1] = 1.0 / T1_(nx, ny, nz);
                iso.data[r2] = 1.0 / T2_(nx, ny, nz);
                //iso.data[r2s] = 1.0 / tissue_t1t2_(idx, 1);
                iso.data[cx] = (nx - dim_[cx] / 2 + 0.5) * res_[cx] + offset_[cx];
                iso.data[cy] = (ny - dim_[cy] / 2 + 0.5) * res_[cy] + offset_[cy];
                iso.data[cz] = (nz - dim_[cz] / 2 + 0.5) * res_[cz] + offset_[cz];
                omp_isochromats[id].push_back(iso);
              }

            }
          }
      }
      //clockEnd = clock();
      //printf("%d\n", clockEnd - clockBegin);
    }
      break;
    default:break;
  }

  isochromats_.clear();
  for (int id = 0; id < omp_core_num; id++) {
    std::copy(omp_isochromats[id].begin(), omp_isochromats[id].end(),
              std::back_inserter(isochromats_));
    omp_isochromats[id].clear();
  }
  #ifdef SSL_OUTPUT_ENABLE
  std::string model;

  if (model_ == mida_brain) model = "MIDA";
  if (model_ == mni_brain) model = "MNI";
  if (model_ == usr_phantom) model = "USR PHANTOM";
  std::string s = str(boost::format("%s %s (%s).\n") % "total isochromats:" % isochromats_.size() % model);
  ssl_color_text("info", s);
#endif

}
void phantom::view(const sol::table &t) const {
  if (!t.valid() || t.empty())
    throw std::runtime_error("invalid phantom 'view' table parameters (nil or empty).");
  std::string axis;
  int slice = -1;
  for (auto &kv : t) {
    axis = kv.first.as<std::string>();
    slice = kv.second.as<int>();
    view(axis, slice);
  }
}
void phantom::view(const std::string &axis, int slice) const {
  int dim = -1; //z
  int max_dim = 0;
  std::string row, col;
  if (axis == "x" || axis == "X") {
    dim = 0;
    max_dim = dim_[_cx];
    row = "z";
    col = "y";
  }
  if (axis == "y" || axis == "Y") {
    dim = 1;
    max_dim = dim_[_cy];
    row = "z";
    col = "x";
  }
  if (axis == "z" || axis == "Z") {
    dim = 2;
    max_dim = dim_[_cz];
    row = "y";
    col = "x";
  }
  if (dim == -1) {
    std::string s = "unkown axis: " + axis + "!\n";
    throw std::runtime_error(s.c_str());
  }

  Eigen::Tensor<double, 2> sub = T2_.chip(slice - 1, dim);
  Eigen::Map<mat> m(sub.data(), sub.dimension(0), sub.dimension(1));
  utility::map transfer_map(m.matrix().cast<double>());
  (*g_lua)["_map"] = transfer_map;
  std::string title = row + col + " plane view (" + axis + "=" + std::to_string(slice) + "/" + std::to_string(max_dim) + ")";
  g_lua->script("plot('title<" + title + "> xlabel<" + row + "> ylabel<" + col + ">', _map)");
}
void phantom::load(const char *filename) {
  model_ = unidentified_phantom;
  H5File file;
  file.openFile(filename, H5F_ACC_RDWR);

  imat dim = h5read_imat(file, "/phantom/dimension");
  dim_[_cx] = dim(0);
  dim_[_cy] = dim(1);
  dim_[_cz] = dim(2);

  mat res = h5read_mat(file, "/phantom/resolution");
  res_[_cx] = res(0);
  res_[_cy] = res(1);
  res_[_cz] = res(2);

  imat id_vec = h5read_imat(file, "/phantom/model_id");
  int id = id_vec(0);

  offset_.setZero();

  dB0_ = cube(dim_[_cx], dim_[_cy], dim_[_cz]);
  dB0_.setZero();

  pd_ = cube(dim_[_cx], dim_[_cy], dim_[_cz]);
  pd_.setConstant(1);

  T1_ = cube(dim_[_cx], dim_[_cy], dim_[_cz]);
  T1_.setZero();

  T2_ = cube(dim_[_cx], dim_[_cy], dim_[_cz]);
  T2_.setZero();

  T2s_ = cube(dim_[_cx], dim_[_cy], dim_[_cz]);
  T2s_.setZero();

  switch (id) {
    case 1: {
      model_ = mida_brain;
      tissue_dist_ = h5read_icube(file, "/phantom/tissue_dist");

      // load T1/T2 parameters.
      std::string path = utility::g_install_dir + "/share/spin-scenario/config/mida_1.5t_relaxation.dat";
      mat par = mat::Zero(59, 3);
      std::ifstream fin(path.c_str());

      if (fin.is_open()) {
        for (size_t row = 0; row < 59; row++)
          for (size_t col = 0; col < 3; col++) {
            double item = 0.0;
            fin >> item;
            par(row, col) = item;
          }
        fin.close();
      }

      //mat par = eigen_read(path);
      tissue_index_ = ivec::Zero(117); // indexing based 1.
      tissue_t1t2_ = mat::Zero(117, 2); // indexing based 1.

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
              int idx = tissue_dist_(nx, ny, nz);
              // do not deal with background and those tissues not assigned.
              if (idx != 50 && tissue_index_(idx) == 1) {
                T1_(nx, ny, nz) = tissue_t1t2_(idx, 0);
                T2_(nx, ny, nz) = tissue_t1t2_(idx, 1);
              }

            }
          }
      }

#ifdef SSL_OUTPUT_ENABLE
      std::string s = str(boost::format("%s%s %s %s).\n") % "MIDA brain loaded, cube dimension x/y/z(" % dim_[cx] % dim_[cy]
                         % dim_[cz]);
      ssl_color_text("info", s);
#endif
    }
      break;
    case 2: {
      model_ = mni_brain;

      T1_ = h5read_cube(file, "/phantom/T1");
      T2_ = h5read_cube(file, "/phantom/T2");
      T2s_ = h5read_cube(file, "/phantom/T2*");
      pd_ = h5read_cube(file, "/phantom/spin density");
#ifdef SSL_OUTPUT_ENABLE
      std::string s = str(boost::format("%s%s %s %s).\n") % "MNI brain loaded, cube dimension x/y/z(" % dim_[cx] % dim_[cy]
                         % dim_[cz]);
      ssl_color_text("info", s);
#endif
    }
      break;
    default: {
      model_ = usr_phantom;
      T1_ = h5read_cube(file, "/phantom/T1");
      T2_ = h5read_cube(file, "/phantom/T2");
#ifdef SSL_OUTPUT_ENABLE
      std::string s = str(boost::format("%s%s %s %s).\n") %
                     "usr phantom loaded, cube dimension x/y/z(" % dim_[cx] %
                     dim_[cy] % dim_[cz]);
      ssl_color_text("info", s);
#endif
    }
		break;
  }

  if (model_ == unidentified_phantom)
    throw std::runtime_error("unidentified phantom!");

  //	boost::mt19937 rng;
  //	boost::uniform_01<boost::mt19937> zeroone(rng);
  //	cube::extent_gen extents;
  //	dB0_.resize(extents[dim_[cx]][dim_[cy]][dim_[cz]]);
  //	//std::fill_n(dB0_.data(), dB0_.num_elements(), 0);
  //	for (auto a : dB0_) for (auto b : a) for (auto&c : b) c = zeroone();

  //	/*std::cout << dB0_[3][5][8] << "\n";
  //	std::cout << dB0_[3][2][6] << "\n";
  //	std::cout << dB0_[3][78][8] << "\n";*/
}
}
}
