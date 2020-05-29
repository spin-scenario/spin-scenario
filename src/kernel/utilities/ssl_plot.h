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

#include "ssl_config.h"
#include <gnuplot-iostream.h>

//#define GP_SCRIPT_OUTPUT 1

namespace ssl {
namespace utility {
extern std::string g_output_terminal;
extern std::string g_output_font;

void set_output_terminal(const sol::table &t);
std::string terminal_cmd(std::string key);

vec linspace(double start, double stop, int num);

sol::object linspace_lua(double start, double stop, int num);

mat random(int rows, int cols);

mat vec_table(const sol::table &t);

struct line {
  vec x;
  vec y;
  std::string file = "";
  std::string line_spec = "";
  bool is_x_only = false;
  bool is_file = false;

  line(std::string a) {
    file = a;
    is_file = true;
  }

  line(std::string a, std::string b) {
    file = a;
    line_spec = b;
    is_file = true;
  }

  line(const vec &a) {
    x = a;
    is_x_only = true;
  }

  line(const vec &a, std::string b) {
    x = a;
    line_spec = b;
    is_x_only = true;
  }

  line(const vec &a, const vec &b) {
    x = a;
    y = b;
  }

  line(const vec &a, const vec &b, std::string c) {
    x = a;
    y = b;
    line_spec = c;
  }

  line(const sol::table &t) {
    init_line(t);
    is_x_only = true;
  }

  line(const sol::table &t, std::string b) {
    init_line(t);
    line_spec = b;
    is_x_only = true;
  }

  line(const sol::table &t0, const sol::table &t1) {
    init_line(t0, t1);
  }

  line(const sol::table &t0, const sol::table &t1, std::string c) {
    init_line(t0, t1);
    line_spec = c;
  }

  void init_line(const sol::table &t0) {
    int np0 = t0.size();
    x = vec(np0);
    for (int i = 0; i < np0; i++) {
      sol::object val = t0[i + 1];
      x[i] = val.as<double>();
    }
  }

  void init_line(const sol::table &t0, const sol::table &t1) {
    int np0 = t0.size();
    int np1 = t1.size();
    x = vec(np0);
    for (int i = 0; i < np0; i++) {
      sol::object val = t0[i + 1];
      x[i] = val.as<double>();
    }
    y = vec(np1);
    for (int i = 0; i < np1; i++) {
      sol::object val = t1[i + 1];
      y[i] = val.as<double>();
    }
  }
};

struct line_series {
  vec x;
  std::vector<vec> y;
  std::vector<std::string> files;
  bool is_y_only = false;
  bool is_file = false;

  line_series(std::string file/*file_list*/) {
    //boost::cmatch what;
    //boost::regex reg("([^\"]*)([^\\d]*)(\\d+):(\\d+)([^\\d]*)([^\"]*)");
    //std::string s1, s2, s3, s4;
    //int n0 = 0, n1 = 0;
    //if (boost::regex_search(file_list.c_str(), what, reg)) {
    //  //std::cout << what[0] << "\n" << what[1] << "\n" << what[2] << "\n" << what[3] << "\n" << what[4] << "\n" << what[5] << "\n" << what[6] << "\n";
    //  s1 = what[1];
    //  s2 = what[2];
    //  n0 = boost::lexical_cast<int>(what[3]);
    //  n1 = boost::lexical_cast<int>(what[4]);
    //  s3 = what[5];
    //  s4 = what[6];
    //}

    //for (int i = 0; i < (n1 - n0 + 1); i++)
    //  files.push_back(s1 + s2 + boost::lexical_cast<std::string>(n0 + i) + s3 + s4);

    //is_file = true;
	  mat data = eigen_read(file);
    for (size_t i = 0; i < data.cols(); i++) 
		y.push_back(data.col(i));
    
	is_y_only = true;
  }

  line_series(const sol::table &t) {
    int np = t.size();
    for (int i = 0; i < np; i++) {
      sol::object val = t[i + 1];
      y.push_back(val.as<vec>());
    }
    is_y_only = true;
  }

  line_series(const vec &a, const sol::table &t) {
    x = a;
    int np = t.size();
    for (int i = 0; i < np; i++) {
      sol::object val = t[i + 1];
      y.push_back(val.as<vec>());
    }
  }
};

struct map {
  mat m;
  std::string file = "";
  std::string map_spec = "";
  bool is_file = false;
  vec2 xrange;
  vec2 yrange;

  map(std::string a) {
    file = a;
    m = eigen_read(a);
    init_range();
    is_file = true;
  }

  map(const mat &a) {
    m = a;
    init_range();
  }
  map(const vec &x, const vec &y, const mat &a) {

  }

  map(const mat &a, std::string b) {
    m = a;
    init_range(b);
  }

  map(std::string a, std::string b) {
    file = a;
    m = eigen_read(a);
    init_range(b);
    map_spec = b;
    is_file = true;
  }

  void init_range(std::string par = "") {
    boost::cmatch what;
    boost::regex reg_xrange("xrange<(\\-?\\d+\\.?\\d*):(\\-?\\d+\\.?\\d*)>{1}");
    if (boost::regex_search(par.c_str(), what, reg_xrange)) {
      xrange[0] = boost::lexical_cast<double>(what[1]);
      xrange[1] = boost::lexical_cast<double>(what[2]);
    } else {
      xrange[0] = 1;
      xrange[1] = m.cols();
    }

    boost::regex reg_yrange("yrange<(\\-?\\d+\\.?\\d*):(\\-?\\d+\\.?\\d*)>{1}");
    if (boost::regex_search(par.c_str(), what, reg_yrange)) {
      yrange[0] = boost::lexical_cast<double>(what[1]);
      yrange[1] = boost::lexical_cast<double>(what[2]);
    } else {
      yrange[0] = 1;
      yrange[1] = m.rows();
    }

    map_spec = par;
  }
};

// plot one or mutil-line in the same figure.
void plot(sol::variadic_args va, const line &);

void plot(std::string fig_info, sol::variadic_args va, const line &);

void plot(const line_series &);

void plot(std::string fig_info, const line_series &);

void plot(std::string fig_info, sol::variadic_args va, const line_series &);

// plot one or mutil-map figures.
void plot(sol::variadic_args va, const map &);

void plot(std::string fig_info, sol::variadic_args va, const map &);

// auxiliary functions.
std::map<std::string, std::string> line_color();

std::map<std::string, std::string> line_style();

std::map<std::string, std::string> line_marker();

std::map<std::string, std::string> line_dash();

std::map<std::string, std::string> map_color();

void write_line(const vec &x, std::string file);

void write_line(const vec &x, const vec &y, std::string file);

void write_map(const mat &x, std::string file);

struct line_spec {
  std::string color; // default line_color.
  std::string style = "l"; // default line_style.
  std::string marker; // default line_marker.
  std::string dash; // default line_marker.
};

struct map_spec {
  std::string color = "Spectral"; // default map_color.
  std::string style = "image";
};

struct fig_spec {
  std::string title;
  std::string xlabel;
  std::string ylabel;
  std::string color = "Paired";
  int ncolor = 8;
  std::string gnu_cmd;
  std::string lw = "2"; // only for line_series (line width)
  bool latex_output = false;
  std::vector<std::string> legend;
  vec2 xrange = vec2::Zero();
  vec2 yrange = vec2::Zero();
};

line_spec parsing_line_spec(std::string par);

map_spec parsing_map_spec(std::string par);

fig_spec parsing_fig_spec(std::string par);

std::string find_color(std::string short_cut);

const std::map<std::string, std::string> g_line_color = line_color();
const std::map<std::string, std::string> g_line_style = line_style();
const std::map<std::string, std::string> g_line_marker = line_marker();
const std::map<std::string, std::string> g_line_dash = line_dash();

const std::map<std::string, std::string> g_map_color = map_color();
}
}
