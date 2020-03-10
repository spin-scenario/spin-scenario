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

#include "ssl_plot.h"

namespace ssl {
namespace utility {

string g_output_terminal = "qt";
string g_output_font = "Arial,12";
void set_output_terminal(const sol::table &t) {
  g_output_terminal = retrieve_table_str("type", t, "gnuplot terminal");
    if(is_retrievable("font", t))
        g_output_font = retrieve_table_str("font", t, "gnuplot terminal");
}
string terminal_cmd(string key) {
	string s;
	if (key == "qt")
    s ="set terminal qt enhanced font '"+ g_output_font + "'\n";
	if (key == "png")
    s= "set terminal pngcairo enhanced font '" + g_output_font + "'\n";
	if (key == "eps")
    s = "set terminal postscript eps enhanced color font '" + g_output_font + "'\n";
    if (key == "svg")
		s = "set terminal svg fname 'Verdana' fsize 10\n";
	if (key == "tex")
    s = "set terminal epslatex color colortext header '\\newcommand{ \\ft }[0]{ \\footnotesize}'\n";
  return s;
}
  vec linspace(double start, double stop, int num) {
  return vec::LinSpaced(num, start, stop);
  //return vec::LinSpaced(num, start, stop).cwiseInverse();
}

sol::object linspace_lua(double start, double stop, int num) {
  vec val = vec::LinSpaced(num, start, stop);
  sol::table t = g_lua->create_table();
  for (size_t i = 0; i < val.size(); i++)
    t.add(val[i]);
  return t;
}

mat random(int rows, int cols) {
  return mat::Random(rows, cols);
}

mat vec_table(const sol::table &t) {
  int cols = t.size();
  sol::object val = t[1];
  int rows = val.as<vec>().size();
  mat m(rows, cols);
  for (int i = 0; i < cols; i++) {
    sol::object val = t[i + 1];
    m.col(i) = val.as<vec>();
  }
  return m;
}

void plot(sol::variadic_args va, const line &b) {
  plot("", va, b);
}

void plot(string fig_info, sol::variadic_args va, const line &) {
#ifdef GP_SCRIPT_OUTPUT
  ofstream gp("plot_line.gnu");
#else
  Gnuplot gp;
#endif
  gp << "reset\n";
  gp << "load '" << g_project_path << "/share/spin-scenario/config/gnuplot/xyborder.cfg'\n";
  gp << "load '" << g_project_path << "/share/spin-scenario/config/gnuplot/grid.cfg'\n";

  fig_spec fig = parsing_fig_spec(fig_info);
  gp << "set title '" << fig.title << "'\n";
  gp << "set xlabel  '" << fig.xlabel << "'\n";
  gp << "set ylabel  '" << fig.ylabel << "'\n";

  gp << "set border back\n";
  gp << "set key opaque\n";

  gp << "set key width 1\n";

  gp << "load '" << g_project_path << "/share/spin-scenario/config/gnuplot/colorbrewer/" << find_color(fig.color) << "'\n";

    string time_s = sys_time();
  gp << terminal_cmd(g_output_terminal);
  if (g_output_terminal != "qt")
    gp << "set output "
       << "'output_" << time_s << "." << g_output_terminal << "'\n";

  size_t nleg = fig.legend.size();
  if (fig.xrange.norm())
    gp << "set xrange [" << fig.xrange[0] << ":" << fig.xrange[1] << "]\n";
  if (fig.yrange.norm())
    gp << "set yrange [" << fig.yrange[0] << ":" << fig.yrange[1] << "]\n";

  gp << fig.gnu_cmd << "\n";

  vector<string> files;
  vector<line_spec> lines;
  size_t i = 0;
  //g_lua->script("os.execute('rm -rf gnuplot')");
  g_lua->script("os.execute('mkdir gnuplot')");
  for (auto v : va) {
    const line &val = v;
    string file;
    if (val.is_file) {
      file = val.file;
      g_lua->script("os.execute('cp " + file + " gnuplot')");
    }
    else {
      file = "dat" + boost::lexical_cast<string>(i) + "_" + time_s;
      if (val.is_x_only)
        write_line(val.x, "gnuplot/" + file);
      else
        write_line(val.x, val.y, "gnuplot/" + file);
    }

    files.push_back(file);
    i++;

    line_spec line = parsing_line_spec(val.line_spec);
    gp << "set style line " << i << " lw "<<fig.lw<<"\n";
    if (!line.dash.empty())
      gp << "set style line " << i << " dt " << line.dash << "\n";
    if (!line.marker.empty())
      gp << "set style line " << i << " pt " << line.marker << " ps 1\n";
    if (!line.color.empty())
      gp << "set style line " << i << " lc rgb " << line.color << "\n";
    lines.push_back(line);
  }
  if (!files.size())
    return;

  gp << "cd 'gnuplot'\n";
  gp << "plot '" << files[0] << "' with " << lines[0].style << " ls 1";
  if (nleg)
    gp << " t '" << fig.legend[0] << "'";
  if (files.size() == 1) {
    gp << "\n";
    return;
  } else
    gp << ",";

  for (i = 1; i < files.size() - 1; i++) {
    gp << " '" << files[i] << "' with " << lines[i].style << " ls " << (i + 1);
    if (nleg) {
      string leg = (i <= nleg - 1) ? fig.legend[i] : fig.legend.back();
      gp << " t '" << leg << "',";
    } else
      gp << ",";
  }

  gp << " '" << files.back() << "' with " << lines.back().style << " ls " << files.size();
  if (nleg) {
    string leg = (i <= nleg - 1) ? fig.legend[i] : fig.legend.back();
    gp << " t '" << leg << "'\n";
  } else
    gp << "\n";
  gp << "set output\n";

#ifdef GP_SCRIPT_OUTPUT
  gp.close();
#endif
}

void plot(const line_series &v) {
  plot("", v);
}

void plot(string fig_info, const line_series &v) {
#ifdef GP_SCRIPT_OUTPUT
  ofstream gp("plot_line.gnu");
#else
  Gnuplot gp;
#endif
  gp << "reset\n";
  gp << "load '" << g_project_path
     << "/share/spin-scenario/config/gnuplot/xyborder.cfg'\n";
  gp << "load '" << g_project_path
     << "/share/spin-scenario/config/gnuplot/grid.cfg'\n";

  fig_spec fig = parsing_fig_spec(fig_info);
  gp << "set title '" << fig.title << "'\n";
  gp << "set xlabel  '" << fig.xlabel << "'\n";
  gp << "set ylabel  '" << fig.ylabel << "'\n";
  gp << "set key width 2\n";
  gp << "set key samplen 2\n";

  int ncolor = fig.ncolor;
  gp << "load '" << g_project_path
     << "/share/spin-scenario/config/gnuplot/colorbrewer/"
     << find_color(fig.color) << "'\n";

  string time_s = sys_time();
  gp << terminal_cmd(g_output_terminal);
  if (g_output_terminal != "qt")
    gp << "set output "
       << "'output_" << time_s << "." << g_output_terminal << "'\n";

  size_t nleg = fig.legend.size();
  if (fig.xrange.norm())
    gp << "set xrange [" << fig.xrange[0] << ":" << fig.xrange[1] << "]\n";
  if (fig.yrange.norm())
    gp << "set yrange [" << fig.yrange[0] << ":" << fig.yrange[1] << "]\n";

  gp << "set border back\n";
  gp << "set key opaque\n";
  gp << fig.gnu_cmd << "\n";

  vector<string> files;
  size_t i = 0;
  g_lua->script("os.execute('mkdir gnuplot')");
  if (v.is_file)
    files = v.files;
  else {
    for (i = 0; i < v.y.size(); i++) {
      string file_i = "dat" + boost::lexical_cast<string>(i + 1) + "_" + time_s;
      if (v.is_y_only)
        write_line(v.y[i], "gnuplot/" + file_i);
      else
         write_line(v.x, v.y[i], "gnuplot/"  + file_i);
      files.push_back(file_i);
    }
  }

  if (!files.size()) return;

  gp << "cd 'gnuplot'\n";
  gp << "plot '" << files[0] << "' with l ls 1 lw "<<fig.lw;
  if (nleg) gp << " t '" << fig.legend[0] << "'";
  if (files.size() == 1) {
    gp << "\n";
    return;
  } else
    gp << ",";

  for (i = 1; i < files.size() - 1; i++) {
    gp << " '" << files[i] << "' with l ls "
       << ((i + 1) % ncolor == 0 ? ncolor : (i + 1) % ncolor) << " lw "<<fig.lw;
    if (nleg) {
      string leg = (i <= nleg - 1) ? fig.legend[i] : fig.legend.back();
      gp << " t '" << leg << "',";
    } else
      gp << ",";
  }

  gp << " '" << files.back() << "' with l ls "
     << ((i + 1) % ncolor == 0 ? ncolor : (i + 1) % ncolor) << " lw " << fig.lw;
  if (nleg) {
    string leg = (i <= nleg - 1) ? fig.legend[i] : fig.legend.back();
    gp << " t '" << leg << "'\n";
  } else
    gp << "\n";

#ifdef GP_SCRIPT_OUTPUT
  gp.close();
#endif
}

void plot(sol::variadic_args va, const map &b) {
  plot("", va, b);
}

void plot(string fig_info, sol::variadic_args va, const map &) {
#ifdef GP_SCRIPT_OUTPUT
  ofstream gp("plot_map.gnu");
#else
  Gnuplot gp;
#endif
  gp << "reset\n";
  gp << "load '" << g_project_path << "/share/spin-scenario/config/gnuplot/xyborder.cfg'\n";
  gp << "load '" << g_project_path << "/share/spin-scenario/config/gnuplot/grid.cfg'\n";
  gp << "unset key\n";

  fig_spec fig = parsing_fig_spec(fig_info);
  gp << "set title '" << fig.title << "'\n";
  gp << "set xlabel  '" << fig.xlabel << "'\n";
  gp << "set ylabel  '" << fig.ylabel << "'\n";

  gp << "load '" << g_project_path << "/share/spin-scenario/config/gnuplot/colorbrewer/" << find_color(fig.color) << "'\n";
  //gp << "set palette negative\n";

    string time_s = sys_time();
  gp << terminal_cmd(g_output_terminal);
  if (g_output_terminal != "qt")
    gp << "set output "
       << "'output_" << time_s << "." << g_output_terminal << "'\n";

  size_t i = 0;
  //g_lua->script("os.execute('rm -rf gnuplot')");
  g_lua->script("os.execute('mkdir gnuplot')");
  for (auto v : va) {
    const map &val = v;
    string file;
    if (val.is_file) {
      file = val.file;
      g_lua->script("os.execute('cp " + file + " gnuplot')");
    }
    else {
      file = "gnu_map_2d_" + boost::lexical_cast<string>(i) + "_" + time_s;
      write_map(val.m, "gnuplot/" + file);
    }

    gp << "cd 'gnuplot'\n";
    map_spec map = parsing_map_spec(val.map_spec);
        string color;
        std::map<string, string>::const_iterator iter;
        iter = g_map_color.find(map.color);
        if (iter != g_map_color.end())
            color = iter->second;

    if (map.style == "3d") {
      gp << R"(
		set pm3d at bs
        set style line 100 lt 5 lw 0.5
		set isosamples 50
		set hidden3d
	)";
    } else if (map.style == "image")
      gp << "set pm3d map\n";
    else if (map.style == "contour") {
      gp << R"(
      set hidden3d
	  set contour 
      #unset surface
      #set view map
      #set pm3d at b
      set style increment default
      set cntrparam levels auto 15 unsorted
      set dgrid3d 128, 128
	)";
    }

    gp << "set xrange [" << val.xrange[0] << ":" << val.xrange[1] << "]\n";
    gp << "set yrange [" << val.yrange[0] << ":" << val.yrange[1] << "]\n";

    if (fig.xrange.norm())
      gp << "set xrange [" << fig.xrange[0] << ":" << fig.xrange[1] << "]\n";
    if (fig.yrange.norm())
      gp << "set yrange [" << fig.yrange[0] << ":" << fig.yrange[1] << "]\n";

    gp << fig.gnu_cmd << "\n";
    //gp << "set palette negative\n";



    if (map.style == "image")
      gp << "plot '" << file << "' u (" << val.xrange[0] << "+$1*"
         << (val.xrange[1] - val.xrange[0]) / (double) (val.m.cols() - 1) << "):" << "(" << val.yrange[0]
         << " + $2*" << (val.yrange[1] - val.yrange[0]) / (double) (val.m.rows() - 1)
         << ") : ($3) matrix with image\n";
    else if (map.style == "3d")
      gp << "splot '" << file << "' u (" << val.xrange[0] << "+$1*"
         << (val.xrange[1] - val.xrange[0]) / (double) (val.m.cols() - 1) << "):" << "(" << val.yrange[0]
         << " + $2*" << (val.yrange[1] - val.yrange[0]) / (double) (val.m.rows() - 1)
         << ") : ($3) matrix with lines\n";
    else if (map.style == "contour")
      gp << "splot '" << file << "' u (" << val.xrange[0] << "+$1*"
         << (val.xrange[1] - val.xrange[0]) / (double)(val.m.cols() - 1) << "):"
         << "(" << val.yrange[0] << " + $2*"
         << (val.yrange[1] - val.yrange[0]) / (double)(val.m.rows() - 1)
         << ") : ($3) matrix with lines\n";
    gp << "set output\n";
  }
#ifdef GP_SCRIPT_OUTPUT
  gp.close();
#endif
}

string find_color(string short_cut) {
  string color;
  std::map<string, string>::const_iterator iter;
  iter = g_map_color.find(short_cut);
  if (iter != g_map_color.end())
    color = iter->second;
  else
    color = "custom/" + short_cut + ".plt";
  return color;
}

fig_spec parsing_fig_spec(string par) {
  fig_spec fig;

  boost::cmatch what;
  boost::regex reg_title("title<(.*?)>{1}");
  if (boost::regex_search(par.c_str(), what, reg_title))
    fig.title = string(what[1]);

  boost::regex reg_xlabel("xlabel<(.*?)>{1}");
  if (boost::regex_search(par.c_str(), what, reg_xlabel))
    fig.xlabel = string(what[1]);

  boost::regex reg_ylabel("ylabel<(.*?)>{1}");
  if (boost::regex_search(par.c_str(), what, reg_ylabel))
    fig.ylabel = string(what[1]);

  boost::regex reg_lw("lw<(.*?)>{1}");
  if (boost::regex_search(par.c_str(), what, reg_lw))
      fig.lw = string(what[1]);

  boost::regex reg_legend("legend<(.*?)>{1}");
  if (boost::regex_search(par.c_str(), what, reg_legend)) {
    string s = string(what[1]);
    boost::split(fig.legend, s, boost::is_any_of(";"), boost::token_compress_on);
  }

  boost::regex reg_color("color<(.*?)>{1}");
  if (boost::regex_search(par.c_str(), what, reg_color)) {
    string s = string(what[1]);
    boost::erase_first(s, " ");
    boost::erase_last(s, " ");
    vector<string> str_vec;
    boost::split(str_vec, s, boost::is_any_of(","), boost::token_compress_on);
    if (str_vec.size() == 1)
      fig.color = s;
    if (str_vec.size() == 2) {
      fig.color = str_vec[0];
      fig.ncolor = boost::lexical_cast<int>(str_vec[1]);
    }
  }

  boost::regex reg_latex("latex<(.*?)>{1}");
  if (boost::regex_search(par.c_str(), what, reg_latex)) {
    string s = string(what[1]);
    boost::erase_first(s, " ");
    boost::erase_last(s, " ");
    boost::to_upper(s);
    if (s == "ON")
      fig.latex_output = true;
    if (s == "OFF")
      fig.latex_output = false;
  }

  boost::regex reg_xrange("xrange<(\\-?\\d+\\.?\\d*):(\\-?\\d+\\.?\\d*)>{1}");
  if (boost::regex_search(par.c_str(), what, reg_xrange)) {
    fig.xrange[0] = boost::lexical_cast<double>(what[1]);
    fig.xrange[1] = boost::lexical_cast<double>(what[2]);
  }

  boost::regex reg_yrange("yrange<(\\-?\\d+\\.?\\d*):(\\-?\\d+\\.?\\d*)>{1}");
  if (boost::regex_search(par.c_str(), what, reg_yrange)) {
    fig.yrange[0] = boost::lexical_cast<double>(what[1]);
    fig.yrange[1] = boost::lexical_cast<double>(what[2]);
  }

  boost::regex reg_cmd("gnuplot<(.*?)>{1}");
  if (boost::regex_search(par.c_str(), what, reg_cmd))
    fig.gnu_cmd = string(what[1]);

  return fig;
}

line_spec parsing_line_spec(string par) {
  line_spec line;
  boost::to_lower(par);
  vector<string> par_vec;
  boost::split(par_vec, par, boost::is_any_of("\t "), boost::token_compress_on);
  std::map<string, string>::const_iterator iter;
  string line_color = "'#0060ad'"; // default line_color.
  string line_style = "lines"; // default line_style.
  string line_marker = "7"; // default line_marker.
  string line_dash = "1"; // default line_marker.
  for (size_t i = 0; i < par_vec.size(); i++) {
    iter = g_line_color.find(par_vec[i]);
    if (iter != g_line_color.end())
      line.color = iter->second;

    iter = g_line_style.find(par_vec[i]);
    if (iter != g_line_style.end())
      line.style = iter->second;

    iter = g_line_marker.find(par_vec[i]);
    if (iter != g_line_marker.end())
      line.marker = iter->second;

    iter = g_line_dash.find(par_vec[i]);
    if (iter != g_line_dash.end())
      line.dash = iter->second;
  }
  return line;
}

map_spec parsing_map_spec(string par) {
  map_spec map;

  boost::cmatch what;

  boost::regex reg_style("style<(.*?)>{1}");
  if (boost::regex_search(par.c_str(), what, reg_style))
    map.style = string(what[1]);

  return map;
}

void write_line(const vec &x, string file) {
  ofstream ofstr(file.c_str());
  ofstr.precision(4);
  ofstr << x;
  ofstr.close();
}

void write_line(const vec &x, const vec &y, string file) {
  if (x.size() != y.size())
    return;
  ofstream ofstr(file.c_str());
  ofstr.precision(4);
  mat xy(x.size(), 2);
  xy.col(0) = x;
  xy.col(1) = y;
  ofstr << xy;
  ofstr.close();
}

void write_map(const mat &x, string file) {
  ofstream ofstr(file.c_str());
  ofstr.precision(4);
  ofstr << x;
  ofstr.close();
}

std::map<string, string> line_color() {
  std::map<string, string> color_map;
  color_map.insert(pair<string, string>("r", "'red'"));
  color_map.insert(pair<string, string>("g", "'green'"));
  //color_map.insert(pair<string, string>("o", "'orange'"));
  color_map.insert(pair<string, string>("y", "'yellow'"));
  color_map.insert(pair<string, string>("b", "'blue'"));
  color_map.insert(pair<string, string>("w", "'white'"));
  color_map.insert(pair<string, string>("k", "'black'"));
  color_map.insert(pair<string, string>("c", "'cyan'"));
  return color_map;
}

std::map<string, string> line_style() {
  std::map<string, string> style_map;
  style_map.insert(pair<string, string>("lp", "linespoints"));
  style_map.insert(pair<string, string>("l", "lines"));
  style_map.insert(pair<string, string>("p", "points"));
  style_map.insert(pair<string, string>("s", "steps"));
  style_map.insert(pair<string, string>("d", "dots"));
  style_map.insert(pair<string, string>("ip", "impulses"));
  style_map.insert(pair<string, string>("hs", "histeps"));
  style_map.insert(pair<string, string>("fs", "fsteps"));
  style_map.insert(pair<string, string>("fc", "filledcurves"));
  style_map.insert(pair<string, string>("hg", "histograms"));
  style_map.insert(pair<string, string>("box", "boxes"));
  return style_map;
}

std::map<string, string> line_marker() {
  std::map<string, string> line_marker_map;
  line_marker_map.insert(pair<string, string>("+", "1"));
  line_marker_map.insert(pair<string, string>("x", "2"));
  line_marker_map.insert(pair<string, string>("*", "3"));
  line_marker_map.insert(pair<string, string>("@-", "4"));
  line_marker_map.insert(pair<string, string>("@", "5"));
  line_marker_map.insert(pair<string, string>("o-", "6"));
  line_marker_map.insert(pair<string, string>("o", "7"));
  line_marker_map.insert(pair<string, string>("^-", "8"));
  line_marker_map.insert(pair<string, string>("^", "9"));
  line_marker_map.insert(pair<string, string>("v-", "10"));
  line_marker_map.insert(pair<string, string>("v", "11"));
  line_marker_map.insert(pair<string, string>("#-", "12"));
  line_marker_map.insert(pair<string, string>("#", "13"));
  return line_marker_map;
}

std::map<string, string> line_dash() {
  std::map<string, string> line_dash_map;
  line_dash_map.insert(pair<string, string>(".", "'.'"));
  line_dash_map.insert(pair<string, string>("-", "'-'"));
  line_dash_map.insert(pair<string, string>("._", "'._'"));
  line_dash_map.insert(pair<string, string>("..-", "'..-'"));
  return line_dash_map;
}

std::map<string, string> map_color() {
  std::map<string, string> color_map;
  // diverging.
  color_map.insert(pair<string, string>("BrBG", "diverging/BrBG.plt"));
  color_map.insert(pair<string, string>("PiYG", "diverging/PiYG.plt"));
  color_map.insert(pair<string, string>("PRGn", "diverging/PRGn.plt"));
  color_map.insert(pair<string, string>("PuOr", "diverging/PuOr.plt"));
  color_map.insert(pair<string, string>("RdBu", "diverging/RdBu.plt"));
  color_map.insert(pair<string, string>("RdGy", "diverging/RdGy.plt"));
  color_map.insert(pair<string, string>("RdYlBu", "diverging/RdYlBu.plt"));
  color_map.insert(pair<string, string>("RdYlGn", "diverging/RdYlGn.plt"));
  color_map.insert(pair<string, string>("Spectral", "diverging/Spectral.plt"));

  // qualitative.
  color_map.insert(pair<string, string>("Accent", "qualitative/Accent.plt"));
  color_map.insert(pair<string, string>("Dark2", "qualitative/Dark2.plt"));
  color_map.insert(pair<string, string>("Paired", "qualitative/Paired.plt"));
  color_map.insert(pair<string, string>("Pastel1", "qualitative/Pastel1.plt"));
  color_map.insert(pair<string, string>("Pastel2", "qualitative/Pastel2.plt"));
  color_map.insert(pair<string, string>("Set1", "qualitative/Set1.plt"));
  color_map.insert(pair<string, string>("Set2", "qualitative/Set2.plt"));
  color_map.insert(pair<string, string>("Set3", "qualitative/Set3.plt"));

  // sequential.
  color_map.insert(pair<string, string>("Blues", "sequential/Blues.plt"));
  color_map.insert(pair<string, string>("BuGn", "sequential/BuGn.plt"));
  color_map.insert(pair<string, string>("BuPu", "sequential/BuPu.plt"));
  color_map.insert(pair<string, string>("GnBu", "sequential/GnBu.plt"));
  color_map.insert(pair<string, string>("Greens", "sequential/Greens.plt"));
  color_map.insert(pair<string, string>("Greys", "sequential/Greys.plt"));
  color_map.insert(pair<string, string>("Oranges", "sequential/Oranges.plt"));
  color_map.insert(pair<string, string>("OrRd", "sequential/OrRd.plt"));

  color_map.insert(pair<string, string>("PuBu", "sequential/PuBu.plt"));
  color_map.insert(pair<string, string>("PuBuGn", "sequential/PuBuGn.plt"));
  color_map.insert(pair<string, string>("PuRd", "sequential/PuRd.plt"));
  color_map.insert(pair<string, string>("Purples", "sequential/Purples.plt"));
  color_map.insert(pair<string, string>("RdPu", "sequential/RdPu.plt"));
  color_map.insert(pair<string, string>("Reds", "sequential/Reds.plt"));
  color_map.insert(pair<string, string>("YlGn", "sequential/YlGn.plt"));
  color_map.insert(pair<string, string>("YlGnBu", "sequential/YlGnBu.plt"));
  color_map.insert(pair<string, string>("YlOrBr", "sequential/YlOrBr.plt"));
  color_map.insert(pair<string, string>("YlOrRd", "sequential/YlOrRd.plt"));

  // usr def.
  color_map.insert(pair<string, string>("YiZhang16", "../usr/YiZhang16.plt"));

  return color_map;
}

}  // namespace utility
}  // namespace ssl
