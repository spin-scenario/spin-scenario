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
#include "seq_block.h"

namespace ssl {
namespace seq {

// write info of seq-blocks into a specific file.
// Lua usage: ssl.write('sb.txt', sb1, sb2, ...)
void write(std::string file, sol::variadic_args va, const seq_block & /*sb*/);

// print info of seq-blocks into the terminal.
// Lua usage: ssl.print(sb1, sb2, ...)
void print(sol::variadic_args va, const seq_block & /*sb*/);

// visualize seq-blocks via Gnuplot.
// Lua usage: ssl.plot(sb1, sb2, ...)
void plot(sol::variadic_args va, const seq_block & /*sb*/);

void plot(const sol::table &t);

// calculate specgram of specific rf pulse using STFT. This is particularly useful for
// the analysis of time-freq characteristic of complex shaped pulses.
// Lua usage: ssl.specgram(rf, {wlen = 16, overlap = 0.9, nfft = 2048, style = 'dB'})
void specgram(const seq_block &rf, const sol::table &t);
// Lua usage: ssl.specgram('coop_pulse.txt', {col ='6 7', fs = 219888, wlen = 32, overlap = 0.9, nfft = 1024})
void specgram(std::string file_name, const sol::table &t);
void specgram(const cx_vec &sig, const sol::table &t, double fs, std::string label = "");

sol::object multi_shaped_rf(const sol::table &t, sol::this_state s);

// each sequence starts with a keyword seq, followed by a pair of braces, containing the seq body.
// Lua usage: seq{sb1, sb2, ...}
seq_block &serial(const sol::table &t);
seq_block &serial(std::vector<seq_block *> sbs);

// running the pulse seq, note 'sb' should be the serial glue block.
sol::object run_seq(const seq_block &sb);
sol::object run_seq_api(const sol::table &t);

// Lua usage: sb#1 ==>sb*1 ==> sb's loop priority to be 1.
seq_block &set_cycle_priority(seq_block &sb, int priority);
// Lua usage: sb~ ==>sb/1 ==> sb's loop style set to be array.
seq_block &set_loop_array(seq_block &sb, int);

seq_block &set_align(seq_block &sb, double ms);
seq_block &set_align(seq_block &sb, std::string label/*c/l/r*/);

// create a new concurrent seq-block. 
// Lua usage: sb1+sb2
seq_block &concurrent(seq_block &sb1, seq_block &sb2);

// only for gradient.
double area(const seq_block &sb);

engine *init_compute_engine(const sol::table &t);
void init_seq_param(const sol::table &t);
}
}
