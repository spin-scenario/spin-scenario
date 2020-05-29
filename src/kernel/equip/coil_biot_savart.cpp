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

#include <kernel/equip/coil_biot_savart.h>
//#include <boost/math/special_functions/ellint_1.hpp>
//#include <boost/math/special_functions/ellint_2.hpp>

namespace ssl {
namespace equip {

CoilBiotSavart::CoilBiotSavart()
    : radius_(0.06),
      azimuth_(0),
      elevation_(0),
      segment_(20),
      current_direction_(1) {
  // TODO Auto-generated constructor stub

}

CoilBiotSavart::~CoilBiotSavart() {
  // TODO Auto-generated destructor stub
}
void CoilBiotSavart::assign() {
  Coil::assign();
  //std::string str = getAttribute("azimuth");
  /*if (!str.empty())
	azimuth_ = boost::lexical_cast<double>(str);

  str = getAttribute("elevation");
  if (!str.empty())
	elevation_ = boost::lexical_cast<double>(str);*/

  azimuth_ *= _pi / 180;
  elevation_ *= _pi / 180;

  //str = getAttribute("pos");
  //if (!str.empty()) {
  // std::vector<double> pars;// = str_cast<double>(str);
  //  pos_(0) = pars[0];
  //  pos_(1) = pars[1];
  //  pos_(2) = pars[2];
  //}
  segments_ = get_sampe_points();
}
cd CoilBiotSavart::get_sensitivity(const vec3 &pos) const {
  cd B1s;
  vec3 B1 = vec3::Zero();
  //for (size_t i = 0; i < segments_.size() - 1; i++) {
  //   vec3 dl = segments_[i + 1] - segments_[i];  // components of segment std::vector dl.
  //   vec3 cp = (segments_[i + 1] + segments_[i]) / 2;  // the location of the midpoint of current segment.
  //   vec r = pos - cp;  // segment on loop and observation point.
  //  double r3 = pow(norm(r, 2), 3);  // r^3 from r std::vector.
  //   vec3 dl_r = cross(dl, r);  // cross product between dl and r.
  //  B1 += (radius_ / k_2pi / r3) * dl_r;  // increment sum of B1 magnetic field.
  //}
  double magn = sqrt(pow(B1(cx), 2) + pow(B1(cy), 2));
  double phase = atan2(B1(cy), B1(cx));
  return cd(magn, phase);
  //double mu0 = 4 * _pi * 1e-7;  // permeability of free space (T*m/A).
  //double I = 1;  // current in the loop (A).
  //double constant = mu0 / (4 * _pi) * I;   // Useful constant.
  // the magnetic field at the center of the coil B0 = Imu0/2*radius.
  // ref from http://www.netdenizen.com/emagnettest/offaxis/?offaxisloop
}
std::vector<vec3> CoilBiotSavart::get_sampe_points() const {
  // vec theta =  linspace< vec>(
  //    0, k_2pi * (double) current_direction_, segment_);
  // rowvec3 N;
  //N << cos(elevation_) * cos(azimuth_) << cos(elevation_) * sin(azimuth_)
  //  << sin(elevation_) <<  endr;
  //N /= norm(N);
  // mat V = null(N);
  //double ang = atan2(dot(cross(V.col(0), V.col(1)), N),
  //                   dot(cross(V.col(0), N), cross(V.col(1), N)));  // determine angle direction.
  //V.col(0) = (ang / abs(ang)) * V.col(0);  // match angle direction.
  std::vector<vec3> points;
  /* for (size_t i = 0; i < theta.size(); i++) {
	  vec3 p = pos_
		 + radius_ * (cos(theta(i)) * V.col(0) + sin(theta(i)) * V.col(1));
	 points.push_back(p);
   }*/
  return points;
}
} /* namespace equip */
} /* namespace ssl */
