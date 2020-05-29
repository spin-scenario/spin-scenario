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

#include <kernel/equip/coil_array.h>
#include <boost/format.hpp>

namespace ssl {
namespace equip {

CoilArray::CoilArray() {
  // TODO Auto-generated constructor stub

}
CoilArray::CoilArray(const char *coil_file) {
  load_via_xml(coil_file);
}

CoilArray::~CoilArray() {
  // TODO Auto-generated destructor stub
}
std::vector<Coil *> CoilArray::get_coil_pointer(CoilMode mode) const {
  std::vector<Coil *> coils;
  for (size_t i = 0; i < array_.size(); i++) {
    if (array_[i]->get_mode() == mode || array_[i]->get_mode() == kTxRx)
      coils.push_back(array_[i]);
  }
  return coils;
}
std::vector<size_t> CoilArray::get_coil_index(CoilMode mode) const {
  std::vector<size_t> index;
  for (size_t i = 0; i < array_.size(); i++) {
    if (array_[i]->get_mode() == mode || array_[i]->get_mode() == kTxRx)
      index.push_back(i);
  }
  return index;
}
std::vector<vec3> CoilArray::get_shape() const {
  std::vector<vec3> data;
  for (size_t i = 0; i < array_.size(); i++) {
    std::vector<vec3> tmp = array_[i]->get_shape();
    std::copy(tmp.begin(), tmp.end(), std::back_inserter(data));
  }
  return data;
}
void CoilArray::graphic_view() const {
  std::string file = "coils.xyza";
  std::ofstream ofstr(file.c_str(), std::ios::binary);
  ofstr.precision(3);
  ofstr.setf(std::ios::fixed);
  std::vector<vec3> data = get_shape();
  if (data.size() == 0)
    return;
  for (size_t i = 0; i < data.size(); i++) {
    vec3 pt = data[i];
    ofstr << pt(cx) << " " << pt(cy) << " " << pt(cz) << " " << 2 << "\n";
  }
  ofstr.close();
#ifdef SSL_ENABLE_VTK
  VTKRenderer vtk;
  vtk.view_xyza(file.c_str());
#endif
}
void CoilArray::load_via_xml(const char *coil_file) {
//  tinyxml2::XMLDocument doc;
//  tinyxml2::XMLError error = doc.LoadFile(coil_file);
//  if (error != tinyxml2::XML_NO_ERROR) {
//    std::cout
//        << boost::format("%s %s %s %s.\n") % "S-S-L error: " % error
//            % "failed to load coil file" % coil_file;
//    exit(0);
//  }
//
//  // print the seq xml tree.
////  XMLPrinter printer;
////  doc.Print(&printer);
////  std::cout << printer.CStr();
//
//  tinyxml2::XMLElement* root = doc.RootElement();
//  /*if (std_string(root->Name()) != "CoilArray") {
//	std::cout
//		<< boost::format("%1% %2% %3%.\n") % "S-S-L error: " % coil_file
//			% "should use CoilArray as root node";
//	exit(0);
//  }*/
//  tinyxml2::XMLElement* child;
//  for (child = root->FirstChildElement(); child != 0;
//      child = child->NextSiblingElement()) {
//    create_coil(child);
//  }
}
void CoilArray::create_coil(/*tinyxml2::XMLElement* node*/) {
  //Coil* coil = factory_coils_.clone_coil(node);
  //if (!coil)  // the input coil does not exist in the coil-factory.
  //{
  // /*std::cout
  //  << boost::format("%1% '%2%' %3%.\n")
  //	  % "S-S-L error: non-existent coil type -" % std_string(node->Name())
  //	  % "in coil array xml file";*/
  //  exit(0);
  //}
  //coil->associate_node(node);
  //coil->assign();
  //array_.push_back(coil);
}
}
/* namespace equip */
} /* namespace ssl */
