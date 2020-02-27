#include "ssl_path.h"

namespace ssl {
namespace utility {

#ifdef WIN32
 std::string g_project_path = ".."; // ONLY for development purpose.
//std::string g_project_path =
    //"C:/Program Files/spin-scenario 1.0.0/bin/..";  // WINDOWS 10 INSTALL.
string g_spin_scenario = "";
#else
string g_spin_scenario = "";
// std::string g_project_path = ".."; //"../spin-scenario DEBUG";
std::string g_project_path = "/usr/bin/..";  //"../spin-scenario"; UBUNTU DEB.
#endif

void set_ssl_usr_dir(string s) { 
	g_spin_scenario = s;
}
}  // namespace utility
}  // namespace ssl