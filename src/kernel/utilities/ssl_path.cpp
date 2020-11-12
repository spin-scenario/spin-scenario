#include "ssl_path.h"
#include <boost/algorithm/string.hpp>

namespace ssl {
namespace utility {
// terminal dir by the user.
std::string g_terminal_dir = "";

#ifdef WIN32
std::string g_install_dir = "../.."; 
#else
std::string g_install_dir = "..";
#endif

void set_terminal_dir(std::string s) {
    g_terminal_dir = s; 
    boost::replace_all(g_terminal_dir, "\\", "\/");
}
void set_install_dir(std::string s) { 
    g_install_dir = s;  // C:\Spin-Scenario\bin\spin-scenario.exe
    boost::replace_all(g_install_dir, "\\", "\/");
#ifdef WIN32
    boost::replace_last(g_install_dir, "spin-scenario.exe", "..");
#else
    boost::replace_last(g_install_dir, "spin-scenario", "..");
#endif
}
 }  // namespace utility
}  // namespace ssl
