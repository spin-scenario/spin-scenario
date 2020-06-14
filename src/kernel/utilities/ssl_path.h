#include <string>
namespace ssl {
namespace utility {

extern std::string g_install_dir;
extern std::string g_terminal_dir;
void set_terminal_dir(std::string s);
void set_install_dir(std::string s);
}
}  // namespace ssl