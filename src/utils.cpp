#include <fmt/core.h>
#include <fmt/color.h>

void print_done(const std::string message = "done\n")
{
    fmt::print(fmt::fg(fmt::color::green), message);
}
