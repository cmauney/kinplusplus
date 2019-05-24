#pragma once

#include <string>

namespace kinpp {

namespace sys_constants {
constexpr size_t max_reactants     = 8;
constexpr size_t max_effenciencies = 16;
constexpr size_t idx_notfound      = -1;
constexpr size_t idx_undefined     = -2;
constexpr double calc_hardeps      = 1.0E-15;
constexpr double calc_softeps      = 1.0E-8;

std::string headbar = std::string( 80, '-' );

}  // namespace sys_constants

}  // namespace kinpp
