#pragma once

#include <cinttypes>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include <eigen3/Eigen/Dense>
#include <autodiff/forward.hpp>
#include <autodiff/forward/eigen.hpp>

#include "system_constants.h"

namespace kinpp {

// typedef std::map<std::string, double> composition_t;
// typedef std::vector<std::string> vector_str_t;

using idx_t = std::vector<int>::size_type;
using scalar_t   = autodiff::dual;

//using vector_t   = std::vector<scalar_t>
using vector_t      = Eigen::VectorXdual;
using vector_ref_t  = Eigen::Ref<vector_t>&;
using vector_cref_t = const Eigen::Ref<const vector_t>&;

using matrix_t      = Eigen::MatrixXd;
using matrix_ref_t  = Eigen::Ref<matrix_t>&;
using matrix_cref_t = const Eigen::Ref<const matrix_t>&;
using idx_v         = std::vector<idx_t>;

using string_v = std::vector<std::string>;
using int_v    = std::vector<int>;

using composition_map_t = std::map<std::string, double>;
using species_map_t = std::map<std::string, idx_t>;

// typedef std::vector<std::array<idx_t, sys_constants::max_reactants>>
//     reaction_map_t;

// static reaction_map_t::value_type emtpy_map_row() {
//   reaction_map_t::value_type ret;
//   for ( idx_t i = 0; i < sys_constants::max_reactants; ++i ) {
//     ret[i] = sys_constants::idx_undefined;
//   }
//   return ret;
// }

static std::string kinpp_msg() {
  std::string m;
  m += "\n";
  m += "********************************\n";
  m += "*   kinetics Plus Plus         *\n";
  m += "********************************\n";
  m += "\n";

  return m;
}
}  // namespace kinpp