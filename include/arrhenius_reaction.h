#pragma once

#include <plog/Log.h>
#include <array>
#include <iostream>
#include <utility>
#include <tuple>

#include <autodiff/forward.hpp>

#include "defs.h"
#include "environment.h"
#include "physical_constants.h"
#include "reaction.h"

//#include "n_body.h"

namespace kinpp {

struct arrhenius_parameters_t
{
  scalar_t _alpha;
  scalar_t _beta;
  scalar_t _gamma;

  arrhenius_parameters_t() :_alpha(0.0), _beta(0.0), _gamma(0.0) {}

  arrhenius_parameters_t(double a, double b, double c){
    _alpha = a;
    _beta = b;
    _gamma = c;
  }

};

struct arrhenius_reaction_t : reaction_t<arrhenius_reaction_t, arrhenius_parameters_t>
{
  using reaction_t::reaction_t;

  inline auto f( vector_cref_t concentrations,
                 const environment_t& env) const
  {
    scalar_t m = 1.0; // constructer is broken, set after
    if ( _depend.size() > 0 )
    {
      m = 0.0;
      for ( const auto& [idx, val] : _depend )
      {
        m += val * concentrations[idx];
      }
    }

    auto k = m * _parameters._alpha
          * pow(env.T, _parameters._beta)
          * exp(-_parameters._gamma * env.iT);

    return std::make_tuple(k, _react, _prod);
  }

};

using arrhenius_reaction_v = std::vector<arrhenius_reaction_t>;
//using rxntypes = decltype(std::tuple_cat(std::declval<rxntypes>(), std::make_tuple(std::declval<arrhenius_reaction_t>())));
//using arrhenius_v = std::vector<arrhenius_reaction_t>;
//using rxntypes = decltype(std::tuple_cat(std::declval<rxntypes>(), std::make_tuple(std::declval<std::vector<arrhenius_reaction_t>>())));
//APPEND_REACTION(arrhenius_reaction)
}  // namespace
