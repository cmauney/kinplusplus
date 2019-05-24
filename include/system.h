#pragma once

#include <algorithm>
#include <array>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <type_traits>

#include <plog/Log.h>
#include <autodiff/forward.hpp>
#include <autodiff/forward/eigen.hpp>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include "arrhenius_reaction.h"
#include "configurator.h"
#include "defs.h"
#include "environment.h"
#include "io_helpers.h"
#include "network.h"
#include "reaction.h"
#include "species.h"
#include "system_constants.h"

namespace detail
{
template< class T >
struct remove_cvref {
    typedef std::remove_cv_t<std::remove_reference_t<T>> type;
};

template< class T >
using remove_cvref_t = typename remove_cvref<T>::type;

}
namespace kinpp {

struct system_t {
  species_collection_t  _species;
  network_t             _net;
  environment_t         _env;

  //idx_t         _sol_size;

  //vector_t   _x, _x0, _dxdt;
  //matrix_t   _jac;

  system_t( const std::string& network_filename )
      : _species( network_filename ),
        _net( network_filename ),
        _env( environment_t::TEST_ENV() )
  {
    _net.map_species_to_reactions(_species);
    LOGI << "system constructed from \"" << network_filename << "\"";
    //LOGI << "vector (matrix) size is " << _sol_size << "(" << _sol_size << "x"
    //     << _sol_size << ")";
  }

  ~system_t() = default;

  inline auto f( const vector_t& x, double t )
  {
    _env.update_environment(t);

    vector_t dxdt = vector_t::Zero( x.size() );
    //ctor_t k    = vector_t::Zero( x.size() );
    vector_t k = vector_t( dxdt );
    std::apply(
        [&]( auto&& rv ) {
      for ( auto&& r : rv ) {
        auto [k, rmap, pmap] = r.f( x, _env );

        scalar_t tmp = k;
        for ( auto& [idx, stoc] : rmap ) { tmp *= x[idx];}
        for ( auto& [idx, stoc] : rmap ) { dxdt[idx] -= tmp; }
        for ( auto& [idx, stoc] : pmap ) { dxdt[idx] += tmp; }
      }
        },
        _net.reactions );

    return dxdt;
  }

  // the parameter needs to be writable, b/c jacobian changes a field value
  inline auto J( vector_t& x, double t )
  {
    vector_t dxdt = this->f( x, t );
    auto fn   = [this, t]( auto&& v ) { return this->f( v, t ); };
    return autodiff::jacobian( fn, dxdt, x );
  }

  inline auto create_solution_vector(const composition_map_t& initial_composition)
  {
    vector_t x = vector_t::Zero(_species.n_species());
    for(const auto& [specname, nfrac] : initial_composition)
    {
      x[_species(specname)] = nfrac;
    }
    return x;
  }

  inline auto create_solution_vector()
  {
    return create_solution_vector({});
  }

};

}