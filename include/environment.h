#pragma once

#include <cmath>
#include <functional>

namespace kinpp {
struct environment_t {
  std::function<double( double )> f_temperature;
  std::function<double( double )> f_density;

  double T, T2, T3, T4;
  double iT, iT2, iT3, iT4;
  // double expT, expmT, expiT, expmiT;

  double rho, irho;

  environment_t()                         = default;
  environment_t( const environment_t& e ) = default;
  ~environment_t()                        = default;

  inline void update_environment( double t ) {
    T   = f_temperature( t );
    rho = f_density( t );

    T2 = T * T;
    T3 = T * T2;
    T4 = T * T3;

    iT  = 1. / T;
    iT2 = iT * iT;
    iT3 = iT2 * iT;
    iT4 = iT3 * iT;

    // expT   = std::exp( T );
    // expmT  = 1. / expT;
    // expiT  = std::exp( iT );
    // expmiT = 1. / expiT;

    irho = 1.0 / rho;
  }

  static environment_t TEST_ENV() {
    environment_t t;

    t.f_temperature = []( double t ) { return 1000.0; };
    t.f_density     = []( double t ) { return 1.0; };

    t.update_environment( 0.0 );

    return t;
  };
};

}  // namespace itch
