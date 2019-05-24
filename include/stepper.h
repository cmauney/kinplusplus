#pragma once

#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <plog/Log.h>
#include <autodiff/forward.hpp>
#include <autodiff/forward/eigen.hpp>
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

// #include <fenv.h>
// int _feenableexcept_status =
//     feenableexcept( FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW );

namespace kinpp {

struct implicit_euler {
  double epsilon;

  Eigen::VectorXd _x, _dxdt;
  Eigen::VectorXd _b;
  Eigen::MatrixXd _jac;
  Eigen::MatrixXd _id;

  Eigen::VectorXd _xp1;

  double en2, en1, en;

  // to be parameterized later
  static constexpr double step_tol  = 1.0E-6;
  static constexpr double min_dt    = 1.0E-12;
  static constexpr double max_dt    = 1.0E0;
  static constexpr idx_t  max_steps = 10000;
  static constexpr double ekP = 0.075, ekI = 0.175, ekD = 0.01;

  implicit_euler( idx_t N, double _epsilon )
      : epsilon( _epsilon ),
        _x( N ),
        _dxdt( N ),
        _b( N ),
        _jac( N, N ),
        _id( Eigen::MatrixXd::Identity( N, N ) ),
        _xp1( N ) {
    en2 = step_tol;
    en1 = step_tol;
    en  = 0;
  }

  template<class System>
  void step( System& s, Eigen::VectorXd& x, double t, double dt ) {
    t += dt;
    s.f( x, _dxdt, t );

    _b = dt * _dxdt;

    // s.J( x, _dxdt, _jac );
    s.J( _jac );

    Eigen::MatrixXd A  = dt * _jac - _id;
    Eigen::VectorXd nx = A.partialPivLu().solve( _b );

    _x = x - nx;

    while ( nx.norm() > epsilon ) {
      s.f( _x, _dxdt, t );
      _b = x - _x + dt * _dxdt;
      nx = A.partialPivLu().solve( _b );

      _x -= nx;
    }
    x = _x;
  }
};

struct generalized_alpha {
  // integration parameters
  double RTOL, MIN_DT, MAX_DT, PCI_MAX;

  // vector storage
  vector_t ydot, ydot0;
  vector_t yp1, ydotp1;
  vector_t yp_alphaf, ydotp_alpham, del_yp_alphaf;
  vector_t g;
  vector_t y_backeuler;
  // matrix storage
  Eigen::MatrixXd K, Jacobin;
  Eigen::MatrixXd II, II_K;

  // fp storage
  double rho_inf;
  double alpha_m, alpha_f;
  double gamma, gammai, gamma_m;
  double idt, err_p1;
  double mf_ratio, mf_ratio_gi, mf_ratio_gi_idt;
  double alpha_m_gi;

  // last iteration storage
  double new_dt;
  double last_dt;
  double last_err;

  generalized_alpha( idx_t  N,
                     double rinf,
                     double step_tolerance,
                     double dt_min,
                     double dt_max,
                     idx_t  max_pci_steps )
      : RTOL( step_tolerance ),
        MIN_DT( dt_min ),
        MAX_DT( dt_max ),
        PCI_MAX( max_pci_steps ),
        ydot( vector_t::Zero( N ) ),
        ydot0( vector_t::Zero( N ) ),
        yp1( vector_t::Zero( N ) ),
        ydotp1( vector_t::Zero( N ) ),
        yp_alphaf( vector_t::Zero( N ) ),
        ydotp_alpham( vector_t::Zero( N ) ),
        del_yp_alphaf( vector_t::Zero( N ) ),
        g( vector_t::Zero( N ) ),
        y_backeuler( vector_t::Zero( N ) ),
        K( matrix_t::Zero( N, N ) ),
        Jacobin( matrix_t::Zero(N, N) ),
        II( matrix_t::Identity( N, N ) ),
        II_K( matrix_t::Identity( N, N ) ),
        rho_inf( rinf ),
        last_dt( MIN_DT ),
        last_err( RTOL ) {
    alpha_m = 0.5 * ( 3. - rho_inf ) / ( rho_inf + 1. );
    alpha_f = 1.0 / ( rho_inf + 1. );
    gamma   = 0.5 + alpha_m - alpha_f;
    gammai  = 1. / gamma;
    gamma_m = ( gamma - 1. ) / gamma;

    mf_ratio    = alpha_m / alpha_f;
    mf_ratio_gi = mf_ratio * gammai;
    alpha_m_gi  = alpha_m * gammai;
  }

  template<class System>
  inline void step( System& s, vector_t& y, double t, double& dt ) {
    ydot = vector_t::Zero( y.size() );
    // get initial y_dot (can be optimized out?)
    ydot0 = s.f( y, t );

    // initialize error
    err_p1 = std::numeric_limits<double>::max();
    // iteration counter
    idx_t dt_step = 0;
    while ( err_p1 > RTOL ) {
      // set dt-dependent values
      idt             = 1. / dt;
      mf_ratio_gi_idt = mf_ratio_gi * idt;
      II_K            = II * mf_ratio_gi_idt;

      // copy initial vectors
      yp1    = y;
      ydotp1 = gamma_m * ydot0;
      // eq 31,32
      yp_alphaf    = y + alpha_f * ( yp1 - y );
      ydotp_alpham = ydot0 + alpha_m * ( ydotp1 - ydot );
      for ( idx_t i = 0; i < PCI_MAX; ++i ) {
        g = s.f( yp_alphaf, t + alpha_f * dt );
        Jacobin = s.J( yp_alphaf, t + alpha_f * dt );
        K = ( II_K - Jacobin );

        g -= ydotp_alpham;

        del_yp_alphaf = K.fullPivLu().solve( g );

        yp_alphaf += del_yp_alphaf;
        ydotp_alpham = ( ( 1. - alpha_m_gi ) * ydot0 ) +
                       ( mf_ratio_gi_idt * ( yp_alphaf - y ) );
      }

      yp1    = y + ( yp_alphaf - y ) / alpha_f;
      ydotp1 = ydot0 + ( ydotp_alpham - ydot0 ) / alpha_m;

      // timestep correction
      y_backeuler = y + dt * ydotp1;

      last_dt = dt;
      err_p1  = autodiff::val(( yp1 - y_backeuler ).norm());
      // constraint - nothing below zero (with tol)
      bool ok = true;

      if ( !ok ) {
        dt *= 0.5;
      }  // otherwise use elementary stepper
      else {
        if ( err_p1 < 1.0E-10 )
          dt *= 2.0;
        else
          dt = std::sqrt( 0.8 * RTOL / ( err_p1 ) ) * dt;

        dt = std::fmin( MAX_DT, std::fmax( MIN_DT, dt ) );
      }

      ++dt_step;
    }
    y = yp1;
  }

};  // namespace itch

}  // namespace itch
