#pragma once

#include <fstream>
#include <iostream>
#include <string>

#include <boost/lexical_cast.hpp>

#include "configurator.h"
#include "defs.h"

namespace kinpp {

template<typename Solution, typename Writer>
struct observer_t {
  const string_v snv;

  // observing_t data;
  std::vector<Solution> solutions;
  std::vector<double>   times, timesteps;

  idx_t store_every, print_every, dump_every;
  idx_t n_stores, n_prints, n_dumps;

  Writer outf;

  observer_t( idx_t    uid,
              idx_t    nstore,
              idx_t    nprint,
              idx_t    ndump,
              string_v snames )
      : snv( std::move( snames ) ),
        store_every( nstore ),
        print_every( nprint ),
        dump_every( ndump ),
        n_stores( 0 ),
        n_prints( 0 ),
        n_dumps( 0 ),
        outf( uid ) {
    solutions.reserve( ndump );
    times.reserve( ndump );
    timesteps.reserve( ndump );
  }

  inline void operator()( idx_t           step,
                          const Solution& s,
                          const double    t,
                          const double    dt ) {
    if ( step % store_every == 0 ) {
      // std::apply(
      //     [&]( const auto&... a ) {
      //       ( static_cast<void>( a.push_back( args ) ), ... );
      //     },
      //     data );
      solutions.push_back( s );
      times.push_back( t );
      timesteps.push_back( dt );
      ++n_stores;
    }

    if ( step % print_every == 0 ) {
      std::cout << "time = " << t << " s, dt = " << dt << " s\n";
      for ( idx_t i = 0; i < s.size(); ++i ) {
        std::cout << "::" << snv[i] << " = " << s[i] << " \n";
      }
      ++n_prints;
    }

    if ( step % dump_every == 0 ) {
      outf.append( solutions, times, timesteps, snv );
      solutions.clear();
      times.clear();
      timesteps.clear();
      ++n_dumps;
    }
  }
};

template<typename Solution>
struct FstreamWriter {
  std::string   filename;
  std::ofstream fobj;

  FstreamWriter( idx_t id ) {
    filename =
        "fstream_observer_" + boost::lexical_cast<std::string>( id ) + ".dat";
    initialize();
  }

  virtual ~FstreamWriter() {
    if ( fobj ) fobj.close();
  }

  virtual void initialize() {
    fobj.open( filename, std::ios::trunc );
    fobj << "itch data file \n";
    fobj.close();
  }

  inline virtual void append( const std::vector<Solution>& s,
                              std::vector<double>&         t,
                              std::vector<double>&         dt,
                              const string_v&              n ) {
    fobj.open( filename, std::ios::app );
    for ( idx_t i = 0; i < s.size(); ++i ) {
      fobj << "time = " << t[i] << "s , dt = " << dt[i] << " s\n";
      for ( idx_t ii = 0; ii < s[i].size(); ++ii )
        fobj << "::" << n[ii] << " = " << s[i][ii] << " \n";
    }
    fobj.close();
  }
};
}  // namespace itch