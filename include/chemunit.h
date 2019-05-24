#pragma once

#include <fstream>
#include <iostream>
#include <utility>
#include <vector>
#include <numeric>

#include <plog/Log.h>

#include "configurator.h"
#include "defs.h"
#include "environment.h"
#include "network.h"
#include "observer.h"
#include "stepper.h"
#include "system.h"

namespace kinpp {

template<typename Stepper>
struct chemunit_t {
  idx_t  uid;
  idx_t  n_maximum_timesteps;
  double t_start, t_end;
  double t, dt;

  system_t sys;
  vector_t          x, x0;

  Stepper                                       stepper;
  //observer_t<Solution, FstreamWriter<Solution>> obs;

  chemunit_t( idx_t id, const configurator_t& config )
      : uid( id ), n_maximum_timesteps(config.n_maximum_timesteps),
        t_start( config.time_start ),
        t_end( config.time_finish ),
        t( t_start ),
        dt( config.initial_timestep ),
        sys( config.network_filename ),
        x( sys.create_solution_vector() ),
        stepper( x.size(),
                 config.ode_rho_inf,
                 config.ode_step_tolerance,
                 config.minimum_timestep,
                 config.maximum_timestep,
                 config.ode_corrector_steps ) {
        /*obs( uid,
             config.n_store_every,
             config.n_print_every,
             config.n_dump_every,
             sys.net.species.get_name_vector() ) {*/

    YAML::Node yf = YAML::LoadFile( config.abundance_filename );
    composition_map_t iacm = yf["initial_composition"].as<composition_map_t>();

    //x  = sys.create_solution_vector();
    x0 = sys.create_solution_vector(iacm);

    LOGI << "chemunit " << uid << " constructed";
  }

  void run() {
    x = x0;
    idx_t n_step = 0;

    idx_t n_dt_same = 0;
    double last_dt = std::numeric_limits<double>::max();

    while ( t < t_end ) {

      if (n_step > n_maximum_timesteps )
      {
        LOGE << "maximum timesteps exceeded";
        break;
      }

      stepper.step( sys, x, t, dt );
      t += dt;

      //obs( n_step, x, t, dt );

      // check for max timestep too low
      if ( std::abs(last_dt - dt) < sys_constants::calc_softeps )
        ++n_dt_same;
      else
        n_dt_same = 0;

      if (n_dt_same > 100)
      {
        LOGW << "timestep at maximum; consider increasing maximum timestep";
      }

      ++n_step;
    }

    LOGI << "t >= t_end [" << t << " >= " << t_end << "]; run over";
  }
};

}  // namespace itch