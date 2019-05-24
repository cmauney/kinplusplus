
#include <plog/Log.h>
#include <chrono>
#include <iostream>
#include <string>

#include "chemunit.h"
#include "configurator.h"
#include "system.h"

const std::string config_file = "data/configs/test_ode.cfg";

int main( int argc, char* argv[] ) {
  plog::init( plog::verbose, "log.txt", 0, 1 );

  LOGI << "kin++ logfile";

  // itch::spill_default( config_file );

  kinpp::configurator_t config;
  config.load( config_file );

  kinpp::chemunit_t<kinpp::generalized_alpha> cu(0, config);

  cu.run();

  auto xan = exp( -1.0E-9 * cu.t );
  std ::cout << cu.x[0] << ", " << xan << std::endl;
  std ::cout << cu.x[1] << ", " << 1.0 - xan << std::endl;
  std ::cout << cu.x0 << std::endl;
  // std::cout << cu.sys._net.print_network() << std::endl;
  //   std::cout << cu.sys._net.print_network() << std::endl;
  //   std ::cout << tt << std::endl;
  // itch::chemunit_t<Eigen::VectorXd, itch::generalized_alpha> cu( 0,
  // config
  // );
  // // itch::solver_t solv( config );

  // // auto x = solv.create_solution_vectors<Eigen::VectorXd>( config );

  // // double t     = config.runtime_cfg.time_start;
  // // double t_end = config.runtime_cfg.time_finish;

  // // double dt = config.ode_cfg.initial_timestep;

  // // itch::generalized_alpha stepper( x.size(), config );

  // auto t1 = std::chrono::high_resolution_clock::now();

  // // while ( t < t_end ) {
  // //   stepper.step( solv, x, t, dt );

  // //   t += dt;
  // cu.run();

  // auto t2 = std::chrono::high_resolution_clock::now();

  // auto duration =
  //     std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1
  //     ).count();
  // std::cout << "ttf = " << duration << " us\n";

  // for ( itch::idx_t i = 0; i < x.size(); ++i )
  //   std ::cout << solv.net.species.get_species_name( i ) << " = " <<
  //   x[i]
  //              << std::endl;
  // std::cout << "dt = " << dt << " s" << std::endl;

  return 0.0;
}
