#pragma once

#include <fstream>
#include <map>
#include <string>
#include <variant>
#include <vector>

#include <plog/Log.h>
#include <boost/algorithm/string/join.hpp>
#include <boost/program_options.hpp>

#include "defs.h"

namespace options = boost::program_options;

namespace kinpp {

struct configurator_t {
  std::string network_filename;
  std::string abundance_filename;
  std::string element_filename;
  std::string thermochem_filename;

  double ode_step_tolerance;
  double ode_rho_inf;
  idx_t  ode_corrector_steps;

  double initial_timestep;
  double minimum_timestep;
  double maximum_timestep;

  double time_start;
  double time_finish;
  idx_t  n_maximum_timesteps;

  string_v follow_species;
  idx_t    n_store_every;
  idx_t    n_print_every;
  idx_t    n_dump_every;

  options::options_description desc;
  options::variables_map       vm;

  configurator_t() : desc( "configuration" ) {
    desc.add_options()( "network_filename",
                        options::value<std::string>( &network_filename ),
                        "YAML file with chemistry data" );

    desc.add_options()( "abundance_filename",
                        options::value<std::string>( &abundance_filename ),
                        "YAML file with chemunit abundance" );

    desc.add_options()( "element_filename",
                        options::value<std::string>( &element_filename ),
                        "YAML file with elemental data" );

    desc.add_options()( "thermochem_filename",
                        options::value<std::string>( &thermochem_filename ),
                        "BURCAT thermochemistry data" );

    desc.add_options()(
        "ode_step_tolerance",
        options::value<double>( &ode_step_tolerance )->default_value( 1.0E-6 ),
        "stepper tolerance" );

    desc.add_options()(
        "ode_rho_inf",
        options::value<double>( &ode_rho_inf )->default_value( 0.5 ),
        "rho_inf for generalized alpha stepper" );

    desc.add_options()(
        "ode_corrector_steps",
        options::value<idx_t>( &ode_corrector_steps )->default_value( 4 ),
        "number of predictor-corrector steps" );

    desc.add_options()(
        "minimum_timestep",
        options::value<double>( &minimum_timestep )->default_value( 1.0E-10 ),
        "minimum allowed timestep (in seconds)" );

    desc.add_options()(
        "maximum_timestep",
        options::value<double>( &maximum_timestep )->default_value( 1.0E4 ),
        "maximum allowed timestep (in seconds)" );

    desc.add_options()(
        "initial_timestep",
        options::value<double>( &initial_timestep )->default_value( 1.0E-2 ),
        "initial timestep (in seconds)" );

    desc.add_options()(
        "time_start",
        options::value<double>( &time_start )->default_value( 0.0 ),
        "initial time value (in seconds)" );

    desc.add_options()(
        "time_finish",
        options::value<double>( &time_finish )->default_value( 1.0 ),
        "final time value (in seconds)" );

    desc.add_options()( "n_maximum_timesteps",
                        options::value<idx_t>( &n_maximum_timesteps )
                            ->default_value( 1000000000ul ),
                        "maximum number of timesteps" );

    desc.add_options()(
        "follow_species",
        options::value<string_v>( &follow_species )->multitoken(),
        "list of species to follow" );

    desc.add_options()(
        "n_store_every",
        options::value<idx_t>( &n_store_every )->default_value( 10 ),
        "store every N step in memory" );

    desc.add_options()(
        "n_print_every",
        options::value<idx_t>( &n_print_every )->default_value( 1000 ),
        "print every N steps" );

    desc.add_options()(
        "n_dump_every",
        options::value<idx_t>( &n_dump_every )->default_value( 1000000 ),
        "dump solution every N steps" );
  }

  void load( const std::string& cfg_file ) {
    std::ifstream config_file( cfg_file.c_str() );

    vm = options::variables_map();

    options::store( options::parse_config_file( config_file, desc ), vm );
    options::notify( vm );

    LOGI << "configuration loaded";

    LOGI << "input\n";
    LOGI << "::network_file = \"" << network_filename << "\"";
    LOGI << "::abundance_file = \"" << thermochem_filename << "\"";
    LOGI << "::element_file = \"" << element_filename << "\"";
    LOGI << "::thermochem_file = \"" << thermochem_filename << "\"";

    LOGI << "ode\n";
    LOGI << "::step_tolerance = " << ode_step_tolerance;
    LOGI << "::corrector_steps = " << ode_corrector_steps;
    LOGI << "::minimum_timestep = " << minimum_timestep;
    LOGI << "::maximum_timestep = " << maximum_timestep;
    LOGI << "::initial_timestep = " << initial_timestep;

    LOGI << "physics";

    LOGI << "runtime";
    LOGI << "::time_start = " << time_start;
    LOGI << "::time_finish = " << time_finish;
    LOGI << "::n_maximum_timesteps = " << n_maximum_timesteps;

    LOGI << "obs";
    LOGI << "::n_store_every = " << n_store_every;
    LOGI << "::n_print_every = " << n_print_every;
    LOGI << "::n_dump_every = " << n_dump_every;
  }
};

}  // namespace kinpp
