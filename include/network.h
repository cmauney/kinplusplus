#pragma once

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "arrhenius_reaction.h"
#include "defs.h"
#include "environment.h"
#include "io_helpers.h"
#include "reaction.h"
#include "species.h"
#include "system_constants.h"

namespace kinpp {

// make sure these corrispond to collection_t
// forward declare here


struct network_t {
using collection_tuple_t = std::tuple<arrhenius_reaction_v>;
//using collection_tuple_t = rxn_tup;
//  using component_idx_pair = std::array<expr_idx_v, 2>;
using component_idx_pair = std::array<idx_v, 2>;
using component_idx_v    = std::vector<component_idx_pair>;

std::string _network_filename;
// species_collection_t species;
collection_tuple_t reactions;

component_idx_v component_idx;

idx_t n_reactions;

network_t( const std::string& network_file ) : _network_filename(network_file), n_reactions( 0 ) {
    load_reactions( network_file );
    LOGI << "network constructed";
}
~network_t() = default;

void load_reactions( const std::string& network_file ) {
  YAML::Node yf = YAML::LoadFile( network_file );

  auto rxn_node = yf["reactions"];

  for ( const auto& rxn : rxn_node ) {
    auto type = rxn["type"].as<std::string>();

    if ( type == "arrhenius" ) {
      std::get<arrhenius_reaction_v>( reactions )
          .emplace_back( rxn.as<arrhenius_reaction_t>() );
    } else {
      LOGI << "\"" << type << "\" not recognized, check input if expected";
      continue;
    }
    LOGI << "\"" << type << "\" reaction imported";
    n_reactions++;
  }
  /*std::apply( [this]( auto&& c ) { c.cement( species ); }, reactions );
  std::apply( [this]( auto&& c ) { n_reactions += c.n_reactions; },
              reactions );*/

  LOGI << "reactions(" << n_reactions << ") loaded";
  }

  inline auto map_species_to_reactions( const species_collection_t& species )
  {
    std::apply( [species](auto&& rv) {
      for(auto& r: rv)
      {
        r.map_to_idx(species);
      }
    }, reactions);
  }

  inline auto print_network() {
    std::ostringstream ss;
    ss << "[NETWORK]" << std::endl;
    ss << sys_constants::headbar << std::endl;
    std::apply(
        [this, &ss]( auto&& rv ) {
          for ( const auto& rxn : rv ) {
            ss << rxn.get_id() << " " << rxn.to_string() << std::endl;
          }
        },
        reactions );
    return ss.str();
  }

};
}  // namespace kinpp
