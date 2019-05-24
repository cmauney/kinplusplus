#pragma once

#include <algorithm>
#include <array>
#include <exception>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "defs.h"
#include "elements.h"

namespace kinpp {
struct NASA9_t {
  std::array<double, 9> coeffs;
  double                Tmin;
  double                Tmax;
};

struct species_t {
  std::string name;

  // std::map<std::string, unsigned int> composition;
  composition_map_t elem_composition;

  std::vector<NASA9_t> nasa9_polynomials;

  species_t( std::string _name, composition_map_t _ec )
      : name( std::move( _name ) ), elem_composition( std::move( _ec ) ) {}
};

struct species_collection_t {
  std::vector<species_t> species;
  species_map_t species_map;

  species_collection_t( const std::string& network_file ) {
    load_species( network_file );
  }

  inline int n_species() { return species.size(); }

  void load_species( const std::string& network_file ) {
    using YAML::Node, YAML::LoadFile;

    Node yf = LoadFile( network_file );

    auto node = yf["species"];

    for ( const auto& s : node ) {
      auto name = s["name"].as<std::string>();

      auto composition = s["composition"].as<composition_map_t>();
      species.emplace_back( name, composition );
    }

    for( idx_t i = 0; i < species.size(); ++i )
    {
      species_map[species[i].name] = i;
    }

  }

  const idx_t get_species_index( const std::string& species_name ) const {
    /*auto sitr = std::find_if(
        std::begin( species ), std::end( species ), [&]( const auto& s ) {
          return s.name == species_name;
        } );

    if ( sitr == species.end() ) {
      std::cout << species_name << std::endl;
      throw new std::runtime_error( "species not found" );
    }
    return ( sitr == species.end()
                 ? sys_constants::idx_notfound
                 : std::distance( std::begin( species ), sitr ) );
                 */
    return species_map.at(species_name);
  }

  const idx_t operator()( const std::string& species_name ) const
  {
    return get_species_index(species_name);
  }

  const std::string get_species_name( const idx_t idx ) const {
    return species[idx].name;
  }

  std::vector<std::string> get_name_vector() {
    std::vector<std::string> name_v;
    for ( const auto& s : species ) name_v.push_back( s.name );
    return name_v;
  }
};

}  // namespace itch
