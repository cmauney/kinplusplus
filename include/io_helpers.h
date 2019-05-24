#pragma once

#include <fstream>
#include <tuple>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/erase.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>

#include <yaml-cpp/yaml.h>

#include "arrhenius_reaction.h"
#include "configurator.h"
#include "defs.h"
#include "network.h"

namespace YAML {

template<>
struct convert<kinpp::arrhenius_reaction_t> {
  static bool decode( const Node& node, kinpp::arrhenius_reaction_t& rhs ) {
    auto eqn    = node["equation"].as<std::string>();
    auto rid    = node["id"].as<std::string>();
    auto params = node["parameters"];
    // auto Mspc   = node["M_species"];
    // auto Meff   = node["M_efficiencies"];
    auto a = params["k"][0].as<double>();
    auto b = params["k"][1].as<double>();
    auto c = params["k"][2].as<double>();

    auto M = kinpp::composition_map_t();
    if ( node["M_eff"] ) auto M = node["M_eff"].as<kinpp::composition_map_t>();

    // auto parms = std::make_tuple(a, b, c );

    // gotta hack it fn
    // rhs.~arrhenius_reaction_t();
    // new ( &rhs ) itch::arrhenius_reaction_t( rid, eqn, a, b, c );
    rhs = kinpp::arrhenius_reaction_t(
        rid, eqn, M, a, b , c );

    return true;
  }
};

template<>
struct convert<kinpp::composition_map_t> {
  static bool decode( const Node& node, kinpp::composition_map_t& rhs ) {
    // if ( !node.IsSequence() ) { return false; }

    // rhs = std::make_tuple( node[0].as<std::string>(), node[1].as<double>() );
    for ( auto it = node.begin(); it != node.end(); ++it )
      rhs[it->first.as<std::string>()] = it->second.as<double>();

    return true;
  }

  static Node encode(kinpp::composition_map_t& rhs)
  {
    Node node;
    for(const auto& [k,v] : rhs)
    {
      node[k] = v;
    }
    return node;
  }
};

}  // namespace YAML
