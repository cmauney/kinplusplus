#pragma once
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>


#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/erase.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/algorithm/string_regex.hpp>
#include <boost/lexical_cast.hpp>

#include "defs.h"
#include "environment.h"
#include "species.h"



namespace kinpp {

static auto parse_eqn( const std::string& input ){
  std::vector<double> stochs;
  string_v labels;
  string_v toks;

  boost::split( toks, input, []( char c ) { return c == '+'; } );

  for ( auto& t : toks ) {
    boost::algorithm::erase_all( t, " " );
    auto it = std::find_if(
        t.begin(), t.end(), []( auto c ) { return !::isdigit( c ); } );
    if ( it == t.begin() )
      stochs.push_back(1.0);
    else {
      std::string nstr( t.begin(), it );
      try {
        stochs.push_back( boost::lexical_cast<int>( nstr ) );
      } catch ( const boost::bad_lexical_cast& ) {
        // dunno, leave?
      }
    }
    labels.emplace_back( it, t.end() );
  }

  return std::make_pair(labels, stochs);
}

using rxn_map = std::pair<idx_t, scalar_t>;

template<typename Derived, typename Parameters>
struct reaction_t
{

  std::string       _id;
  std::string       _eqn;

  composition_map_t _spec_depend;

  Parameters        _parameters;

  //composition_map_t     _react_map,   _prod_map;
  std::vector<std::string> _reactlbl, _prodlbl;

  std::vector<rxn_map> _react, _prod, _depend;

  reaction_t() : _id(""), _eqn("") {}
  template<class ...Ps>
  reaction_t( std::string id,
              std::string eqn,
              composition_map_t spec_depend,
              Ps&&... params)
            : _id (std::move(id)),
              _eqn(std::move(eqn)),
              _spec_depend(std::move(spec_depend)),
              _parameters(std::forward<Ps>(params)...)
  {}

  template<class Fn>
  inline auto map_to_idx( const Fn& fnm )
  {
    string_v eqn_toks;
    boost::split(eqn_toks, _eqn, [](char c) { return c == '='; });

    auto do_eqnmap = [&fnm](const auto& eqns, auto& mp, auto& sp){
      auto [tlbl, tstoc] = parse_eqn(eqns);
      for (idx_t i = 0; i < tlbl.size(); ++i)
      {
        mp.emplace_back(fnm(tlbl[i]), tstoc[i]);
        sp.emplace_back( tlbl[i] );
      }
    };

    do_eqnmap(eqn_toks[0], _react, _reactlbl);
    do_eqnmap(eqn_toks[1], _prod, _prodlbl);

    for(const auto& [k, v] : _spec_depend)
    {
      idx_t idx = fnm(k);
      _depend.emplace_back(idx, v);
    }
  }

  inline auto to_string() const {
    std::ostringstream ss;
    auto               pf = [this, &ss]( auto&& x, auto&& y, bool del ) {
      for ( idx_t i = 0; i < x.size(); ++i )
      {
        ss << std::get<scalar_t>( x[i] ) << y[i] << "["<<std::get<idx_t>(x[i])<<"]";
        if ( i == x.size() - 1 ){
          if (!del)
            ss << " = ";
          else
            ss<<std::endl;
        }
        else {
          ss << " + ";
        }
      }
    };

    pf( _react, _reactlbl, false );
    pf( _prod, _prodlbl, true );

    return ss.str();
  }
  constexpr inline auto get_id() const { return _id; }
};


//using rxn_tup = std::tuple<>;
/*#define APPEND_REACTION(rxn)  using rxn##_t = rxn; \
                              using rxn##_v = std::vector<rxn>; \
                              using rxntypes = decltype(std::tuple_cat(std::declval<rxn_tup>(), std::make_tuple(std::declval<rxn##_v>()))); \
*/
}  // namespace kinpp
