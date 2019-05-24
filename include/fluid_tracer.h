#ifndef FLUID_TRACER_H
#define FLUID_TRACER_H

#include <boost/math/interpolators/barycentric_rational.hpp>
#include <string>
#include <valarray>
#include <vector>

#include "defs.h"

namespace itch {
typedef boost::math::barycentric_rational<double> interpolate_t;

struct fluid_tracer_t {
  interpolate_t temperature_spline;
  interpolate_t density_spline;

  template<class time_itr, class temperature_itr, class density_itr>
  fluid_tracer_t( time_itr        start_time,
                  time_itr        end_time,
                  temperature_itr temp_start,
                  density_itr     rho_start )
      : temperature_spline( start_time, end_time, temp_start ),
        density_spline( start_time, end_time, rho_start ) {}
  ~fluid_tracer_t() {}
};

}  // namespace itch

#endif  // FLUID_TRACER_H
