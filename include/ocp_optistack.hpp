#ifndef OCP_OPTISTACK_HPP
#define OCP_OPTISTACK_HPP

#include <casadi/casadi.hpp>
#include <string>

std::vector<casadi::DM> ocp_optistack(std::string vgd,
                                      bool locality,
                                      float k_weight=1.0f,
                                      float p_weight=1.0f,
                                      float f_weight=1.0f,
                                      float thrust_limit=0.2f,
                                      float min_station_distance=1.0f,
                                      float knot_cost_normalization=100.0f,
                                      float fuel_cost_normalization=0.00001f
                                      );

#endif // OCP_OPTISTACK_HPP