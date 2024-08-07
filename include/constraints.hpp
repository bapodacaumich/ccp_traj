#ifndef CONSTRAINTS_HPP
#define CONSTRAINTS_HPP

#include <casadi/casadi.hpp>
#include <vector>

void enforce_convex_hull_from_points(
    const std::vector<std::vector<float>>& normals,
    const std::vector<std::vector<float>>& points,
    casadi::Opti& opti,
    std::vector<std::vector<casadi::MX>>& X,
    float min_station_distance
    );

void enforce_convex_hull_from_points(
    const std::vector<std::vector<float>>& normals,
    const std::vector<std::vector<float>>& points,
    casadi::Opti& opti,
    casadi::MX& X,
    float min_station_distance
    );

void enforce_station_keepout(
    const std::vector<std::vector<std::vector<float>>>& station_normals,
    const std::vector<std::vector<std::vector<float>>>& station_points,
    casadi::Opti& opti,
    const std::vector<std::vector<float>>& knot_points,
    casadi::MX& T,
    float min_station_distance=0.0f,
    size_t nT=3
    );

void enforce_station_keepout(
    const std::vector<std::vector<std::vector<float>>>& station_normals,
    const std::vector<std::vector<std::vector<float>>>& station_points,
    casadi::Opti& opti,
    const std::vector<std::vector<float>>& knot_points,
    casadi::MX& X,
    casadi::MX& T,
    float min_station_distance=0.0f,
    size_t nT=3
    );

void enforce_station_keepout(
    const std::vector<std::vector<std::vector<float>>>& station_normals,
    const std::vector<std::vector<std::vector<float>>>& station_points,
    casadi::Opti& opti,
    casadi::MX& X,
    float min_station_distance=0.0f
    );

// define ode
casadi::MX f(const casadi::MX& x, const casadi::MX& u);

// integrate dynamics
void integrate_runge_kutta(const casadi::MX& X, const casadi::MX& U, const std::vector<float>& dt, casadi::Opti& opti);

// compute fuel cost
casadi::MX compute_fuel_cost(const casadi::MX& U, const std::vector<float>& dt);

// compute knot cost (closest)
casadi::MX compute_knot_cost(const casadi::MX& X, const std::vector<std::vector<float>>& knots, const std::vector<int>& knot_idx);

// compute path cost (squared path length)
casadi::MX compute_path_cost(const casadi::MX& X);

#endif // CONSTRAINTS_HPP