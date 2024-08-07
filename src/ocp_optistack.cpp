#include "constraints.hpp"
#include "ocp_optistack.hpp"
#include "time_interval_struct.hpp"
#include "utils.hpp"

#include <casadi/casadi.hpp>
#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


std::vector<casadi::DM> ocp_optistack(std::string vgd,
                                      bool locality,
                                      float k_weight,
                                      float p_weight,
                                      float f_weight,
                                      float thrust_limit,
                                      float min_station_distance,
                                      float knot_cost_normalization,
                                      float fuel_cost_normalization
                                      ) {

    std::string locality_str = "";
    if (locality) {locality_str = "_local";}

    // initialize time recording variables
    auto start_time = std::chrono::high_resolution_clock::now();

    // find file that matches vgd, locality, and knot directory
    std::string knotdir = "../data/knot_points/";
    std::string matching_file = firstFileWithPrefix(knotdir, vgd+locality_str);
    if (matching_file == "") {
        std::cerr << "No matching file found in directory: " << knotdir << std::endl;
        return {};
    }

    // load file contents into knot_points
    std::vector<std::vector<float>> knot_points;
    loadCSV(matching_file, knot_points);

    // load in station normals and centroids
    std::vector<std::vector<std::vector<float>>> station_normals;
    std::vector<std::vector<std::vector<float>>> station_centroids;

    size_t num_obstacles = 15;
    for (size_t i = 0; i < num_obstacles; i++){
        std::string normal_file = "../data/model/normal_centroid/obs" + std::to_string(i+1) + "_normals.txt";
        std::string centroid_file = "../data/model/normal_centroid/obs" + std::to_string(i+1) + "_points.txt";

        std::vector<std::vector<float>> normals;
        std::vector<std::vector<float>> centroids;

        loadCSV(normal_file, normals);
        loadCSV(centroid_file, centroids);

        station_normals.push_back(normals);
        station_centroids.push_back(centroids);
    }

    // display time to load knot points and station normals and centroids
    auto load_files_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = load_files_time - start_time;
    std::cout << "Time to load knot points and station: " << elapsed.count() << " seconds " << std::endl;

    // hyperparameters
    float velocity = 0.2f;
    size_t n_timesteps = 400;

    //compute time intervals for each optimized state
    TimeInterval time_intervals = compute_time_intervals(knot_points, velocity, n_timesteps);
    n_timesteps = time_intervals.dt.size();
    casadi::DM initial_path = linear_initial_path(knot_points, time_intervals);

    // initialize optistack
    auto opti = casadi::Opti();

    // optimization variables
    casadi::MX X = opti.variable(n_timesteps+1, 6);
    casadi::MX U = opti.variable(n_timesteps, 3);

    // CONSTRAINTS
    // integrate runge kutta using clohessy wiltshire ode
    integrate_runge_kutta(X, U, time_intervals.dt, opti);

    // constrain start pose
    float min_start_distance = 0.5f; // meters
    casadi::MX start_error = (X(0,0) - knot_points[0][0]) * (X(0,0) - knot_points[0][0]) + (X(0,1) - knot_points[0][1]) * (X(0,1) - knot_points[0][1]) + (X(0,2) - knot_points[0][2]) * (X(0,2) - knot_points[0][2]);
    opti.subject_to(start_error <= min_start_distance * min_start_distance);

    // impose thrust limit constraints
    for (size_t i = 0; i < n_timesteps; i++) {
        opti.subject_to( U(i,0) * U(i,0) + U(i,1) * U(i,1) + U(i,2) * U(i,2) <= thrust_limit * thrust_limit);
    }

    // compute fuel cost
    casadi::MX fuel_cost = compute_fuel_cost(U, time_intervals.dt);

    // compute knot cost
    casadi::MX knot_cost = compute_knot_cost(X, knot_points, time_intervals.knot_idx);

    // compute path cost (squared path length)
    casadi::MX path_cost = compute_path_cost(X);

    // build objective function
    casadi::MX cost = f_weight * fuel_cost / fuel_cost_normalization + k_weight * knot_cost / knot_cost_normalization + p_weight * path_cost;

    // direct optistack to minimize objective cost above
    opti.minimize(cost);

    // convex hull obstacle
    std::cout << "Enforce station keepout" << std::endl;
    enforce_station_keepout(station_normals, station_centroids, opti, X, min_station_distance);

    // set initial X
    opti.set_initial(X, initial_path);

    // set initial U
    opti.set_initial(U, casadi::DM::zeros(n_timesteps, 3));

    // -- solver --
    casadi::Dict p_opts;
    casadi::Dict s_opts;
    s_opts["print_level"] = 8;
    s_opts["tol"] = 1e-3;
    s_opts["max_iter"] = 10000;

    std::cout << "Solving OCP..." << std::endl;
    opti.solver("ipopt", p_opts, s_opts);

    casadi::OptiSol sol = opti.solve();

    // display solution
    std::cout << "Solution States =" << sol.value(X) << std::endl;
    std::cout << "Solution Thrust Input =" << sol.value(U) << std::endl;
    std::cout << "Final Solution Cost =" << sol.value(cost) << std::endl;

    return {sol.value(X), sol.value(U)};
}