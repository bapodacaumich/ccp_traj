#include "constraints.hpp"
#include "ocp_optistack.hpp"
#include "utils.hpp"

#include <casadi/casadi.hpp>
#include <chrono>
#include <iostream>
#include <sstring>
#include <string>
#include <vector


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
        return 0;
    }

    // load file contents into knot_points
    std::vector<std::vector<float>> knot_points;
    loadCSV(matching_file, knot_points);

    size_t n_knots = knot_points.size();

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

    // initialize optistack
    auto opti = casadi::Opti();

    // variables
    size_t 
}