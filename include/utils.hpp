#ifndef UTILS_HPP
#define UTILS_HPP

#include <casadi/casadi.hpp>
#include <string>
#include <vector>
#include "time_interval_struct.hpp"

// Clohessy Wiltshire computation for casadi variables -- overload to handle casadi::MX for start and end
casadi::MX cw_v_init(const std::vector<float>& start, casadi::MX& dx, casadi::MX& dy, casadi::MX& dz, casadi::MX& t);
casadi::MX cw_v_init(casadi::MX& x, casadi::MX& y, casadi::MX& z, const std::vector<float>& end, casadi::MX& t);
casadi::MX cw_v_init(const std::vector<float> &start, const std::vector<float> &end, casadi::MX& t);
casadi::DM cw_v_init(const std::vector<float> &start, const std::vector<float> &end, casadi::DM& t);

casadi::MX cw_v_end(const std::vector<float> &start, const std::vector<float> &end, casadi::MX& t);

std::vector<casadi::MX> cw_pose(const std::vector<float>& start, casadi::MX& v_start, casadi::MX& t);
std::vector<casadi::DM> cw_pose(const std::vector<float>& start, casadi::DM& v_start, casadi::DM& t);
std::vector<casadi::MX> cw_pose(casadi::MX& x, casadi::MX& y, casadi::MX& z, casadi::MX& v_start, casadi::MX& t);

casadi::MX compute_path_cost(casadi::MX& T, const std::vector<std::vector<float>> &knot_points, bool square = true);

// loading csv files for knotpoints
bool loadCSV(const std::string& filename, std::vector<std::vector<float>>& data, int rowlen = 3);

// save out solutions to csv files
void saveCSV(const std::string& filename, const casadi::DM& data);
void saveCSV(const std::string& filename, const std::vector<std::vector<float>>& data);

// find first file with prefix in directory
std::string firstFileWithPrefix(const std::string& directory, const std::string& prefix);

// print 2d vector (debugging)
void print2DVector(const std::vector<std::vector<float>>& vec);

// display progress bar
void displayProgressBar(double progress, int width = 150);

// get initial intermediate points for warm start
casadi::DM get_intermediate_points_init(const std::vector<std::vector<float>>& knot_points, const casadi::DM& T_init);
casadi::DM linear_initial_path(const std::vector<std::vector<float>>& knot_points, const TimeInterval& time_intervals);

// get time intervals for each segment and state point
TimeInterval compute_time_intervals(const std::vector<std::vector<float>>& knot_points, float velocity, size_t n_timesteps);


#endif // UTILS_HPP