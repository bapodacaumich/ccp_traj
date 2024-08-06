#include "constraints.hpp"
#include "utils.hpp"

#include <casadi/casadi.hpp>
#include <vector>

void enforce_convex_hull_from_points(const std::vector<std::vector<float>>& normals,
                                     const std::vector<std::vector<float>>& points,
                                     casadi::Opti& opti,
                                     std::vector<std::vector<casadi::MX>>& X,
                                     float min_station_distance
                                     ) {

    // enforce convex hull constraints
    size_t n_normals = normals.size();
    size_t n_timesteps = X.size();

    for (size_t j = 0; j < n_timesteps; j++) {
        // get point from states to evaluate constraint
        std::vector<casadi::MX> x = X[j]; // state at timestamp j

        // initialize distance vector
        casadi::MX dists(n_normals,1);

        for (size_t i = 0; i < n_normals; i++) {
            casadi::MX normal = casadi::MX::horzcat({normals[i][0], normals[i][1], normals[i][2]});
            casadi::MX point = casadi::MX::horzcat({points[i][0], points[i][1], points[i][2]});

            dists(i) = (x[0] - points[i][0]) * normals[i][0] + (x[1] - points[i][1]) * normals[i][1] + (x[2] - points[i][2]) * normals[i][2];
        }

        opti.subject_to(mmax(dists) > min_station_distance);
    }
}

void enforce_station_keepout(const std::vector<std::vector<std::vector<float>>>& station_normals,
                             const std::vector<std::vector<std::vector<float>>>& station_points,
                             casadi::Opti& opti,
                             const std::vector<std::vector<float>>& knot_points,
                             casadi::MX& T,
                             float min_station_distance,
                             size_t nT
                             ) {
    // enforce station keepout constraints
    size_t n_obs = station_normals.size();
    size_t n_drift = knot_points.size() - 1;

    // message user about enforcing station keepout constraints
    std::cout << "Enforcing station keepout constraints for " << n_obs << " obstacles for " << n_drift << " drift segments" << std::endl;
    for (size_t i = 0; i < n_drift; i++) {
        // progress bar
        double progress = static_cast<double>(i) / static_cast<double>(n_drift);
        displayProgressBar(progress);

        // get start and end points of drift segment
        std::vector<float> start = knot_points[i];
        std::vector<float> end = knot_points[i+1];

        // get velocity at start of drift segment
        casadi::MX t = T(i,0);
        casadi::MX v_start = cw_v_init(start, end, t);

        // get poses along trajectory
        std::vector<std::vector<casadi::MX>> poses;

        for (size_t j = 0; j < nT; j++) {
            casadi::MX dt = T(i,0) * (j+1) / (nT + 1);
            poses.push_back(cw_pose(start, v_start, dt));
        }

        for (size_t j = 0; j < n_obs; j++) {
            // enforce convex hull constraints
            progress += 1.0 / static_cast<double>(n_obs) / static_cast<double>(n_drift);
            displayProgressBar(progress);
            enforce_convex_hull_from_points(station_normals[j], station_points[j], opti, poses, min_station_distance);
        }
    }
    displayProgressBar(1.0);
    std::cout << std::endl;
}

void enforce_station_keepout(const std::vector<std::vector<std::vector<float>>>& station_normals,
                             const std::vector<std::vector<std::vector<float>>>& station_points,
                             casadi::Opti& opti,
                             const std::vector<std::vector<float>>& knot_points,
                             casadi::MX& X,
                             casadi::MX& T,
                             float min_station_distance,
                             size_t nT
                             ) {
    // enforce station keepout constraints
    size_t n_obs = station_normals.size();
    size_t n_drift = knot_points.size() - 1;

    // message user about enforcing station keepout constraints
    std::cout << "Enforcing station keepout constraints for " << n_obs << " obstacles for " << n_drift * (nT * 2 + 1) << " drift segments" << std::endl;
    for (size_t i = 0; i < n_drift; i++) {
        // progress bar
        double progress = static_cast<double>(i) / static_cast<double>(n_drift);
        displayProgressBar(progress);

        // get start and end points of drift segment
        std::vector<float> start0 = knot_points[i];

        // get velocity at start of drift segment
        casadi::MX t0 = T(2*i,0);
        casadi::MX x0 = X(i,0);
        casadi::MX x1 = X(i,1);
        casadi::MX x2 = X(i,2);
        casadi::MX v_start = cw_v_init(knot_points[i], x0, x1, x2, t0);

        // get poses along trajectory
        std::vector<std::vector<casadi::MX>> poses;
        std::vector<casadi::MX> x = {x0, x1, x2};
        poses.push_back(x);

        for (size_t j = 0; j < nT; j++) {
            casadi::MX dt = t0 * (j+1) / (nT + 1);
            poses.push_back(cw_pose(knot_points[i], v_start, dt));
        }

        // get start and end points of drift segment
        casadi::MX t1 = T(i*2+1,0);
        std::vector<float> end = knot_points[i+1];

        // get velocity at start of drift segment
        v_start = cw_v_init(x0, x1, x2, end, t1);

        for (size_t j = 0; j < nT; j++) {
            casadi::MX dt = t1 * (j+1) / (nT + 1);
            poses.push_back(cw_pose(x0, x1, x2, v_start, dt));
        }

        for (size_t j = 0; j < n_obs; j++) {
            // enforce convex hull constraints
            progress += 1.0 / static_cast<double>(n_obs) / static_cast<double>(n_drift);
            displayProgressBar(progress);
            enforce_convex_hull_from_points(station_normals[j], station_points[j], opti, poses, min_station_distance);
        }
    }
    displayProgressBar(1.0);
    std::cout << std::endl;
}

casadi::MX f(const casadi::MX& x, const casadi::MX& u) {
    float m_I = 5.75f; // inspector mass (kg)
    float mu = 3.986e14f; // standard gravitational parameter
    float a = 6.6e6f; // semi-major axis of international space station orbit
    float n = sqrt(mu / (a * a * a)); // mean motion of ISS -- orbital rate of ISS

    return casadi::MX::vertcat({
        x(3),
        x(4),
        x(5),
        3 * n * n * x(0) + 2 * n * x(4) + u(0) / m_I,
        -2 * n * x(3) + u(1) / m_I,
        -n * n * x(2) + u(2) / m_I
    });
}

void integrate_runge_kutta(const casadi::MX& X, const casadi::MX& U, const std::vector<float>& dt, casadi::Opti& opti) {
    size_t n_timesteps = U.size1();

    casadi::Slice all;
    for (size_t k = 0; k < n_timesteps; k++) {
        casadi::MX k1 = f(X(k, all),                    U(k, all));
        casadi::MX k2 = f(X(k, all) + dt[k] / 2 * k1,   U(k, all));
        casadi::MX k3 = f(X(k, all) + dt[k] / 2 * k2,   U(k, all));
        casadi::MX k4 = f(X(k, all) + dt[k] * k3,       U(k, all));

        casadi::MX x_next = X(all, k) + dt[k] / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
        opti.subject_to(X(k+1, all) == x_next);
    }
}

casadi::MX compute_fuel_cost(const casadi::MX& U, const std::vector<float>& dt) {
    float g0 = 9.81f; // m/s^2
    float Isp = 80.0f; // s

    size_t n_timesteps = U.size1();

    casadi::MX total_impulse = 0.0f;
    casadi::Slice all;
    for (size_t k = 0; k < n_timesteps; k++) {
        total_impulse += sumsqr(U(k, all) * dt[k] * dt[k]);
    }
    casadi::MX fuel_cost = total_impulse / g0 / g0 / Isp / Isp;
    return fuel_cost;
}