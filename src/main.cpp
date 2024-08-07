#include "ocp_optistack.hpp"
#include "utils.hpp"

#include <chrono>

int main(int argc, char* argv[]) {
    // argc: number of arguments
    // argv: array of arguments

    std::string vgd = "2.0m";
    bool locality = true;

    if (argc > 1) {
        vgd = std::string(argv[1]);

        if (argc > 2) {
            locality = (argv[2][0] == 't');
        }
    }

    std::string locality_str = "";
    if (locality) {locality_str = "_local";}

    auto start_time = std::chrono::high_resolution_clock::now();
    std::vector<casadi::DM> sol = ocp_optistack(vgd, locality);

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "Time to set up and solve OCP: " << elapsed.count() << " seconds " << std::endl;

    // save out solution to csv files
    saveCSV("../data/solution/" + vgd + locality_str + "_x.csv", sol[0]);
    saveCSV("../data/solution/" + vgd + locality_str + "_u.csv", sol[1]);
}