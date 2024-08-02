#ifndef TIME_INTERVAL_STRUCT_HPP
#define TIME_INTERVAL_STRUCT_HPP

#include <vector>

struct TimeInterval {
    std::vector<float> dt;
    std::vector<int> knot_idx;
};

#endif // TIME_INTERVAL_STRUCT_HPP