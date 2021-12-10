#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "parameters.hpp"
#include "simulation.hpp"

#include <array>
#include <string>

std::array<size_t, 5> do_analysis(Param all_parameters);
std::string get_outcome(const std::array<size_t, 5>& cell_counts);

#endif // ANALYSIS_H
