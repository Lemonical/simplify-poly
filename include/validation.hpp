#pragma once

#include <string>

#include "geometry.hpp"

namespace atpps {

// This validates ring cardinality and edge intersection topology constraints
bool ValidatePolygonTopology(const Polygon& polygon, std::string& errorMessage);

}  // namespace atpps
