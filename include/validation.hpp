#pragma once

#include <string>

#include "geometry.hpp"

namespace atpps {

// This validates baseline polygon structure and will expand to full topology checks later
bool ValidatePolygonTopology(const Polygon& polygon, std::string& errorMessage);

}  // namespace atpps
