#pragma once

#include <iosfwd>
#include <string>

#include "geometry.hpp"

namespace atpps {

// This loads ring_id,vertex_id,x,y rows and groups them into ring containers
bool LoadPolygonCsv(const std::string& filePath, Polygon& polygon, std::string& errorMessage);

// This writes polygon rows with normalized sequential vertex_id values per ring
void WritePolygonCsv(std::ostream& output, const Polygon& polygon);

}  // namespace atpps
