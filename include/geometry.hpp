#pragma once

#include <cstddef>
#include <vector>

namespace atpps {

// This stores a 2D coordinate using double to match input precision
struct Point {
    double x = 0.0;
    double y = 0.0;
};

// This represents a single ring where points are ordered and implicit closure is used
struct Ring {
    int ringId = -1;
    std::vector<Point> vertices;
};

// This represents a polygon-like dataset that can include multiple rings
struct Polygon {
    std::vector<Ring> rings;
};

// This returns the signed area where orientation controls the sign
double ComputeSignedArea(const Ring& ring);

// This accumulates signed area across every ring in the polygon
double ComputeTotalSignedArea(const Polygon& polygon);

// This is a placeholder for assignment areal displacement and will be implemented in a later commit
double ComputeTotalArealDisplacement(const Polygon& inputPolygon, const Polygon& outputPolygon);

// This counts all vertices across all rings for target loop control
std::size_t CountTotalVertices(const Polygon& polygon);

}  // namespace atpps
