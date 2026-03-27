#include "geometry.hpp"

#include <cmath>

namespace atpps {

double ComputeSignedArea(const Ring& ring) {
    // This guards against degenerate rings so stats stay predictable
    if (ring.vertices.size() < 3U) {
        return 0.0;
    }

    // This uses the shoelace sum with implicit edge from last back to first
    double twiceArea = 0.0;
    const std::size_t count = ring.vertices.size();
    for (std::size_t i = 0; i < count; ++i) {
        const Point& current = ring.vertices[i];
        const Point& next = ring.vertices[(i + 1U) % count];
        twiceArea += (current.x * next.y) - (next.x * current.y);
    }

    // This converts from twice-area to true area while preserving sign
    return 0.5 * twiceArea;
}

double ComputeTotalSignedArea(const Polygon& polygon) {
    // This folds per-ring signed area into one total metric
    double total = 0.0;
    for (const Ring& ring : polygon.rings) {
        total += ComputeSignedArea(ring);
    }
    return total;
}

double ComputeTotalArealDisplacement(const Polygon& inputPolygon, const Polygon& outputPolygon) {
    // This stays as a stable placeholder until displacement logic is added in a focused commit
    (void)inputPolygon;
    (void)outputPolygon;
    return 0.0;
}

std::size_t CountTotalVertices(const Polygon& polygon) {
    // This sums all ring vertex counts for global target handling
    std::size_t total = 0U;
    for (const Ring& ring : polygon.rings) {
        total += ring.vertices.size();
    }
    return total;
}

}  // namespace atpps
