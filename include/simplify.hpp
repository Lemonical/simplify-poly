#pragma once

#include <string>

#include "geometry.hpp"

namespace atpps {

// This carries simplification output and target-feasibility metadata
struct SimplificationResult {
    Polygon polygon;
    bool reachedExactTarget = false;
    std::size_t requestedTarget = 0U;
    std::size_t finalVertexCount = 0U;
};

// This performs simplification and is stubbed in commit 1 for pipeline wiring
SimplificationResult SimplifyPolygonToTarget(
    const Polygon& inputPolygon,
    std::size_t targetVertices,
    std::string& noteMessage
);

}  // namespace atpps
