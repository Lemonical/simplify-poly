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
    double totalArealDisplacement = 0.0;
};

// This performs deterministic topology-safe simplification toward a target vertex count
SimplificationResult SimplifyPolygonToTarget(
    const Polygon& inputPolygon,
    std::size_t targetVertices,
    std::string& noteMessage
);

}  // namespace atpps
