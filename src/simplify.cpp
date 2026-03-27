#include "simplify.hpp"

#include "geometry.hpp"

namespace atpps {

SimplificationResult SimplifyPolygonToTarget(
    const Polygon& inputPolygon,
    const std::size_t targetVertices,
    std::string& noteMessage
) {
    // This keeps commit 1 deterministic by returning the input unchanged
    SimplificationResult result;
    result.polygon = inputPolygon;
    result.requestedTarget = targetVertices;
    result.finalVertexCount = CountTotalVertices(result.polygon);
    result.reachedExactTarget = (result.finalVertexCount == targetVertices);

    // This note is routed to stderr by the CLI so stdout format stays assignment-safe
    noteMessage = "simplification stub is active and no vertices were removed";
    return result;
}

}  // namespace atpps
