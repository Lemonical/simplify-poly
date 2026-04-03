#include "geometry.hpp"

#include <cmath>
#include <limits>
#include <map>

namespace {

// This computes squared point-to-point distance as a base for segment projection distance
double SquaredPointDistance(const atpps::Point& lhs, const atpps::Point& rhs) {
    const double dx = lhs.x - rhs.x;
    const double dy = lhs.y - rhs.y;
    return dx * dx + dy * dy;
}

// This computes squared distance from point p to segment a-b with endpoint clamping
double SquaredDistanceToSegment(const atpps::Point& p, const atpps::Point& a, const atpps::Point& b) {
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    const double segmentLengthSquared = dx * dx + dy * dy;
    if (segmentLengthSquared <= 0.0) {
        return SquaredPointDistance(p, a);
    }

    const double projection = ((p.x - a.x) * dx + (p.y - a.y) * dy) / segmentLengthSquared;
    const double t = (projection < 0.0) ? 0.0 : ((projection > 1.0) ? 1.0 : projection);
    const atpps::Point nearest{a.x + t * dx, a.y + t * dy};
    return SquaredPointDistance(p, nearest);
}

// This sums per-vertex nearest-segment distances from source ring to reference ring boundary
double SumVertexToBoundaryDistances(const atpps::Ring& source, const atpps::Ring& reference) {
    if (source.vertices.empty() || reference.vertices.size() < 2U) {
        return 0.0;
    }

    double total = 0.0;
    const std::size_t referenceCount = reference.vertices.size();
    for (const atpps::Point& vertex : source.vertices) {
        double bestSquared = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0U; i < referenceCount; ++i) {
            const atpps::Point& a = reference.vertices[i];
            const atpps::Point& b = reference.vertices[(i + 1U) % referenceCount];
            const double candidateSquared = SquaredDistanceToSegment(vertex, a, b);
            if (candidateSquared < bestSquared) {
                bestSquared = candidateSquared;
            }
        }
        total += std::sqrt(bestSquared);
    }

    return total;
}

}  // namespace

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

double ComputeTotalRingAreaDrift(const Polygon& inputPolygon, const Polygon& outputPolygon) {
    // This maps ring id to signed area so we can compare ring-wise area drift deterministically
    std::map<int, double> inputAreaByRing;
    for (const Ring& ring : inputPolygon.rings) {
        inputAreaByRing[ring.ringId] = ComputeSignedArea(ring);
    }

    // This maps output ring areas with the same ring id key space
    std::map<int, double> outputAreaByRing;
    for (const Ring& ring : outputPolygon.rings) {
        outputAreaByRing[ring.ringId] = ComputeSignedArea(ring);
    }

    // This accumulates absolute ring-area drift per ring as a stable diagnostic metric
    double totalDrift = 0.0;
    for (const auto& [ringId, inputArea] : inputAreaByRing) {
        const auto outputIt = outputAreaByRing.find(ringId);
        const double outputArea = (outputIt != outputAreaByRing.end()) ? outputIt->second : 0.0;
        totalDrift += std::abs(inputArea - outputArea);
    }

    // This also accounts for rings that only exist in output data
    for (const auto& [ringId, outputArea] : outputAreaByRing) {
        if (inputAreaByRing.find(ringId) == inputAreaByRing.end()) {
            totalDrift += std::abs(outputArea);
        }
    }

    return totalDrift;
}

double ComputeBidirectionalVertexDisplacementProxy(const Polygon& inputPolygon, const Polygon& outputPolygon) {
    // This maps rings by id so displacement proxy stays deterministic under ring ordering
    std::map<int, const Ring*> inputByRingId;
    for (const Ring& ring : inputPolygon.rings) {
        inputByRingId[ring.ringId] = &ring;
    }

    std::map<int, const Ring*> outputByRingId;
    for (const Ring& ring : outputPolygon.rings) {
        outputByRingId[ring.ringId] = &ring;
    }

    // This computes a symmetric proxy by summing output->input and input->output nearest-boundary distances
    double total = 0.0;
    for (const auto& [ringId, inputRing] : inputByRingId) {
        const auto outputIt = outputByRingId.find(ringId);
        if (outputIt == outputByRingId.end()) {
            continue;
        }

        const Ring* outputRing = outputIt->second;
        total += SumVertexToBoundaryDistances(*inputRing, *outputRing);
        total += SumVertexToBoundaryDistances(*outputRing, *inputRing);
    }

    return total;
}

double ComputeTotalArealDisplacement(const Polygon& inputPolygon, const Polygon& outputPolygon) {
    // This keeps compatibility with historical call sites while returning ring-area drift semantics
    return ComputeTotalRingAreaDrift(inputPolygon, outputPolygon);
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
