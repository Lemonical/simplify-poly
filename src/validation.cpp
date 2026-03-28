#include "validation.hpp"

#include <cmath>

namespace {

// This epsilon keeps orientation and equality checks stable under floating-point noise
constexpr double kEpsilon = 1e-9;

// This compares scalars with a symmetric tolerance window
bool NearlyEqual(const double lhs, const double rhs) {
    return std::abs(lhs - rhs) <= kEpsilon;
}

// This compares coordinates so geometric predicates can treat near-identical points as equal
bool PointsEqual(const atpps::Point& lhs, const atpps::Point& rhs) {
    return NearlyEqual(lhs.x, rhs.x) && NearlyEqual(lhs.y, rhs.y);
}

// This computes the oriented area factor for triangle a-b-c
double OrientationValue(const atpps::Point& a, const atpps::Point& b, const atpps::Point& c) {
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

// This classifies orientation sign with tolerance to avoid unstable branch flips
int OrientationSign(const atpps::Point& a, const atpps::Point& b, const atpps::Point& c) {
    const double value = OrientationValue(a, b, c);
    if (value > kEpsilon) {
        return 1;
    }
    if (value < -kEpsilon) {
        return -1;
    }
    return 0;
}

// This checks whether p lies on segment a-b including endpoints under tolerance
bool IsPointOnSegment(const atpps::Point& p, const atpps::Point& a, const atpps::Point& b) {
    if (OrientationSign(a, b, p) != 0) {
        return false;
    }

    const double minX = (a.x < b.x) ? a.x : b.x;
    const double maxX = (a.x > b.x) ? a.x : b.x;
    const double minY = (a.y < b.y) ? a.y : b.y;
    const double maxY = (a.y > b.y) ? a.y : b.y;

    const bool inX = p.x >= (minX - kEpsilon) && p.x <= (maxX + kEpsilon);
    const bool inY = p.y >= (minY - kEpsilon) && p.y <= (maxY + kEpsilon);
    return inX && inY;
}

// This returns true for proper intersections and collinear overlap/touch cases
bool SegmentsIntersect(
    const atpps::Point& a1,
    const atpps::Point& a2,
    const atpps::Point& b1,
    const atpps::Point& b2
) {
    const int o1 = OrientationSign(a1, a2, b1);
    const int o2 = OrientationSign(a1, a2, b2);
    const int o3 = OrientationSign(b1, b2, a1);
    const int o4 = OrientationSign(b1, b2, a2);

    if (o1 != o2 && o3 != o4) {
        return true;
    }

    if (o1 == 0 && IsPointOnSegment(b1, a1, a2)) {
        return true;
    }
    if (o2 == 0 && IsPointOnSegment(b2, a1, a2)) {
        return true;
    }
    if (o3 == 0 && IsPointOnSegment(a1, b1, b2)) {
        return true;
    }
    if (o4 == 0 && IsPointOnSegment(a2, b1, b2)) {
        return true;
    }

    return false;
}

// This identifies adjacency where shared endpoints are expected in a simple closed ring
bool AreAdjacentEdgeIndices(const std::size_t edgeA, const std::size_t edgeB, const std::size_t edgeCount) {
    if (edgeA == edgeB) {
        return true;
    }

    if ((edgeA + 1U) % edgeCount == edgeB) {
        return true;
    }

    if ((edgeB + 1U) % edgeCount == edgeA) {
        return true;
    }

    return false;
}

// This validates local ring hygiene before expensive intersection passes
bool ValidateRingPrimitiveShape(const atpps::Ring& ring, std::string& errorMessage) {
    // This enforces minimum ring size so closed-edge logic remains valid
    if (ring.vertices.size() < 3U) {
        errorMessage = "ring has fewer than 3 vertices";
        return false;
    }

    const std::size_t count = ring.vertices.size();

    // This rejects zero-length edges, including closure from last to first
    for (std::size_t i = 0U; i < count; ++i) {
        const atpps::Point& current = ring.vertices[i];
        const atpps::Point& next = ring.vertices[(i + 1U) % count];
        if (PointsEqual(current, next)) {
            errorMessage = "ring contains a zero-length edge";
            return false;
        }
    }

    // This rejects repeated points that are not consecutive because they imply a pinched ring
    for (std::size_t i = 0U; i < count; ++i) {
        for (std::size_t j = i + 1U; j < count; ++j) {
            if (!PointsEqual(ring.vertices[i], ring.vertices[j])) {
                continue;
            }

            const bool neighbors = (j == i + 1U) || (i == 0U && j + 1U == count);
            if (!neighbors) {
                errorMessage = "ring contains duplicate non-adjacent vertices";
                return false;
            }
        }
    }

    return true;
}

// This checks one ring against itself while ignoring expected adjacent-edge endpoint touches
bool ValidateRingSelfIntersection(const atpps::Ring& ring, std::string& errorMessage) {
    const std::size_t edgeCount = ring.vertices.size();

    for (std::size_t edgeA = 0U; edgeA < edgeCount; ++edgeA) {
        const atpps::Point& a1 = ring.vertices[edgeA];
        const atpps::Point& a2 = ring.vertices[(edgeA + 1U) % edgeCount];

        for (std::size_t edgeB = edgeA + 1U; edgeB < edgeCount; ++edgeB) {
            if (AreAdjacentEdgeIndices(edgeA, edgeB, edgeCount)) {
                continue;
            }

            const atpps::Point& b1 = ring.vertices[edgeB];
            const atpps::Point& b2 = ring.vertices[(edgeB + 1U) % edgeCount];

            if (SegmentsIntersect(a1, a2, b1, b2)) {
                errorMessage = "ring self-intersection detected";
                return false;
            }
        }
    }

    return true;
}

// This checks edge intersections between two different rings
bool ValidateInterRingIntersection(const atpps::Ring& lhs, const atpps::Ring& rhs, std::string& errorMessage) {
    const std::size_t lhsCount = lhs.vertices.size();
    const std::size_t rhsCount = rhs.vertices.size();

    for (std::size_t lhsEdge = 0U; lhsEdge < lhsCount; ++lhsEdge) {
        const atpps::Point& a1 = lhs.vertices[lhsEdge];
        const atpps::Point& a2 = lhs.vertices[(lhsEdge + 1U) % lhsCount];

        for (std::size_t rhsEdge = 0U; rhsEdge < rhsCount; ++rhsEdge) {
            const atpps::Point& b1 = rhs.vertices[rhsEdge];
            const atpps::Point& b2 = rhs.vertices[(rhsEdge + 1U) % rhsCount];

            if (SegmentsIntersect(a1, a2, b1, b2)) {
                errorMessage = "intersection detected between different rings";
                return false;
            }
        }
    }

    return true;
}

}  // namespace

namespace atpps {

bool ValidatePolygonTopology(const Polygon& polygon, std::string& errorMessage) {
    // This validates each ring primitive shape and self-intersection constraints first
    for (const Ring& ring : polygon.rings) {
        if (!ValidateRingPrimitiveShape(ring, errorMessage)) {
            errorMessage = "ring " + std::to_string(ring.ringId) + ": " + errorMessage;
            return false;
        }

        if (!ValidateRingSelfIntersection(ring, errorMessage)) {
            errorMessage = "ring " + std::to_string(ring.ringId) + ": " + errorMessage;
            return false;
        }
    }

    // This then checks pairwise ring edge interactions to preserve multi-ring topology
    for (std::size_t i = 0U; i < polygon.rings.size(); ++i) {
        for (std::size_t j = i + 1U; j < polygon.rings.size(); ++j) {
            if (!ValidateInterRingIntersection(polygon.rings[i], polygon.rings[j], errorMessage)) {
                errorMessage = "ring " + std::to_string(polygon.rings[i].ringId)
                    + " vs ring " + std::to_string(polygon.rings[j].ringId)
                    + ": " + errorMessage;
                return false;
            }
        }
    }

    // This clears any previous message once every topology rule passes
    errorMessage.clear();
    return true;
}

}  // namespace atpps
