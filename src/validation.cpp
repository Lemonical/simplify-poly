#include "validation.hpp"

namespace atpps {

bool ValidatePolygonTopology(const Polygon& polygon, std::string& errorMessage) {
    // This only enforces minimum ring cardinality in commit 1 to keep bootstrap simple
    for (const Ring& ring : polygon.rings) {
        if (ring.vertices.size() < 3U) {
            errorMessage = "ring has fewer than 3 vertices";
            return false;
        }
    }

    // This returns true for now and full intersection checks will be added later
    errorMessage.clear();
    return true;
}

}  // namespace atpps
