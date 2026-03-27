#include "io_csv.hpp"

#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace atpps {

namespace {

// This splits one CSV line without handling quoted commas because test data is numeric only
std::vector<std::string> SplitCsvLine(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) {
        fields.push_back(field);
    }
    return fields;
}

}  // namespace

bool LoadPolygonCsv(const std::string& filePath, Polygon& polygon, std::string& errorMessage) {
    // This resets output so callers never observe partial state on failures
    polygon.rings.clear();

    std::ifstream file(filePath);
    if (!file.is_open()) {
        errorMessage = "failed to open input file";
        return false;
    }

    std::string header;
    if (!std::getline(file, header)) {
        errorMessage = "input file is empty";
        return false;
    }

    // This preserves ring order by ring_id using std::map for deterministic iteration
    std::map<int, std::vector<Point>> ringVertices;

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }

        const std::vector<std::string> fields = SplitCsvLine(line);
        if (fields.size() != 4U) {
            errorMessage = "invalid CSV row field count";
            return false;
        }

        try {
            const int ringId = std::stoi(fields[0]);
            const double x = std::stod(fields[2]);
            const double y = std::stod(fields[3]);
            ringVertices[ringId].push_back(Point{x, y});
        } catch (...) {
            errorMessage = "failed to parse numeric CSV value";
            return false;
        }
    }

    for (const auto& [ringId, vertices] : ringVertices) {
        Ring ring;
        ring.ringId = ringId;
        ring.vertices = vertices;
        polygon.rings.push_back(ring);
    }

    return true;
}

void WritePolygonCsv(std::ostream& output, const Polygon& polygon) {
    // This prints the assignment-required CSV header first
    output << "ring_id,vertex_id,x,y\n";

    // This emits rows in deterministic ring order and contiguous vertex ids per ring
    for (const Ring& ring : polygon.rings) {
        for (std::size_t vertexId = 0U; vertexId < ring.vertices.size(); ++vertexId) {
            const Point& point = ring.vertices[vertexId];
            output << ring.ringId << ',' << vertexId << ',' << point.x << ',' << point.y << '\n';
        }
    }
}

}  // namespace atpps
