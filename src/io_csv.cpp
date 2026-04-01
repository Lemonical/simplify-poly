#include "io_csv.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace atpps {

namespace {

struct CsvVertexRow {
    int vertexId = -1;
    Point point;
};

// This trims whitespace and trailing carriage returns from CSV fields
std::string TrimCopy(const std::string& value) {
    std::size_t start = 0U;
    while (start < value.size() && std::isspace(static_cast<unsigned char>(value[start])) != 0) {
        ++start;
    }

    std::size_t end = value.size();
    while (end > start && std::isspace(static_cast<unsigned char>(value[end - 1U])) != 0) {
        --end;
    }

    return value.substr(start, end - start);
}

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

// This validates the known assignment header fields while allowing whitespace around commas
bool IsExpectedHeader(const std::string& headerLine) {
    const std::vector<std::string> fields = SplitCsvLine(headerLine);
    if (fields.size() != 4U) {
        return false;
    }

    return TrimCopy(fields[0]) == "ring_id"
        && TrimCopy(fields[1]) == "vertex_id"
        && TrimCopy(fields[2]) == "x"
        && TrimCopy(fields[3]) == "y";
}

// This prints coordinates with expected-style significant digits while keeping compact decimal text
std::string FormatCoordinate(const double value) {
    std::ostringstream stream;
    stream << std::setprecision(10) << std::defaultfloat << value;
    std::string text = stream.str();

    if (text == "-0" || text == "-0.0") {
        return "0";
    }

    return text;
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

    if (!IsExpectedHeader(header)) {
        errorMessage = "CSV header must be ring_id,vertex_id,x,y";
        return false;
    }

    // This preserves ring order by ring_id using std::map for deterministic iteration
    std::map<int, std::vector<CsvVertexRow>> ringVertices;
    std::map<int, std::unordered_set<int>> seenVertexIds;

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
            const int ringId = std::stoi(TrimCopy(fields[0]));
            const int vertexId = std::stoi(TrimCopy(fields[1]));
            const double x = std::stod(TrimCopy(fields[2]));
            const double y = std::stod(TrimCopy(fields[3]));

            if (ringId < 0) {
                errorMessage = "ring_id must be non-negative";
                return false;
            }

            if (vertexId < 0) {
                errorMessage = "vertex_id must be non-negative";
                return false;
            }

            auto& ringSeenVertexIds = seenVertexIds[ringId];
            const bool inserted = ringSeenVertexIds.insert(vertexId).second;
            if (!inserted) {
                errorMessage = "duplicate vertex_id found in ring";
                return false;
            }

            ringVertices[ringId].push_back(CsvVertexRow{vertexId, Point{x, y}});
        } catch (const std::invalid_argument&) {
            errorMessage = "failed to parse numeric CSV value";
            return false;
        } catch (const std::out_of_range&) {
            errorMessage = "numeric CSV value out of range";
            return false;
        }
    }

    for (auto& [ringId, vertices] : ringVertices) {
        std::sort(vertices.begin(), vertices.end(), [](const CsvVertexRow& lhs, const CsvVertexRow& rhs) {
            return lhs.vertexId < rhs.vertexId;
        });

        Ring ring;
        ring.ringId = ringId;
        ring.vertices.reserve(vertices.size());
        for (const CsvVertexRow& row : vertices) {
            ring.vertices.push_back(row.point);
        }

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
            output << ring.ringId << ',' << vertexId << ','
                   << FormatCoordinate(point.x) << ','
                   << FormatCoordinate(point.y) << '\n';
        }
    }
}

}  // namespace atpps
