#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>

#include "geometry.hpp"
#include "io_csv.hpp"
#include "simplify.hpp"
#include "validation.hpp"

namespace {

// This caps expensive global validation on extremely large single-ring datasets
constexpr std::size_t kValidationVertexLimit = 50000U;

// This env toggle allows opting back into old timeout-focused behavior for very large single-ring startup validation
bool ShouldSkipLargeSingleRingInputValidation() {
    const char* raw = std::getenv("ATPPS_SKIP_LARGE_SINGLE_RING_INPUT_VALIDATION");
    if (raw == nullptr || raw[0] == '\0') {
        return false;
    }

    const std::string value(raw);
    return value == "1" || value == "true" || value == "TRUE" || value == "on" || value == "ON";
}

// This centralizes usage text so error and help paths stay consistent
void PrintUsage(const char* executableName) {
    std::cerr << "Usage: " << executableName << " <input_file> <target_vertices>\n";
}

}  // namespace

int main(int argc, char** argv) {
    // This keeps argument handling strict so automation has predictable behavior
    if (argc != 3) {
        PrintUsage(argv[0]);
        return EXIT_FAILURE;
    }

    const std::string inputFilePath = argv[1];

    std::size_t targetVertices = 0U;
    try {
        const int parsedTarget = std::stoi(argv[2]);
        if (parsedTarget <= 0) {
            std::cerr << "target_vertices must be a positive integer\n";
            return EXIT_FAILURE;
        }
        targetVertices = static_cast<std::size_t>(parsedTarget);
    } catch (const std::exception&) {
        std::cerr << "failed to parse target_vertices\n";
        return EXIT_FAILURE;
    }

    atpps::Polygon inputPolygon;
    std::string ioError;
    if (!atpps::LoadPolygonCsv(inputFilePath, inputPolygon, ioError)) {
        std::cerr << "input error: " << ioError << '\n';
        return EXIT_FAILURE;
    }

    const std::size_t inputVertexCount = atpps::CountTotalVertices(inputPolygon);
    const bool skipLargeSingleRingInputValidation = ShouldSkipLargeSingleRingInputValidation();
    const bool skipGlobalValidation =
        skipLargeSingleRingInputValidation
        && inputPolygon.rings.size() == 1U
        && inputVertexCount >= kValidationVertexLimit;

    if (!skipGlobalValidation) {
        std::string topologyError;
        if (!atpps::ValidatePolygonTopology(inputPolygon, topologyError)) {
            std::cerr << "topology error: " << topologyError << '\n';
            return EXIT_FAILURE;
        }
    }

    std::string simplifyNote;
    const atpps::SimplificationResult result =
        atpps::SimplifyPolygonToTarget(inputPolygon, targetVertices, simplifyNote);

    // This writes the simplified polygon CSV block first as required
    atpps::WritePolygonCsv(std::cout, result.polygon);

    // This computes assignment metrics even in stub mode so output shape is stable
    const double inputArea = atpps::ComputeTotalSignedArea(inputPolygon);
    const double outputArea = atpps::ComputeTotalSignedArea(result.polygon);
    const double arealDisplacement = result.totalArealDisplacement;

    // This uses scientific output format to match expected test file style
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Total signed area in input: " << inputArea << '\n';
    std::cout << "Total signed area in output: " << outputArea << '\n';
    std::cout << "Total areal displacement: " << arealDisplacement << '\n';

    // This intentionally suppresses optional diagnostic notes so output comparisons stay clean

    return EXIT_SUCCESS;
}
