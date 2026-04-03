#include "simplify.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <sstream>
#include <vector>

#include "geometry.hpp"
#include "validation.hpp"

namespace {

// This epsilon keeps tie-break comparisons stable when doubles are nearly equal
constexpr double kEpsilon = 1e-12;

// This tolerance controls ring-area restoration stop criteria in relative terms
constexpr double kAreaRestoreRelativeTolerance = 1e-9;

// This tolerance guards very small absolute area differences near zero-area limits
constexpr double kAreaRestoreAbsoluteTolerance = 1e-9;

// This cap keeps the restoration pass deterministic and bounded
constexpr std::size_t kMaxAreaRestoreIterations = 12U;

// This captures one removable-vertex option with deterministic ordering fields
struct RemovalCandidate {
    std::size_t ringIndex = 0U;
    std::size_t vertexIndex = 0U;
    int ringId = -1;
    atpps::Point point;
    double shapeError = std::numeric_limits<double>::infinity();
    double areaDelta = std::numeric_limits<double>::infinity();
};

// This stores one pair-collapse option where two neighbors are replaced by one computed point
struct PairCollapseCandidate {
    std::size_t ringIndex = 0U;
    std::size_t indexB = 0U;
    std::size_t indexC = 0U;
    int ringId = -1;
    int containmentDepth = 0;
    atpps::Point replacement;
    double cost = std::numeric_limits<double>::infinity();
    double displacement = 0.0;
};

// This carries one lazy heap entry for the single-ring fast simplification mode
struct HeapCandidate {
    std::size_t vertexIndex = 0U;
    std::size_t generation = 0U;
    atpps::Point point;
    double shapeError = std::numeric_limits<double>::infinity();
    double areaDelta = std::numeric_limits<double>::infinity();
};

// This carries one lazy heap entry for single-ring pair-collapse simplification
struct SingleRingPairCandidate {
    std::size_t indexB = 0U;
    std::size_t generationB = 0U;
    std::size_t generationC = 0U;
    std::size_t recentAge = std::numeric_limits<std::size_t>::max();
    std::size_t bucketIndex = 0U;
    atpps::Point replacement;
    double displacement = std::numeric_limits<double>::infinity();
    double cornerSharpness = std::numeric_limits<double>::infinity();
    double cost = std::numeric_limits<double>::infinity();
};

// This stores optional trace output state so normal runs remain unchanged
struct TraceState {
    bool initialized = false;
    bool enabled = false;
    std::ofstream file;
};

// This lazily initializes tracing from ATPPS_TRACE_FILE and keeps one shared trace sink
TraceState& GetTraceState() {
    static TraceState state;
    if (state.initialized) {
        return state;
    }

    state.initialized = true;
    const char* path = std::getenv("ATPPS_TRACE_FILE");
    if (path == nullptr || path[0] == '\0') {
        return state;
    }

    state.file.open(path, std::ios::out | std::ios::trunc);
    state.enabled = state.file.is_open();
    return state;
}

// This reports whether trace logging is currently active
bool TraceEnabled() {
    return GetTraceState().enabled;
}

// This writes one trace line only when tracing is enabled
void TraceWriteLine(const std::string& line) {
    TraceState& state = GetTraceState();
    if (!state.enabled) {
        return;
    }

    state.file << line << '\n';
}

// This formats one point in compact text for decision-trace output
std::string FormatPointForTrace(const atpps::Point& point) {
    std::ostringstream output;
    output << point.x << "," << point.y;
    return output.str();
}

// This formats one pair-collapse candidate so traces are easy to scan
std::string FormatPairCandidateForTrace(const PairCollapseCandidate& candidate) {
    std::ostringstream output;
    output << "ring=" << candidate.ringId
           << " depth=" << candidate.containmentDepth
           << " b=" << candidate.indexB
           << " c=" << candidate.indexC
           << " E=(" << FormatPointForTrace(candidate.replacement) << ')'
           << " cost=" << candidate.cost
           << " disp=" << candidate.displacement;
    return output.str();
}

// This formats one ring vertex list so final trace snapshots can be compared line-by-line
std::string FormatRingForTrace(const atpps::Ring& ring) {
    std::ostringstream output;
    output << "ring=" << ring.ringId << " n=" << ring.vertices.size();
    for (std::size_t i = 0U; i < ring.vertices.size(); ++i) {
        output << " | " << i << ':' << ring.vertices[i].x << ',' << ring.vertices[i].y;
    }
    return output.str();
}

// This computes squared distance from p to segment a-b for local shape-error scoring
double SquaredDistanceToSegment(const atpps::Point& p, const atpps::Point& a, const atpps::Point& b) {
    const double dx = b.x - a.x;
    const double dy = b.y - a.y;
    const double lengthSquared = (dx * dx) + (dy * dy);

    if (lengthSquared <= kEpsilon) {
        const double px = p.x - a.x;
        const double py = p.y - a.y;
        return (px * px) + (py * py);
    }

    const double projection = ((p.x - a.x) * dx + (p.y - a.y) * dy) / lengthSquared;
    const double t = (projection < 0.0) ? 0.0 : ((projection > 1.0) ? 1.0 : projection);
    const double nearestX = a.x + (t * dx);
    const double nearestY = a.y + (t * dy);
    const double diffX = p.x - nearestX;
    const double diffY = p.y - nearestY;
    return (diffX * diffX) + (diffY * diffY);
}

// This returns the signed area of triangle a-b-c used for local area-change scoring
double SignedTriangleArea(const atpps::Point& a, const atpps::Point& b, const atpps::Point& c) {
    const double cross = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    return 0.5 * cross;
}

// This returns squared distance between two points for scale and movement normalization
double SquaredDistanceBetweenPoints(const atpps::Point& lhs, const atpps::Point& rhs) {
    const double dx = lhs.x - rhs.x;
    const double dy = lhs.y - rhs.y;
    return dx * dx + dy * dy;
}

// This estimates local squared-length scale from the four-edge neighborhood around one pair collapse
double ComputeSingleRingLocalLengthScaleSq(
    const atpps::Point& A,
    const atpps::Point& B,
    const atpps::Point& C,
    const atpps::Point& D
) {
    const double ab = SquaredDistanceBetweenPoints(A, B);
    const double bc = SquaredDistanceBetweenPoints(B, C);
    const double cd = SquaredDistanceBetweenPoints(C, D);
    const double da = SquaredDistanceBetweenPoints(D, A);

    // This uses average edge scale so movement normalization is stable across varied ring sizes
    const double averageEdgeLengthSq = 0.25 * (ab + bc + cd + da);
    return std::max(averageEdgeLengthSq, kEpsilon);
}

// This estimates local area scale for displacement normalization with a length-based fallback on near-collinear neighborhoods
double ComputeSingleRingLocalAreaScale(
    const atpps::Point& A,
    const atpps::Point& B,
    const atpps::Point& C,
    const atpps::Point& D,
    const double localLengthScaleSq
) {
    const double areaABC = std::abs(SignedTriangleArea(A, B, C));
    const double areaACD = std::abs(SignedTriangleArea(A, C, D));

    // This keeps displacement normalization well-defined even when the local configuration is almost collinear
    const double localAreaScale = areaABC + areaACD;
    return std::max(localAreaScale, 0.5 * localLengthScaleSq);
}

// This computes one ring-level squared-length scale from average edge lengths for cross-candidate normalization
double ComputeRingLengthScaleSq(const atpps::Ring& ring) {
    if (ring.vertices.empty()) {
        return 1.0;
    }

    double sumSquaredEdgeLengths = 0.0;
    for (std::size_t i = 0U; i < ring.vertices.size(); ++i) {
        const atpps::Point& current = ring.vertices[i];
        const atpps::Point& next = ring.vertices[(i + 1U) % ring.vertices.size()];
        sumSquaredEdgeLengths += SquaredDistanceBetweenPoints(current, next);
    }

    const double averageSquaredEdgeLength = sumSquaredEdgeLengths / static_cast<double>(ring.vertices.size());
    return std::max(averageSquaredEdgeLength, kEpsilon);
}

// This computes one ring-level area scale with deterministic fallback for near-collinear rings
double ComputeRingAreaScale(const atpps::Ring& ring, const double ringLengthScaleSq) {
    const double ringAreaMagnitude = std::abs(atpps::ComputeSignedArea(ring));
    return std::max(ringAreaMagnitude, 0.5 * ringLengthScaleSq);
}

// This computes local corner sharpness as abs(sin(theta)) so 0 means collinear and 1 means orthogonal
double ComputeCornerSharpness(const atpps::Point& prev, const atpps::Point& pivot, const atpps::Point& next) {
    const double v1x = prev.x - pivot.x;
    const double v1y = prev.y - pivot.y;
    const double v2x = next.x - pivot.x;
    const double v2y = next.y - pivot.y;

    const double len1 = std::hypot(v1x, v1y);
    const double len2 = std::hypot(v2x, v2y);
    if (len1 <= 1e-12 || len2 <= 1e-12) {
        return 0.0;
    }

    const double crossMagnitude = std::abs(v1x * v2y - v1y * v2x);
    return crossMagnitude / (len1 * len2);
}

// This maps simplification progress to [0,1] so early stages can preserve corners more aggressively
double ComputeSingleRingProgress(
    const std::size_t initialCount,
    const std::size_t aliveCount,
    const std::size_t targetVertices
) {
    if (initialCount <= targetVertices) {
        return 1.0;
    }

    const double removed = static_cast<double>(initialCount - aliveCount);
    const double required = static_cast<double>(initialCount - targetVertices);
    if (required <= 0.0) {
        return 1.0;
    }

    const double ratio = removed / required;
    if (ratio < 0.0) {
        return 0.0;
    }
    if (ratio > 1.0) {
        return 1.0;
    }
    return ratio;
}

// This maps remaining simplification work to [0,1] so schedule bands can target early, mid, and late phases
double ComputeSingleRingRemainingRatio(
    const std::size_t initialCount,
    const std::size_t aliveCount,
    const std::size_t targetVertices
) {
    if (initialCount <= targetVertices) {
        return 0.0;
    }

    const double remaining = static_cast<double>(aliveCount - targetVertices);
    const double required = static_cast<double>(initialCount - targetVertices);
    if (required <= 0.0) {
        return 0.0;
    }

    const double ratio = remaining / required;
    if (ratio < 0.0) {
        return 0.0;
    }
    if (ratio > 1.0) {
        return 1.0;
    }
    return ratio;
}

// This picks one fixed schedule-band value by remaining-ratio to keep policy deterministic across runs
double SelectSingleRingProgressBandValue(
    const double remainingRatio,
    const double earlyValue,
    const double midValue,
    const double lateValue
) {
    if (remainingRatio > (2.0 / 3.0)) {
        return earlyValue;
    }
    if (remainingRatio > (1.0 / 3.0)) {
        return midValue;
    }
    return lateValue;
}

// This builds deterministic local metrics for one candidate removal
RemovalCandidate BuildCandidate(
    const atpps::Ring& ring,
    const std::size_t ringIndex,
    const std::size_t vertexIndex
) {
    const std::size_t count = ring.vertices.size();
    const std::size_t prevIndex = (vertexIndex + count - 1U) % count;
    const std::size_t nextIndex = (vertexIndex + 1U) % count;

    const atpps::Point& prev = ring.vertices[prevIndex];
    const atpps::Point& current = ring.vertices[vertexIndex];
    const atpps::Point& next = ring.vertices[nextIndex];

    RemovalCandidate candidate;
    candidate.ringIndex = ringIndex;
    candidate.vertexIndex = vertexIndex;
    candidate.ringId = ring.ringId;
    candidate.point = current;

    // This favors removals that keep the local curve closest to its neighbor chord
    candidate.shapeError = SquaredDistanceToSegment(current, prev, next);

    // This favors removals that perturb local signed area as little as possible
    candidate.areaDelta = std::abs(SignedTriangleArea(prev, current, next));
    return candidate;
}

// This compares doubles with epsilon while keeping strict weak ordering behavior
bool LessWithTolerance(const double lhs, const double rhs) {
    if (lhs + kEpsilon < rhs) {
        return true;
    }
    return false;
}

// This gives deterministic strict ordering for heap candidates while preserving epsilon-aware ties
bool IsBetterHeapCandidate(const HeapCandidate& lhs, const HeapCandidate& rhs) {
    if (LessWithTolerance(lhs.shapeError, rhs.shapeError)) {
        return true;
    }
    if (LessWithTolerance(rhs.shapeError, lhs.shapeError)) {
        return false;
    }

    if (LessWithTolerance(lhs.areaDelta, rhs.areaDelta)) {
        return true;
    }
    if (LessWithTolerance(rhs.areaDelta, lhs.areaDelta)) {
        return false;
    }

    if (lhs.vertexIndex != rhs.vertexIndex) {
        return lhs.vertexIndex < rhs.vertexIndex;
    }

    if (LessWithTolerance(lhs.point.x, rhs.point.x)) {
        return true;
    }
    if (LessWithTolerance(rhs.point.x, lhs.point.x)) {
        return false;
    }

    if (LessWithTolerance(lhs.point.y, rhs.point.y)) {
        return true;
    }
    if (LessWithTolerance(rhs.point.y, lhs.point.y)) {
        return false;
    }

    return false;
}

// This comparator makes std::priority_queue behave as a min-heap by candidate quality
struct HeapCandidateWorse {
    bool operator()(const HeapCandidate& lhs, const HeapCandidate& rhs) const {
        return IsBetterHeapCandidate(rhs, lhs);
    }
};

// This enforces deterministic ordering for single-ring pair-collapse candidates
bool IsBetterSingleRingPairCandidate(const SingleRingPairCandidate& lhs, const SingleRingPairCandidate& rhs) {
    if (LessWithTolerance(lhs.cost, rhs.cost)) {
        return true;
    }
    if (LessWithTolerance(rhs.cost, lhs.cost)) {
        return false;
    }

    return lhs.indexB < rhs.indexB;
}

// This comparator gives std::priority_queue min-heap behavior for pair-collapse cost
struct SingleRingPairCandidateWorse {
    bool operator()(const SingleRingPairCandidate& lhs, const SingleRingPairCandidate& rhs) const {
        return IsBetterSingleRingPairCandidate(rhs, lhs);
    }
};

/****************************
*      CORE THRESHOLDS      *
*****************************/

constexpr std::size_t kHugeSingleRingVertexThreshold = 50000U;  // This threshold switches extremely large single-ring inputs to linear-time deterministic sampling
constexpr std::size_t kMultiRingBeamWidth = 6U;                 // This bounds multi-ring sequence search width for deterministic runtime
constexpr std::size_t kMultiRingBranchFactor = 6U;              // This bounds per-state candidate expansion during multi-ring sequence search
constexpr double kRingRoleComparableCostMultiplier = 1.0;       // This multiplier defines when two pair-collapse costs are close enough to let ring-role priority decide
constexpr double kMultiRingMovementWeight = 1e-12;              // This scales movement influence in multi-ring pair-collapse scoring

/****************************
*    SINGLE RING SCORING    *
*****************************/

constexpr double kSingleRingMovementWeight = 0.1;                 // This scales movement influence in single-ring heap scoring to calibrate sequence behavior
constexpr double kSingleRingMovementStartFactor = 1.0;            // This scales single-ring movement influence at simplification start so early-stage behavior can be calibrated
constexpr double kSingleRingMovementEndFactor = 1.0;              // This scales single-ring movement influence near target so late-stage behavior can be calibrated
constexpr double kSingleRingFeaturePenaltyMultiplier = 0.75;      // This scales early-stage corner-preservation pressure in single-ring collapse scoring
constexpr int kSingleRingStagedPolicyMode = 0;                    // This controls whether staged single-ring corner-preservation policy is enabled
constexpr int kSingleRingScoreMode = 0;                           // This controls which single-ring base score formula is used before optional staged corner penalties
constexpr int kSingleRingDispersionMode = 0;                      // This controls whether deterministic single-ring dispersion policy is disabled, hard-filtered, or soft-penalized
constexpr std::size_t kSingleRingDispersionWindow = 8U;           // This defines how many collapse steps count as locally recent for dispersion logic
constexpr double kSingleRingDispersionPenaltyMultiplier = 1.5;    // This scales soft dispersion penalties when nearby collapses happened very recently

/****************************
*   SINGLE RING SELECTION   *
*****************************/

// Mode 4 keeps cost-first ordering and only applies bucket-age preference when costs are comparable
// Mode 5 adds corner-smoothness preference inside comparable-cost cohorts
// Mode 6 adaptively switches between mode 3 and mode 4 based on one-step window cost spread
// Mode 7 adaptively switches between mode 1 and mode 3 based on initial ring vertex count
constexpr int kSingleRingSelectionMode = 7;                        // This controls whether single-ring collapse always takes heap-best candidate or uses windowed dispersion-aware tie selection
constexpr std::size_t kSingleRingAdaptiveSizeThreshold = 2000U;    // This defines the initial-vertex threshold where adaptive mode 7 switches from mode 1 (small rings) to mode 3 (larger rings)
constexpr std::size_t kSingleRingSelectionWindow = 8U;             // This bounds how many valid heap candidates are considered in one windowed single-ring selection step
constexpr int kSingleRingWindowPoolMode = 1;                       // This controls whether single-ring window candidates are composed from an expanded pooled frontier before prefiltering
constexpr std::size_t kSingleRingWindowPoolExpandFactor = 2U;      // This expands pooled frontier size before deterministic downselection into the final single-ring window
constexpr std::size_t kSingleRingWindowPoolMinCount = 8U;          // This avoids pool-composition logic when too few pooled candidates are available

/****************************
*  PROGRESS SCHEDULE BANDS  *
****************************/

constexpr int kSingleRingProgressScheduleMode = 0;              // This controls single-ring progress schedules: 0 off, 1 pool-only, 2 pool plus prefilter thresholds
constexpr double kSingleRingWindowPoolExpandEarly = 3.0;        // This sets single-ring pool expansion in the early phase when most removals remain
constexpr double kSingleRingWindowPoolExpandMid = 2.0;          // This sets single-ring pool expansion in the mid phase when roughly half removals remain
constexpr double kSingleRingWindowPoolExpandLate = 1.5;         // This sets single-ring pool expansion in the late phase near the target vertex count

/****************************
*    WINDOW PREFILTERING    *
*****************************/

constexpr int kSingleRingWindowPrefilterMode = 1;                      // This controls whether deterministic single-ring window prefiltering is disabled or displacement-percentile based
constexpr double kSingleRingWindowPrefilterPercentile = 0.25;          // This keeps only low-displacement candidates up to this percentile inside one single-ring selection window
constexpr double kSingleRingWindowPrefilterEarly = 0.35;               // This sets displacement prefilter percentile in the early phase when broad shape pruning is still acceptable
constexpr double kSingleRingWindowPrefilterMid = 0.25;                 // This sets displacement prefilter percentile in the mid phase
constexpr double kSingleRingWindowPrefilterLate = 0.15;                // This sets displacement prefilter percentile in the late phase when picks should be tighter
constexpr std::size_t kSingleRingWindowPrefilterMinCount = 4U;         // This avoids applying window prefilter logic on very small candidate windows
constexpr int kSingleRingWindowCostPrefilterMode = 1;                  // This controls whether deterministic single-ring window cost-prefiltering is disabled or enabled
constexpr double kSingleRingWindowCostPrefilterPercentile = 0.5;       // This keeps only low-cost candidates up to this percentile after displacement prefiltering
constexpr double kSingleRingWindowCostPrefilterEarly = 0.70;           // This sets cost prefilter percentile in the early phase to keep more alternatives alive
constexpr double kSingleRingWindowCostPrefilterMid = 0.50;             // This sets cost prefilter percentile in the mid phase
constexpr double kSingleRingWindowCostPrefilterLate = 0.35;            // This sets cost prefilter percentile in the late phase for tighter final edits
constexpr std::size_t kSingleRingWindowCostPrefilterMinCount = 3U;     // This avoids applying window cost-prefilter logic on very small already-prefiltered windows

/****************************
*      AGE TIE HANDLING     *
*****************************/

constexpr int kSingleRingWindowAgePrefilterMode = 0;                   // This controls whether deterministic single-ring recent-age prefiltering is disabled, always enabled, or enabled only for comparable-cost pools
constexpr double kSingleRingWindowAgePrefilterPercentile = 0.5;        // This keeps only a top fraction of high recent-age candidates from the post-cost pool
constexpr std::size_t kSingleRingWindowAgePrefilterMinCount = 3U;      // This avoids applying window recent-age prefilter logic on very small pools
constexpr std::size_t kSingleRingBucketCount = 32U;                    // This sets how many deterministic arc buckets are used when bucketed single-ring selection mode is active
constexpr std::size_t kSingleRingTraceSteps = 200U;                    // This bounds how many single-ring collapse steps are written to trace so long runs remain inspectable
constexpr double kSingleRingSelectionTieCostMultiplier = 0.1;          // This defines when single-ring candidate costs are comparable enough for recent-age tie preference
constexpr double kSingleRingAdaptiveModeCostSpreadThreshold = 0.08;    // This defines when adaptive single-ring selection mode treats one window as tight-cost and switches to mode-4 behavior

/****************************
*      TOPOLOGY GUARDS      *
*****************************/

constexpr int kSingleRingTopologyGuardMode = 1;                     // This controls whether accepted single-ring pair collapses are filtered by a local self-intersection guard
constexpr std::size_t kSingleRingTopologyGuardVertexLimit = 20000U; // This caps the ring size where the single-ring topology guard runs so very large cases stay tractable
constexpr int kMultiRingDepthQuotaMode = 0;                         // This controls whether multi-ring candidate truncation enforces depth-parity quotas
constexpr int kMultiRingNormalizedScoreMode = 0;                    // This controls whether multi-ring pair-collapse scoring uses local scale-normalized terms

/****************************
*      ENV READ HELPERS     *
*****************************/

// This reads one unsigned integer env value once and clamps it to a safe deterministic range
std::size_t ReadEnvSizeTClamped(
    const char* envName,
    const std::size_t fallback,
    const std::size_t minimum,
    const std::size_t maximum
) {
    const char* raw = std::getenv(envName);
    if (raw == nullptr || raw[0] == '\0') {
        return fallback;
    }

    char* parseEnd = nullptr;
    const unsigned long long parsed = std::strtoull(raw, &parseEnd, 10);
    if (parseEnd == raw || (parseEnd != nullptr && parseEnd[0] != '\0')) {
        return fallback;
    }

    std::size_t value = static_cast<std::size_t>(parsed);
    if (value < minimum) {
        value = minimum;
    }
    if (value > maximum) {
        value = maximum;
    }
    return value;
}

// This reads one boolean env flag with common truthy spellings and a deterministic fallback
bool ReadEnvBool(const char* envName, const bool fallback) {
    const char* raw = std::getenv(envName);
    if (raw == nullptr || raw[0] == '\0') {
        return fallback;
    }

    const std::string value(raw);
    if (value == "1" || value == "true" || value == "TRUE" || value == "on" || value == "ON") {
        return true;
    }
    if (value == "0" || value == "false" || value == "FALSE" || value == "off" || value == "OFF") {
        return false;
    }

    return fallback;
}

// This reads one integer env value once and clamps it to a safe deterministic range
int ReadEnvIntClamped(const char* envName, const int fallback, const int minimum, const int maximum) {
    const char* raw = std::getenv(envName);
    if (raw == nullptr || raw[0] == '\0') {
        return fallback;
    }

    char* parseEnd = nullptr;
    const long parsed = std::strtol(raw, &parseEnd, 10);
    if (parseEnd == raw || (parseEnd != nullptr && parseEnd[0] != '\0')) {
        return fallback;
    }

    int value = static_cast<int>(parsed);
    if (value < minimum) {
        value = minimum;
    }
    if (value > maximum) {
        value = maximum;
    }
    return value;
}

// This reads one floating-point env value once and clamps it to a safe deterministic range
double ReadEnvDoubleClamped(
    const char* envName,
    const double fallback,
    const double minimum,
    const double maximum
) {
    const char* raw = std::getenv(envName);
    if (raw == nullptr || raw[0] == '\0') {
        return fallback;
    }

    char* parseEnd = nullptr;
    const double parsed = std::strtod(raw, &parseEnd);
    if (parseEnd == raw || (parseEnd != nullptr && parseEnd[0] != '\0')) {
        return fallback;
    }

    double value = parsed;
    if (value < minimum) {
        value = minimum;
    }
    if (value > maximum) {
        value = maximum;
    }
    return value;
}

// This allows runtime beam-width calibration while preserving default deterministic behavior
std::size_t GetConfiguredMultiRingBeamWidth() {
    static const std::size_t configured =
        ReadEnvSizeTClamped("ATPPS_MULTI_RING_BEAM_WIDTH", kMultiRingBeamWidth, 1U, 64U);
    return configured;
}

// This allows runtime branch-factor calibration while preserving default deterministic behavior
std::size_t GetConfiguredMultiRingBranchFactor() {
    static const std::size_t configured =
        ReadEnvSizeTClamped("ATPPS_MULTI_RING_BRANCH_FACTOR", kMultiRingBranchFactor, 1U, 64U);
    return configured;
}

// This allows runtime tuning of huge single-ring fallback threshold while preserving default behavior
std::size_t GetConfiguredHugeSingleRingVertexThreshold() {
    static const std::size_t configured =
        ReadEnvSizeTClamped("ATPPS_SINGLE_RING_HUGE_THRESHOLD", kHugeSingleRingVertexThreshold, 1000U, 5000000U);
    return configured;
}

// This allows runtime tuning of ring-role comparability without changing default behavior
double GetConfiguredRingRoleComparableCostMultiplier() {
    static const double configured =
        ReadEnvDoubleClamped("ATPPS_RING_ROLE_COST_MULTIPLIER", kRingRoleComparableCostMultiplier, 0.0, 10.0);
    return configured;
}

// This allows runtime tuning of multi-ring movement influence without changing default behavior
double GetConfiguredMultiRingMovementWeight() {
    static const double configured =
        ReadEnvDoubleClamped("ATPPS_MULTI_RING_MOVEMENT_WEIGHT", kMultiRingMovementWeight, 0.0, 10.0);
    return configured;
}

// This allows runtime tuning of single-ring corner-preservation pressure without changing default behavior
double GetConfiguredSingleRingFeaturePenaltyMultiplier() {
    static const double configured =
        ReadEnvDoubleClamped(
            "ATPPS_SINGLE_RING_FEATURE_PENALTY_MULT",
            kSingleRingFeaturePenaltyMultiplier,
            0.0,
            10.0
        );
    return configured;
}

// This allows runtime tuning of normalized single-ring movement influence without changing default behavior
double GetConfiguredSingleRingMovementWeight() {
    static const double configured =
        ReadEnvDoubleClamped("ATPPS_SINGLE_RING_MOVEMENT_WEIGHT", kSingleRingMovementWeight, 0.0, 100.0);
    return configured;
}

// This allows runtime tuning of movement-weight schedule at simplification start
double GetConfiguredSingleRingMovementStartFactor() {
    static const double configured =
        ReadEnvDoubleClamped(
            "ATPPS_SINGLE_RING_MOVEMENT_START_FACTOR",
            kSingleRingMovementStartFactor,
            0.0,
            10.0
        );
    return configured;
}

// This allows runtime tuning of movement-weight schedule near target
double GetConfiguredSingleRingMovementEndFactor() {
    static const double configured =
        ReadEnvDoubleClamped(
            "ATPPS_SINGLE_RING_MOVEMENT_END_FACTOR",
            kSingleRingMovementEndFactor,
            0.0,
            10.0
        );
    return configured;
}

// This allows runtime toggling of staged single-ring policy without changing default behavior
int GetConfiguredSingleRingStagedPolicyMode() {
    static const int configured =
        ReadEnvIntClamped("ATPPS_SINGLE_RING_STAGED_POLICY_MODE", kSingleRingStagedPolicyMode, 0, 1);
    return configured;
}

// This allows runtime switching between deterministic single-ring base scoring formulas
int GetConfiguredSingleRingScoreMode() {
    static const int configured =
        ReadEnvIntClamped("ATPPS_SINGLE_RING_SCORE_MODE", kSingleRingScoreMode, 0, 2);
    return configured;
}

// This allows runtime switching of deterministic single-ring dispersion policy behavior
int GetConfiguredSingleRingDispersionMode() {
    static const int configured =
        ReadEnvIntClamped("ATPPS_SINGLE_RING_DISPERSION_MODE", kSingleRingDispersionMode, 0, 2);
    return configured;
}

// This allows runtime tuning of the recent-collapse window used by single-ring dispersion policy
std::size_t GetConfiguredSingleRingDispersionWindow() {
    static const std::size_t configured =
        ReadEnvSizeTClamped("ATPPS_SINGLE_RING_DISPERSION_WINDOW", kSingleRingDispersionWindow, 1U, 2048U);
    return configured;
}

// This allows runtime tuning of soft single-ring dispersion penalty strength
double GetConfiguredSingleRingDispersionPenaltyMultiplier() {
    static const double configured =
        ReadEnvDoubleClamped(
            "ATPPS_SINGLE_RING_DISPERSION_PENALTY_MULT",
            kSingleRingDispersionPenaltyMultiplier,
            0.0,
            10.0
        );
    return configured;
}

// This allows runtime switching between strict heap-best and windowed single-ring candidate selection
int GetConfiguredSingleRingSelectionMode() {
    static const int configured =
        ReadEnvIntClamped("ATPPS_SINGLE_RING_SELECTION_MODE", kSingleRingSelectionMode, 0, 7);
    return configured;
}

// This allows runtime tuning of adaptive single-ring size threshold used by mode 7
std::size_t GetConfiguredSingleRingAdaptiveSizeThreshold() {
    static const std::size_t configured =
        ReadEnvSizeTClamped(
            "ATPPS_SINGLE_RING_ADAPTIVE_SIZE_THRESHOLD",
            kSingleRingAdaptiveSizeThreshold,
            100U,
            1000000U
        );
    return configured;
}

// This allows runtime tuning of window size used by deterministic single-ring candidate selection
std::size_t GetConfiguredSingleRingSelectionWindow() {
    static const std::size_t configured =
        ReadEnvSizeTClamped("ATPPS_SINGLE_RING_SELECTION_WINDOW", kSingleRingSelectionWindow, 1U, 64U);
    return configured;
}

// This allows runtime switching of pre-window pool composition before single-ring prefilter stages
int GetConfiguredSingleRingWindowPoolMode() {
    static const int configured =
        ReadEnvIntClamped("ATPPS_SINGLE_RING_WINDOW_POOL_MODE", kSingleRingWindowPoolMode, 0, 1);
    return configured;
}

// This allows runtime tuning of pooled frontier expansion before single-ring window downselection
std::size_t GetConfiguredSingleRingWindowPoolExpandFactor() {
    static const std::size_t configured =
        ReadEnvSizeTClamped(
            "ATPPS_SINGLE_RING_WINDOW_POOL_EXPAND_FACTOR",
            kSingleRingWindowPoolExpandFactor,
            1U,
            16U
        );
    return configured;
}

// This allows runtime tuning of minimum pooled size where window-composition logic is applied
std::size_t GetConfiguredSingleRingWindowPoolMinCount() {
    static const std::size_t configured =
        ReadEnvSizeTClamped(
            "ATPPS_SINGLE_RING_WINDOW_POOL_MIN_COUNT",
            kSingleRingWindowPoolMinCount,
            2U,
            256U
        );
    return configured;
}

// This allows runtime toggling of progress-band schedules for single-ring pool and optional prefilter parameters
int GetConfiguredSingleRingProgressScheduleMode() {
    static const int configured =
        ReadEnvIntClamped("ATPPS_SINGLE_RING_PROGRESS_SCHEDULE_MODE", kSingleRingProgressScheduleMode, 0, 2);
    return configured;
}

// This allows runtime switching of deterministic single-ring window prefilter behavior
int GetConfiguredSingleRingWindowPrefilterMode() {
    static const int configured =
        ReadEnvIntClamped("ATPPS_SINGLE_RING_WINDOW_PREFILTER_MODE", kSingleRingWindowPrefilterMode, 0, 1);
    return configured;
}

// This allows runtime tuning of displacement percentile used by single-ring window prefiltering
double GetConfiguredSingleRingWindowPrefilterPercentile() {
    static const double configured =
        ReadEnvDoubleClamped(
            "ATPPS_SINGLE_RING_WINDOW_PREFILTER_PERCENTILE",
            kSingleRingWindowPrefilterPercentile,
            0.0,
            1.0
        );
    return configured;
}

// This allows runtime tuning of minimum window size where single-ring prefiltering is applied
std::size_t GetConfiguredSingleRingWindowPrefilterMinCount() {
    static const std::size_t configured =
        ReadEnvSizeTClamped(
            "ATPPS_SINGLE_RING_WINDOW_PREFILTER_MIN_COUNT",
            kSingleRingWindowPrefilterMinCount,
            2U,
            64U
        );
    return configured;
}

// This allows runtime switching of deterministic single-ring cost prefilter behavior
int GetConfiguredSingleRingWindowCostPrefilterMode() {
    static const int configured =
        ReadEnvIntClamped(
            "ATPPS_SINGLE_RING_WINDOW_COST_PREFILTER_MODE",
            kSingleRingWindowCostPrefilterMode,
            0,
            1
        );
    return configured;
}

// This allows runtime tuning of cost percentile used by single-ring cost prefiltering
double GetConfiguredSingleRingWindowCostPrefilterPercentile() {
    static const double configured =
        ReadEnvDoubleClamped(
            "ATPPS_SINGLE_RING_WINDOW_COST_PREFILTER_PERCENTILE",
            kSingleRingWindowCostPrefilterPercentile,
            0.0,
            1.0
        );
    return configured;
}

// This allows runtime tuning of minimum post-displacement pool size where cost prefiltering is applied
std::size_t GetConfiguredSingleRingWindowCostPrefilterMinCount() {
    static const std::size_t configured =
        ReadEnvSizeTClamped(
            "ATPPS_SINGLE_RING_WINDOW_COST_PREFILTER_MIN_COUNT",
            kSingleRingWindowCostPrefilterMinCount,
            2U,
            64U
        );
    return configured;
}

// This allows runtime switching of deterministic single-ring recent-age prefilter behavior
int GetConfiguredSingleRingWindowAgePrefilterMode() {
    static const int configured =
        ReadEnvIntClamped(
            "ATPPS_SINGLE_RING_WINDOW_AGE_PREFILTER_MODE",
            kSingleRingWindowAgePrefilterMode,
            0,
            2
        );
    return configured;
}

// This allows runtime tuning of retained high-age fraction used by single-ring age prefiltering
double GetConfiguredSingleRingWindowAgePrefilterPercentile() {
    static const double configured =
        ReadEnvDoubleClamped(
            "ATPPS_SINGLE_RING_WINDOW_AGE_PREFILTER_PERCENTILE",
            kSingleRingWindowAgePrefilterPercentile,
            0.0,
            1.0
        );
    return configured;
}

// This allows runtime tuning of minimum post-cost pool size where age prefiltering is applied
std::size_t GetConfiguredSingleRingWindowAgePrefilterMinCount() {
    static const std::size_t configured =
        ReadEnvSizeTClamped(
            "ATPPS_SINGLE_RING_WINDOW_AGE_PREFILTER_MIN_COUNT",
            kSingleRingWindowAgePrefilterMinCount,
            2U,
            64U
        );
    return configured;
}

// This allows runtime tuning of cost-comparability for single-ring windowed tie selection
double GetConfiguredSingleRingSelectionTieCostMultiplier() {
    static const double configured =
        ReadEnvDoubleClamped(
            "ATPPS_SINGLE_RING_SELECTION_TIE_MULT",
            kSingleRingSelectionTieCostMultiplier,
            0.0,
            5.0
        );
    return configured;
}

// This allows runtime tuning of adaptive single-ring selector spread threshold that switches between mode 3 and mode 4
double GetConfiguredSingleRingAdaptiveModeCostSpreadThreshold() {
    static const double configured =
        ReadEnvDoubleClamped(
            "ATPPS_SINGLE_RING_ADAPTIVE_SPREAD_THRESHOLD",
            kSingleRingAdaptiveModeCostSpreadThreshold,
            0.0,
            10.0
        );
    return configured;
}

// This allows runtime tuning of arc-bucket count used by bucketed single-ring selection mode
std::size_t GetConfiguredSingleRingBucketCount() {
    static const std::size_t configured =
        ReadEnvSizeTClamped("ATPPS_SINGLE_RING_BUCKET_COUNT", kSingleRingBucketCount, 1U, 256U);
    return configured;
}

// This allows runtime tuning of how many single-ring steps are emitted into trace logs
std::size_t GetConfiguredSingleRingTraceSteps() {
    static const std::size_t configured =
        ReadEnvSizeTClamped("ATPPS_SINGLE_RING_TRACE_STEPS", kSingleRingTraceSteps, 0U, 1000000U);
    return configured;
}

// This allows runtime toggling of local topology-safe acceptance for single-ring pair collapses
int GetConfiguredSingleRingTopologyGuardMode() {
    static const int configured =
        ReadEnvIntClamped("ATPPS_SINGLE_RING_TOPO_GUARD_MODE", kSingleRingTopologyGuardMode, 0, 1);
    return configured;
}

// This allows runtime tuning of ring-size cap where the single-ring topology guard is applied
std::size_t GetConfiguredSingleRingTopologyGuardVertexLimit() {
    static const std::size_t configured =
        ReadEnvSizeTClamped(
            "ATPPS_SINGLE_RING_TOPO_GUARD_VERTEX_LIMIT",
            kSingleRingTopologyGuardVertexLimit,
            1000U,
            500000U
        );
    return configured;
}

// This allows runtime toggling of depth-conditioned multi-ring quota policy without changing default behavior
int GetConfiguredMultiRingDepthQuotaMode() {
    static const int configured =
        ReadEnvIntClamped("ATPPS_MULTI_RING_DEPTH_QUOTA_MODE", kMultiRingDepthQuotaMode, 0, 1);
    return configured;
}

// This allows runtime switching across legacy, local-normalized, and ring-normalized multi-ring scoring modes
int GetConfiguredMultiRingNormalizedScoreMode() {
    static const int configured =
        ReadEnvIntClamped("ATPPS_MULTI_RING_NORMALIZED_SCORE_MODE", kMultiRingNormalizedScoreMode, 0, 2);
    return configured;
}

// This chooses same-side AB/CD selection mode: 0 none, 1 all, 2 holes-only, 3 outers-only
int GetConfiguredApscSameSideMode() {
    static const int configured = []() {
        const char* explicitMode = std::getenv("ATPPS_APSC_SAME_SIDE_MODE");
        if (explicitMode != nullptr && explicitMode[0] != '\0') {
            return ReadEnvIntClamped("ATPPS_APSC_SAME_SIDE_MODE", 3, 0, 3);
        }

        const bool legacyFlip = ReadEnvBool("ATPPS_APSC_SAME_SIDE_FLIP", false);
        return legacyFlip ? 1 : 3;
    }();
    return configured;
}

// This maps configured same-side branch mode to one ring by containment-depth parity
bool ShouldFlipApscSameSideForContainmentDepth(const int containmentDepth) {
    const int mode = GetConfiguredApscSameSideMode();
    const bool outerLike = (containmentDepth % 2) == 0;

    if (mode == 1) {
        return true;
    }
    if (mode == 2) {
        return !outerLike;
    }
    if (mode == 3) {
        return outerLike;
    }
    return false;
}

// This allows independent single-ring control of same-side AB/CD branch selection during calibration
bool GetConfiguredSingleRingApscFlip() {
    static const bool configured = ReadEnvBool("ATPPS_SINGLE_RING_APSC_FLIP", false);
    return configured;
}

// This enforces deterministic candidate ordering for global best-choice selection
bool IsBetterCandidate(const RemovalCandidate& lhs, const RemovalCandidate& rhs) {
    if (LessWithTolerance(lhs.shapeError, rhs.shapeError)) {
        return true;
    }
    if (LessWithTolerance(rhs.shapeError, lhs.shapeError)) {
        return false;
    }

    if (LessWithTolerance(lhs.areaDelta, rhs.areaDelta)) {
        return true;
    }
    if (LessWithTolerance(rhs.areaDelta, lhs.areaDelta)) {
        return false;
    }

    if (lhs.ringId != rhs.ringId) {
        return lhs.ringId < rhs.ringId;
    }

    if (lhs.vertexIndex != rhs.vertexIndex) {
        return lhs.vertexIndex < rhs.vertexIndex;
    }

    if (LessWithTolerance(lhs.point.x, rhs.point.x)) {
        return true;
    }
    if (LessWithTolerance(rhs.point.x, lhs.point.x)) {
        return false;
    }

    if (LessWithTolerance(lhs.point.y, rhs.point.y)) {
        return true;
    }
    if (LessWithTolerance(rhs.point.y, lhs.point.y)) {
        return false;
    }

    return false;
}

// This checks scalar equality with a tight epsilon used by local topology predicates
bool NearlyEqual(const double lhs, const double rhs) {
    return std::abs(lhs - rhs) <= 1e-9;
}

// This checks coordinate equality so segment predicates can robustly reason about shared endpoints
bool PointsEqual(const atpps::Point& lhs, const atpps::Point& rhs) {
    return NearlyEqual(lhs.x, rhs.x) && NearlyEqual(lhs.y, rhs.y);
}

// This computes orientation sign for segment intersection tests with tolerance handling
int OrientationSign(const atpps::Point& a, const atpps::Point& b, const atpps::Point& c) {
    const double value = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
    if (value > 1e-9) {
        return 1;
    }
    if (value < -1e-9) {
        return -1;
    }
    return 0;
}

// This checks whether p lies on segment a-b under tolerance for collinear overlap cases
bool IsPointOnSegment(const atpps::Point& p, const atpps::Point& a, const atpps::Point& b) {
    if (OrientationSign(a, b, p) != 0) {
        return false;
    }

    const double minX = std::min(a.x, b.x);
    const double maxX = std::max(a.x, b.x);
    const double minY = std::min(a.y, b.y);
    const double maxY = std::max(a.y, b.y);
    return p.x >= (minX - 1e-9) && p.x <= (maxX + 1e-9)
        && p.y >= (minY - 1e-9) && p.y <= (maxY + 1e-9);
}

// This checks segment intersection including endpoint touches and collinear overlap
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

// This maps one original ring index to a deterministic arc bucket for distributed selection policies
std::size_t ComputeSingleRingBucketIndex(
    const std::size_t indexB,
    const std::size_t initialCount,
    const std::size_t bucketCount
) {
    if (bucketCount <= 1U || initialCount == 0U) {
        return 0U;
    }

    std::size_t bucket = (indexB * bucketCount) / initialCount;
    if (bucket >= bucketCount) {
        bucket = bucketCount - 1U;
    }
    return bucket;
}

// This computes signed area for a raw point loop used by local displacement geometry
double ComputeLoopSignedArea(const std::vector<atpps::Point>& loop) {
    if (loop.size() < 3U) {
        return 0.0;
    }

    double twiceArea = 0.0;
    for (std::size_t i = 0U; i < loop.size(); ++i) {
        const atpps::Point& a = loop[i];
        const atpps::Point& b = loop[(i + 1U) % loop.size()];
        twiceArea += a.x * b.y - b.x * a.y;
    }
    return 0.5 * twiceArea;
}

// This computes intersection between segment p1-p2 and infinite line a-b for clipping
atpps::Point IntersectSegmentWithLine(
    const atpps::Point& p1,
    const atpps::Point& p2,
    const atpps::Point& a,
    const atpps::Point& b
) {
    const double A1 = p2.y - p1.y;
    const double B1 = p1.x - p2.x;
    const double C1 = A1 * p1.x + B1 * p1.y;

    const double A2 = b.y - a.y;
    const double B2 = a.x - b.x;
    const double C2 = A2 * a.x + B2 * a.y;

    const double det = A1 * B2 - A2 * B1;
    if (std::abs(det) <= 1e-12) {
        return p2;
    }

    return atpps::Point{
        (B2 * C1 - B1 * C2) / det,
        (A1 * C2 - A2 * C1) / det
    };
}

// This clips one subject polygon by one oriented clip edge half-plane
std::vector<atpps::Point> ClipByEdge(
    const std::vector<atpps::Point>& subject,
    const atpps::Point& clipA,
    const atpps::Point& clipB,
    const bool clipIsCCW
) {
    std::vector<atpps::Point> output;
    if (subject.empty()) {
        return output;
    }

    auto inside = [&](const atpps::Point& p) {
        const double cross = (clipB.x - clipA.x) * (p.y - clipA.y) - (clipB.y - clipA.y) * (p.x - clipA.x);
        return clipIsCCW ? (cross >= -1e-12) : (cross <= 1e-12);
    };

    atpps::Point prev = subject.back();
    bool prevInside = inside(prev);

    for (const atpps::Point& cur : subject) {
        const bool curInside = inside(cur);

        if (curInside) {
            if (!prevInside) {
                output.push_back(IntersectSegmentWithLine(prev, cur, clipA, clipB));
            }
            output.push_back(cur);
        } else if (prevInside) {
            output.push_back(IntersectSegmentWithLine(prev, cur, clipA, clipB));
        }

        prev = cur;
        prevInside = curInside;
    }

    return output;
}

// This computes intersection area between arbitrary polygon and triangle via Sutherland-Hodgman clipping
double ComputePolygonTriangleIntersectionArea(
    const std::vector<atpps::Point>& polygon,
    const std::vector<atpps::Point>& triangle
) {
    if (polygon.size() < 3U || triangle.size() != 3U) {
        return 0.0;
    }

    std::vector<atpps::Point> clipped = polygon;
    const bool clipIsCCW = ComputeLoopSignedArea(triangle) >= 0.0;

    for (std::size_t i = 0U; i < 3U; ++i) {
        const atpps::Point& a = triangle[i];
        const atpps::Point& b = triangle[(i + 1U) % 3U];
        clipped = ClipByEdge(clipped, a, b, clipIsCCW);
        if (clipped.size() < 3U) {
            return 0.0;
        }
    }

    return std::abs(ComputeLoopSignedArea(clipped));
}

// This computes local areal displacement from replacing chain A-B-C-D with A-E-D
double ComputeLocalCollapseDisplacement(
    const atpps::Point& A,
    const atpps::Point& B,
    const atpps::Point& C,
    const atpps::Point& D,
    const atpps::Point& E
) {
    const std::vector<atpps::Point> oldLoop = {A, B, C, D};
    const std::vector<atpps::Point> newLoop = {A, E, D};

    const double oldArea = std::abs(ComputeLoopSignedArea(oldLoop));
    const double newArea = std::abs(ComputeLoopSignedArea(newLoop));
    const double intersectionArea = ComputePolygonTriangleIntersectionArea(oldLoop, newLoop);
    const double symmetricDiffArea = oldArea + newArea - 2.0 * intersectionArea;
    return (symmetricDiffArea < 0.0) ? 0.0 : symmetricDiffArea;
}

// This evaluates point-in-ring with boundary detection so containment checks stay robust
int PointInRing(const atpps::Point& p, const atpps::Ring& ring) {
    bool inside = false;
    const std::size_t n = ring.vertices.size();
    for (std::size_t i = 0U, j = n - 1U; i < n; j = i++) {
        const atpps::Point& a = ring.vertices[j];
        const atpps::Point& b = ring.vertices[i];

        if (IsPointOnSegment(p, a, b)) {
            return 0;
        }

        const bool intersects = ((a.y > p.y) != (b.y > p.y))
            && (p.x < (b.x - a.x) * (p.y - a.y) / ((b.y - a.y) + 1e-20) + a.x);
        if (intersects) {
            inside = !inside;
        }
    }

    return inside ? 1 : -1;
}

// This records strict containment relation for each ring pair in one deterministic matrix
std::vector<std::vector<bool>> BuildContainmentMatrix(const atpps::Polygon& polygon) {
    const std::size_t ringCount = polygon.rings.size();
    std::vector<std::vector<bool>> matrix(ringCount, std::vector<bool>(ringCount, false));

    for (std::size_t i = 0U; i < ringCount; ++i) {
        for (std::size_t j = 0U; j < ringCount; ++j) {
            if (i == j || polygon.rings[j].vertices.empty()) {
                continue;
            }

            const int relation = PointInRing(polygon.rings[j].vertices[0], polygon.rings[i]);
            matrix[i][j] = (relation > 0);
        }
    }

    return matrix;
}

// This ensures tentative simplification preserves original nesting between every ring pair
bool ContainmentPreserved(
    const atpps::Polygon& edited,
    const std::vector<std::vector<bool>>& baselineContainment
) {
    const auto editedContainment = BuildContainmentMatrix(edited);
    if (editedContainment.size() != baselineContainment.size()) {
        return false;
    }

    for (std::size_t i = 0U; i < baselineContainment.size(); ++i) {
        for (std::size_t j = 0U; j < baselineContainment[i].size(); ++j) {
            if (editedContainment[i][j] != baselineContainment[i][j]) {
                return false;
            }
        }
    }

    return true;
}

// This counts how many rings strictly contain one ring so we can derive nesting depth parity
int ComputeRingContainmentDepth(const std::vector<std::vector<bool>>& containmentMatrix, const std::size_t ringIndex) {
    int depth = 0;
    for (std::size_t container = 0U; container < containmentMatrix.size(); ++container) {
        if (containmentMatrix[container][ringIndex]) {
            ++depth;
        }
    }
    return depth;
}

// This marks even-depth rings as outer-like boundaries and odd-depth rings as hole-like boundaries
bool IsOuterLikeContainmentDepth(const int containmentDepth) {
    return (containmentDepth % 2) == 0;
}

// This keeps ring-role priority as a tie-break only when costs are near each other
bool CostsComparableForRingRolePriority(const double lhsCost, const double rhsCost) {
    const double smallerCost = std::min(lhsCost, rhsCost);
    const double largerCost = std::max(lhsCost, rhsCost);
    const double maxAllowedGap = GetConfiguredRingRoleComparableCostMultiplier() * (smallerCost + 1e-12);
    return (largerCost - smallerCost) <= maxAllowedGap;
}

// This keeps single-ring windowed selection age tie-break active only when candidate costs are close enough
bool CostsComparableForSingleRingSelection(const double lhsCost, const double rhsCost) {
    const double smallerCost = std::min(lhsCost, rhsCost);
    const double largerCost = std::max(lhsCost, rhsCost);
    const double maxAllowedGap = GetConfiguredSingleRingSelectionTieCostMultiplier() * (smallerCost + 1e-12);
    return (largerCost - smallerCost) <= maxAllowedGap;
}

// This computes signed side value for one point against directed line A->D
double SignedSideValue(const atpps::Point& A, const atpps::Point& D, const atpps::Point& point) {
    return (D.x - A.x) * (point.y - A.y) - (D.y - A.y) * (point.x - A.x);
}

// This maps signed side values to deterministic left or right codes
int SideCode(const double signedSideValue) {
    if (signedSideValue > 1e-12) {
        return 1;
    }
    if (signedSideValue < -1e-12) {
        return -1;
    }
    return 0;
}

// This builds one concrete point on implicit line a*x + b*y + c = 0 for side tests
atpps::Point BuildPointOnImplicitLine(const double a, const double b, const double c) {
    if (std::abs(b) > std::abs(a)) {
        return atpps::Point{0.0, -c / b};
    }

    if (std::abs(a) > 1e-12) {
        return atpps::Point{-c / a, 0.0};
    }

    return atpps::Point{0.0, 0.0};
}

// This intersects implicit line a*x + b*y + c = 0 with infinite line through segment start-end
bool IntersectImplicitLineWithTwoPointLine(
    const double a,
    const double b,
    const double c,
    const atpps::Point& start,
    const atpps::Point& end,
    atpps::Point& intersectionOut
) {
    const double dx = end.x - start.x;
    const double dy = end.y - start.y;
    const double denominator = a * dx + b * dy;
    if (std::abs(denominator) <= 1e-12) {
        return false;
    }

    const double t = -(a * start.x + b * start.y + c) / denominator;
    intersectionOut = atpps::Point{start.x + t * dx, start.y + t * dy};
    return true;
}

// This computes APSC replacement point using E-line and AB or CD intersection rules
atpps::Point ComputeReplacementPoint(
    const atpps::Point& A,
    const atpps::Point& B,
    const atpps::Point& C,
    const atpps::Point& D,
    const bool flipSameSideSelection
) {
    // This handles the collinear B-C-D singular case where removing C at B gives zero displacement
    if (std::abs(SignedTriangleArea(B, C, D)) <= 1e-12) {
        return B;
    }

    // This builds the paper's implicit E-line equation that preserves area exactly
    const double a = D.y - A.y;
    const double b = A.x - D.x;
    const double c =
        -B.y * A.x +
        (A.y - C.y) * B.x +
        (B.y - D.y) * C.x +
        C.y * D.x;

    atpps::Point intersectionWithAB;
    atpps::Point intersectionWithCD;
    const bool hasIntersectionWithAB = IntersectImplicitLineWithTwoPointLine(a, b, c, A, B, intersectionWithAB);
    const bool hasIntersectionWithCD = IntersectImplicitLineWithTwoPointLine(a, b, c, C, D, intersectionWithCD);

    // This keeps deterministic behavior for rare degenerate line-pair combinations
    if (!hasIntersectionWithAB && !hasIntersectionWithCD) {
        return atpps::Point{0.5 * (B.x + C.x), 0.5 * (B.y + C.y)};
    }
    if (!hasIntersectionWithAB) {
        return intersectionWithCD;
    }
    if (!hasIntersectionWithCD) {
        return intersectionWithAB;
    }

    const double sideBValue = SignedSideValue(A, D, B);
    const double sideCValue = SignedSideValue(A, D, C);
    const int sideB = SideCode(sideBValue);
    const int sideC = SideCode(sideCValue);

    const atpps::Point pointOnELine = BuildPointOnImplicitLine(a, b, c);
    const int sideELine = SideCode(SignedSideValue(A, D, pointOnELine));

    const double distanceB = std::abs(sideBValue);
    const double distanceC = std::abs(sideCValue);

    // This branch follows the paper rule for B and C on the same side of AD
    if (sideB == sideC) {
        if (distanceB > distanceC + 1e-12) {
            return flipSameSideSelection ? intersectionWithCD : intersectionWithAB;
        }
        if (distanceC > distanceB + 1e-12) {
            return flipSameSideSelection ? intersectionWithAB : intersectionWithCD;
        }

        if (sideB == sideELine) {
            return flipSameSideSelection ? intersectionWithCD : intersectionWithAB;
        }
        return flipSameSideSelection ? intersectionWithAB : intersectionWithCD;
    }

    // This branch follows the paper rule for B and C on opposite sides of AD
    if (sideB == sideELine) {
        return intersectionWithAB;
    }
    return intersectionWithCD;
}

// This creates one edited polygon state by replacing B/C with one replacement point
atpps::Polygon ApplyPairCollapse(
    const atpps::Polygon& polygon,
    const std::size_t ringIndex,
    const std::size_t indexB,
    const std::size_t indexC,
    const atpps::Point& replacement
) {
    atpps::Polygon edited = polygon;
    auto& vertices = edited.rings[ringIndex].vertices;

    // This path handles only forward-neighbor collapses in current candidate generation
    vertices[indexB] = replacement;
    vertices.erase(vertices.begin() + static_cast<std::ptrdiff_t>(indexC));
    return edited;
}

// This compares pair-collapse candidates deterministically by score and stable tie-breaks
bool IsBetterPairCollapseCandidate(const PairCollapseCandidate& lhs, const PairCollapseCandidate& rhs) {
    const bool lhsOuterLike = IsOuterLikeContainmentDepth(lhs.containmentDepth);
    const bool rhsOuterLike = IsOuterLikeContainmentDepth(rhs.containmentDepth);
    if (lhsOuterLike != rhsOuterLike && CostsComparableForRingRolePriority(lhs.cost, rhs.cost)) {
        return lhsOuterLike;
    }

    if (LessWithTolerance(lhs.cost, rhs.cost)) {
        return true;
    }
    if (LessWithTolerance(rhs.cost, lhs.cost)) {
        return false;
    }

    if (lhs.containmentDepth != rhs.containmentDepth) {
        return lhs.containmentDepth < rhs.containmentDepth;
    }

    if (lhs.ringId != rhs.ringId) {
        return lhs.ringId < rhs.ringId;
    }

    if (lhs.indexB != rhs.indexB) {
        return lhs.indexB < rhs.indexB;
    }

    if (LessWithTolerance(lhs.replacement.x, rhs.replacement.x)) {
        return true;
    }
    if (LessWithTolerance(rhs.replacement.x, lhs.replacement.x)) {
        return false;
    }

    if (LessWithTolerance(lhs.replacement.y, rhs.replacement.y)) {
        return true;
    }
    if (LessWithTolerance(rhs.replacement.y, lhs.replacement.y)) {
        return false;
    }

    return false;
}

// This compares pair-collapse keys so quota fill can avoid duplicate expansions
bool PairCollapseSameKey(const PairCollapseCandidate& lhs, const PairCollapseCandidate& rhs) {
    return lhs.ringIndex == rhs.ringIndex
        && lhs.indexB == rhs.indexB
        && lhs.indexC == rhs.indexC;
}

// This checks whether one candidate key is already present in a small selected set
bool PairCollapseKeyAlreadySelected(
    const std::vector<PairCollapseCandidate>& selected,
    const PairCollapseCandidate& candidate
) {
    for (const PairCollapseCandidate& existing : selected) {
        if (PairCollapseSameKey(existing, candidate)) {
            return true;
        }
    }
    return false;
}

// This collects top legal pair-collapse candidates for one polygon state
void CollectPairCollapseCandidates(
    const atpps::Polygon& polygon,
    const std::size_t maxCandidates,
    std::vector<PairCollapseCandidate>& candidatesOut
) {
    candidatesOut.clear();
    const int multiRingNormalizedScoreMode = GetConfiguredMultiRingNormalizedScoreMode();
    const double multiRingMovementWeight = GetConfiguredMultiRingMovementWeight();
    const auto baselineContainment = BuildContainmentMatrix(polygon);
    std::vector<int> containmentDepthByRing(polygon.rings.size(), 0);
    std::vector<double> ringLengthScaleSqByRing(polygon.rings.size(), 1.0);
    std::vector<double> ringAreaScaleByRing(polygon.rings.size(), 1.0);
    for (std::size_t ringIndex = 0U; ringIndex < polygon.rings.size(); ++ringIndex) {
        containmentDepthByRing[ringIndex] = ComputeRingContainmentDepth(baselineContainment, ringIndex);

        if (multiRingNormalizedScoreMode == 2) {
            const atpps::Ring& ring = polygon.rings[ringIndex];
            const double ringLengthScaleSq = ComputeRingLengthScaleSq(ring);
            ringLengthScaleSqByRing[ringIndex] = ringLengthScaleSq;
            ringAreaScaleByRing[ringIndex] = ComputeRingAreaScale(ring, ringLengthScaleSq);
        }
    }

    for (std::size_t ringIndex = 0U; ringIndex < polygon.rings.size(); ++ringIndex) {
        const atpps::Ring& ring = polygon.rings[ringIndex];
        const std::size_t n = ring.vertices.size();
        if (n <= 4U) {
            continue;
        }

        for (std::size_t indexB = 1U; indexB < n; ++indexB) {
            const std::size_t indexC = (indexB + 1U) % n;

            // This keeps anchor vertex 0 alive by never removing index 0 as C
            if (indexC == 0U) {
                continue;
            }

            const std::size_t indexA = (indexB + n - 1U) % n;
            const std::size_t indexD = (indexC + 1U) % n;

            const atpps::Point& A = ring.vertices[indexA];
            const atpps::Point& B = ring.vertices[indexB];
            const atpps::Point& C = ring.vertices[indexC];
            const atpps::Point& D = ring.vertices[indexD];

            const bool flipSameSide = ShouldFlipApscSameSideForContainmentDepth(containmentDepthByRing[ringIndex]);
            const atpps::Point E = ComputeReplacementPoint(A, B, C, D, flipSameSide);

            const atpps::Polygon edited = ApplyPairCollapse(polygon, ringIndex, indexB, indexC, E);
            std::string topologyError;
            if (!atpps::ValidatePolygonTopology(edited, topologyError)) {
                continue;
            }

            if (!ContainmentPreserved(edited, baselineContainment)) {
                continue;
            }

            const double dxB = E.x - B.x;
            const double dyB = E.y - B.y;
            const double dxC = E.x - C.x;
            const double dyC = E.y - C.y;
            const double movementCost = (dxB * dxB + dyB * dyB) + (dxC * dxC + dyC * dyC);
            const double localDisplacement = ComputeLocalCollapseDisplacement(A, B, C, D, E);

            double cost = 0.0;
            if (multiRingNormalizedScoreMode == 0) {
                cost = localDisplacement + multiRingMovementWeight * movementCost;
            } else if (multiRingNormalizedScoreMode == 1) {
                const double localLengthScaleSq = ComputeSingleRingLocalLengthScaleSq(A, B, C, D);
                const double localAreaScale = ComputeSingleRingLocalAreaScale(A, B, C, D, localLengthScaleSq);

                // This keeps multi-ring ranking dimensionless across datasets with different coordinate scales
                const double normalizedDisplacement = localDisplacement / localAreaScale;
                const double normalizedMovement = movementCost / localLengthScaleSq;
                cost = normalizedDisplacement + multiRingMovementWeight * normalizedMovement;
            } else {
                // This reuses one ring-level scale per candidate ring so ranking stays consistent across local neighborhoods
                const double normalizedDisplacement = localDisplacement / ringAreaScaleByRing[ringIndex];
                const double normalizedMovement = movementCost / ringLengthScaleSqByRing[ringIndex];
                cost = normalizedDisplacement + multiRingMovementWeight * normalizedMovement;
            }

            PairCollapseCandidate candidate;
            candidate.ringIndex = ringIndex;
            candidate.indexB = indexB;
            candidate.indexC = indexC;
            candidate.ringId = ring.ringId;
            candidate.containmentDepth = containmentDepthByRing[ringIndex];
            candidate.replacement = E;
            candidate.cost = cost;
            candidate.displacement = localDisplacement;

            candidatesOut.push_back(candidate);
        }
    }

    std::sort(
        candidatesOut.begin(),
        candidatesOut.end(),
        [](const PairCollapseCandidate& lhs, const PairCollapseCandidate& rhs) {
            return IsBetterPairCollapseCandidate(lhs, rhs);
        }
    );

    if (candidatesOut.size() <= maxCandidates) {
        return;
    }

    if (GetConfiguredMultiRingDepthQuotaMode() == 0 || maxCandidates <= 1U) {
        candidatesOut.resize(maxCandidates);
        return;
    }

    std::size_t outerLikeRingCount = 0U;
    std::size_t holeLikeRingCount = 0U;
    for (const int depth : containmentDepthByRing) {
        if (IsOuterLikeContainmentDepth(depth)) {
            ++outerLikeRingCount;
        } else {
            ++holeLikeRingCount;
        }
    }

    std::vector<PairCollapseCandidate> outerLikeCandidates;
    std::vector<PairCollapseCandidate> holeLikeCandidates;
    outerLikeCandidates.reserve(candidatesOut.size());
    holeLikeCandidates.reserve(candidatesOut.size());
    for (const PairCollapseCandidate& candidate : candidatesOut) {
        if (IsOuterLikeContainmentDepth(candidate.containmentDepth)) {
            outerLikeCandidates.push_back(candidate);
        } else {
            holeLikeCandidates.push_back(candidate);
        }
    }

    if (outerLikeCandidates.empty() || holeLikeCandidates.empty()) {
        candidatesOut.resize(maxCandidates);
        return;
    }

    const std::size_t parityDenominator = outerLikeRingCount + holeLikeRingCount;
    if (parityDenominator == 0U) {
        candidatesOut.resize(maxCandidates);
        return;
    }

    std::size_t desiredOuter = (maxCandidates * outerLikeRingCount + parityDenominator - 1U) / parityDenominator;
    if (desiredOuter == 0U) {
        desiredOuter = 1U;
    }
    if (desiredOuter >= maxCandidates) {
        desiredOuter = maxCandidates - 1U;
    }

    std::size_t desiredHole = maxCandidates - desiredOuter;
    if (desiredHole == 0U) {
        desiredHole = 1U;
        if (desiredOuter > 0U) {
            --desiredOuter;
        }
    }

    std::vector<PairCollapseCandidate> selected;
    selected.reserve(maxCandidates);

    for (std::size_t i = 0U; i < outerLikeCandidates.size() && selected.size() < desiredOuter; ++i) {
        selected.push_back(outerLikeCandidates[i]);
    }

    std::size_t holeAdded = 0U;
    for (std::size_t i = 0U; i < holeLikeCandidates.size() && holeAdded < desiredHole; ++i) {
        selected.push_back(holeLikeCandidates[i]);
        ++holeAdded;
    }

    if (selected.size() < maxCandidates) {
        for (const PairCollapseCandidate& candidate : candidatesOut) {
            if (selected.size() >= maxCandidates) {
                break;
            }

            if (PairCollapseKeyAlreadySelected(selected, candidate)) {
                continue;
            }

            selected.push_back(candidate);
        }
    }

    std::sort(
        selected.begin(),
        selected.end(),
        [](const PairCollapseCandidate& lhs, const PairCollapseCandidate& rhs) {
            return IsBetterPairCollapseCandidate(lhs, rhs);
        }
    );

    if (selected.size() > maxCandidates) {
        selected.resize(maxCandidates);
    }

    candidatesOut = std::move(selected);
}

// This validates the new stitched edge against all unaffected edges of the edited ring
bool NewEdgeIsRingSafe(const atpps::Ring& ring, const std::size_t vertexIndex) {
    const std::size_t count = ring.vertices.size();
    const std::size_t prevIndex = (vertexIndex + count - 1U) % count;
    const std::size_t nextIndex = (vertexIndex + 1U) % count;
    const std::size_t beforePrev = (prevIndex + count - 1U) % count;

    const atpps::Point& newA = ring.vertices[prevIndex];
    const atpps::Point& newB = ring.vertices[nextIndex];

    // This rejects collapsing to a zero-length edge which would violate ring primitive shape
    if (PointsEqual(newA, newB)) {
        return false;
    }

    for (std::size_t edge = 0U; edge < count; ++edge) {
        // This skips the two removed edges and the two edges adjacent to the stitched edge
        if (edge == prevIndex || edge == vertexIndex || edge == beforePrev || edge == nextIndex) {
            continue;
        }

        const atpps::Point& otherA = ring.vertices[edge];
        const atpps::Point& otherB = ring.vertices[(edge + 1U) % count];
        if (SegmentsIntersect(newA, newB, otherA, otherB)) {
            return false;
        }
    }

    return true;
}

// This validates the stitched edge against every edge of every other ring
bool NewEdgeIsInterRingSafe(
    const atpps::Polygon& polygon,
    const std::size_t ringIndex,
    const std::size_t vertexIndex
) {
    const atpps::Ring& ring = polygon.rings[ringIndex];
    const std::size_t count = ring.vertices.size();
    const std::size_t prevIndex = (vertexIndex + count - 1U) % count;
    const std::size_t nextIndex = (vertexIndex + 1U) % count;

    const atpps::Point& newA = ring.vertices[prevIndex];
    const atpps::Point& newB = ring.vertices[nextIndex];

    for (std::size_t otherRing = 0U; otherRing < polygon.rings.size(); ++otherRing) {
        if (otherRing == ringIndex) {
            continue;
        }

        const atpps::Ring& other = polygon.rings[otherRing];
        for (std::size_t edge = 0U; edge < other.vertices.size(); ++edge) {
            const atpps::Point& otherA = other.vertices[edge];
            const atpps::Point& otherB = other.vertices[(edge + 1U) % other.vertices.size()];
            if (SegmentsIntersect(newA, newB, otherA, otherB)) {
                return false;
            }
        }
    }

    return true;
}

// This runs fast local legality checks that are equivalent to full validation for one removal edit
bool IsRemovalLocallyLegal(const atpps::Polygon& polygon, const RemovalCandidate& candidate) {
    const atpps::Ring& ring = polygon.rings[candidate.ringIndex];

    // This keeps each ring polygonal under simplification
    if (ring.vertices.size() <= 3U) {
        return false;
    }

    // This fast path keeps large single-ring lake datasets tractable under strict runtime caps
    // These inputs have no inter-ring topology constraints, so we prioritize deterministic progress
    if (polygon.rings.size() == 1U) {
        return true;
    }

    if (!NewEdgeIsRingSafe(ring, candidate.vertexIndex)) {
        return false;
    }

    if (!NewEdgeIsInterRingSafe(polygon, candidate.ringIndex, candidate.vertexIndex)) {
        return false;
    }

    return true;
}

// This checks whether one single-ring pair collapse keeps the ring locally simple under current linked topology state
bool IsSingleRingPairCollapseTopologicallySafe(
    const atpps::Ring& ring,
    const std::vector<std::size_t>& prev,
    const std::vector<std::size_t>& next,
    const std::vector<bool>& alive,
    const std::size_t indexB,
    const atpps::Point& replacement
) {
    if (indexB == 0U || !alive[indexB]) {
        return false;
    }

    const std::size_t indexC = next[indexB];
    if (indexC == 0U || !alive[indexC]) {
        return false;
    }

    const std::size_t indexA = prev[indexB];
    const std::size_t indexD = next[indexC];
    if (!alive[indexA] || !alive[indexD]) {
        return false;
    }

    const atpps::Point& pointA = ring.vertices[indexA];
    const atpps::Point& pointD = ring.vertices[indexD];

    // This rejects zero-length stitched edges so accepted collapses stay polygonal
    if (PointsEqual(pointA, replacement) || PointsEqual(replacement, pointD) || PointsEqual(pointA, pointD)) {
        return false;
    }

    for (std::size_t edgeStart = 0U; edgeStart < alive.size(); ++edgeStart) {
        if (!alive[edgeStart]) {
            continue;
        }

        const std::size_t edgeEnd = next[edgeStart];
        if (!alive[edgeEnd]) {
            continue;
        }

        // This skips the three removed edges and both stitched-edge adjacency endpoints
        if (
            edgeStart == indexA || edgeStart == indexB || edgeStart == indexC || edgeStart == indexD ||
            edgeEnd == indexA || edgeEnd == indexB || edgeEnd == indexC || edgeEnd == indexD
        ) {
            continue;
        }

        const atpps::Point& otherA = ring.vertices[edgeStart];
        const atpps::Point& otherB = ring.vertices[edgeEnd];

        if (SegmentsIntersect(pointA, replacement, otherA, otherB)) {
            return false;
        }

        if (SegmentsIntersect(replacement, pointD, otherA, otherB)) {
            return false;
        }
    }

    return true;
}

// This simplifies one single ring in O(n log n) style by lazy heap updates around edited neighborhoods
void SimplifySingleRingFast(atpps::Ring& ring, const std::size_t targetVertices, double& totalDisplacement) {
    const std::size_t initialCount = ring.vertices.size();
    if (initialCount <= targetVertices || initialCount <= 3U) {
        return;
    }

    std::vector<std::size_t> prev(initialCount);
    std::vector<std::size_t> next(initialCount);
    std::vector<std::size_t> generation(initialCount, 0U);
    std::vector<std::size_t> lastTouchStep(initialCount, std::numeric_limits<std::size_t>::max());
    std::vector<bool> alive(initialCount, true);
    std::size_t aliveCount = initialCount;
    std::size_t collapseStep = 0U;

    // This snapshots staged-policy knobs once so one run stays internally consistent
    const int stagedPolicyMode = GetConfiguredSingleRingStagedPolicyMode();
    const int singleRingScoreMode = GetConfiguredSingleRingScoreMode();
    const int singleRingDispersionMode = GetConfiguredSingleRingDispersionMode();
    const int singleRingSelectionMode = GetConfiguredSingleRingSelectionMode();
    const double singleRingAdaptiveSpreadThreshold = GetConfiguredSingleRingAdaptiveModeCostSpreadThreshold();
    const std::size_t singleRingAdaptiveSizeThreshold = GetConfiguredSingleRingAdaptiveSizeThreshold();
    const std::size_t singleRingSelectionWindow = GetConfiguredSingleRingSelectionWindow();
    const int singleRingWindowPoolMode = GetConfiguredSingleRingWindowPoolMode();
    const std::size_t singleRingWindowPoolExpandFactor = GetConfiguredSingleRingWindowPoolExpandFactor();
    const std::size_t singleRingWindowPoolMinCount = GetConfiguredSingleRingWindowPoolMinCount();
    const int singleRingProgressScheduleMode = GetConfiguredSingleRingProgressScheduleMode();
    const int singleRingWindowPrefilterMode = GetConfiguredSingleRingWindowPrefilterMode();
    const double singleRingWindowPrefilterPercentile = GetConfiguredSingleRingWindowPrefilterPercentile();
    const std::size_t singleRingWindowPrefilterMinCount = GetConfiguredSingleRingWindowPrefilterMinCount();
    const int singleRingWindowCostPrefilterMode = GetConfiguredSingleRingWindowCostPrefilterMode();
    const double singleRingWindowCostPrefilterPercentile = GetConfiguredSingleRingWindowCostPrefilterPercentile();
    const std::size_t singleRingWindowCostPrefilterMinCount = GetConfiguredSingleRingWindowCostPrefilterMinCount();
    const int singleRingWindowAgePrefilterMode = GetConfiguredSingleRingWindowAgePrefilterMode();
    const double singleRingWindowAgePrefilterPercentile = GetConfiguredSingleRingWindowAgePrefilterPercentile();
    const std::size_t singleRingWindowAgePrefilterMinCount = GetConfiguredSingleRingWindowAgePrefilterMinCount();
    const std::size_t singleRingBucketCount = GetConfiguredSingleRingBucketCount();
    const std::size_t singleRingDispersionWindow = GetConfiguredSingleRingDispersionWindow();
    const std::size_t singleRingTraceSteps = GetConfiguredSingleRingTraceSteps();
    const double singleRingDispersionPenaltyMultiplier = GetConfiguredSingleRingDispersionPenaltyMultiplier();
    const double featurePenaltyMultiplier = GetConfiguredSingleRingFeaturePenaltyMultiplier();
    const double singleRingMovementWeight = GetConfiguredSingleRingMovementWeight();
    const double singleRingMovementStartFactor = GetConfiguredSingleRingMovementStartFactor();
    const double singleRingMovementEndFactor = GetConfiguredSingleRingMovementEndFactor();
    const int singleRingTopologyGuardMode = GetConfiguredSingleRingTopologyGuardMode();
    const std::size_t singleRingTopologyGuardVertexLimit = GetConfiguredSingleRingTopologyGuardVertexLimit();
    const bool singleRingTopologyGuardEnabled =
        singleRingTopologyGuardMode != 0 && initialCount <= singleRingTopologyGuardVertexLimit;
    std::vector<std::size_t> bucketLastTouch(singleRingBucketCount, std::numeric_limits<std::size_t>::max());
    std::size_t topologyGuardRejectCount = 0U;

    if (TraceEnabled()) {
        std::ostringstream line;
        line << "TRACE_SINGLE_BEGIN"
             << " initial_vertices=" << initialCount
             << " target=" << targetVertices
             << " score_mode=" << singleRingScoreMode
             << " selection_mode=" << singleRingSelectionMode
             << " adaptive_spread_threshold=" << singleRingAdaptiveSpreadThreshold
             << " adaptive_size_threshold=" << singleRingAdaptiveSizeThreshold
             << " selection_window=" << singleRingSelectionWindow
             << " pool_mode=" << singleRingWindowPoolMode
             << " pool_expand_factor=" << singleRingWindowPoolExpandFactor
             << " pool_min_count=" << singleRingWindowPoolMinCount
             << " progress_schedule_mode=" << singleRingProgressScheduleMode
             << " prefilter_mode=" << singleRingWindowPrefilterMode
             << " prefilter_percentile=" << singleRingWindowPrefilterPercentile
             << " prefilter_min_count=" << singleRingWindowPrefilterMinCount
             << " cost_prefilter_mode=" << singleRingWindowCostPrefilterMode
             << " cost_prefilter_percentile=" << singleRingWindowCostPrefilterPercentile
             << " cost_prefilter_min_count=" << singleRingWindowCostPrefilterMinCount
             << " age_prefilter_mode=" << singleRingWindowAgePrefilterMode
             << " age_prefilter_percentile=" << singleRingWindowAgePrefilterPercentile
             << " age_prefilter_min_count=" << singleRingWindowAgePrefilterMinCount
             << " bucket_count=" << singleRingBucketCount
             << " staged_mode=" << stagedPolicyMode
             << " dispersion_mode=" << singleRingDispersionMode
             << " move_weight=" << singleRingMovementWeight
             << " move_start_factor=" << singleRingMovementStartFactor
             << " move_end_factor=" << singleRingMovementEndFactor
             << " topo_guard_mode=" << singleRingTopologyGuardMode
             << " topo_guard_limit=" << singleRingTopologyGuardVertexLimit
             << " topo_guard_enabled=" << (singleRingTopologyGuardEnabled ? 1 : 0)
             << " trace_steps=" << singleRingTraceSteps;
        TraceWriteLine(line.str());
        TraceWriteLine("TRACE_SINGLE_INPUT_" + FormatRingForTrace(ring));
    }

    for (std::size_t i = 0U; i < initialCount; ++i) {
        prev[i] = (i + initialCount - 1U) % initialCount;
        next[i] = (i + 1U) % initialCount;
    }

    auto buildPairCandidate = [&](const std::size_t indexB, SingleRingPairCandidate& outCandidate) {
        if (indexB == 0U || !alive[indexB]) {
            return false;
        }

        const std::size_t indexC = next[indexB];
        if (indexC == 0U || !alive[indexC]) {
            return false;
        }

        const std::size_t indexA = prev[indexB];
        const std::size_t indexD = next[indexC];
        if (!alive[indexA] || !alive[indexD]) {
            return false;
        }

        const atpps::Point& A = ring.vertices[indexA];
        const atpps::Point& B = ring.vertices[indexB];
        const atpps::Point& C = ring.vertices[indexC];
        const atpps::Point& D = ring.vertices[indexD];

        const atpps::Point E = ComputeReplacementPoint(A, B, C, D, GetConfiguredSingleRingApscFlip());
        const double displacement = ComputeLocalCollapseDisplacement(A, B, C, D, E);

        const double dxB = E.x - B.x;
        const double dyB = E.y - B.y;
        const double dxC = E.x - C.x;
        const double dyC = E.y - C.y;
        const double movementCost = (dxB * dxB + dyB * dyB) + (dxC * dxC + dyC * dyC);

        const double localLengthScaleSq = ComputeSingleRingLocalLengthScaleSq(A, B, C, D);
        const double localAreaScale = ComputeSingleRingLocalAreaScale(A, B, C, D, localLengthScaleSq);

        // This makes both score terms dimensionless so one global lambda behaves consistently across coordinate scales
        const double normalizedDisplacement = displacement / localAreaScale;
        const double normalizedMovement = movementCost / localLengthScaleSq;

        outCandidate.indexB = indexB;
        outCandidate.generationB = generation[indexB];
        outCandidate.generationC = generation[indexC];
        outCandidate.bucketIndex = ComputeSingleRingBucketIndex(indexB, initialCount, singleRingBucketCount);
        outCandidate.replacement = E;
        outCandidate.displacement = displacement;

        const double progress = ComputeSingleRingProgress(initialCount, aliveCount, targetVertices);

        // This schedules movement influence across progress so early and late collapse behavior can be calibrated independently
        const double movementFactor =
            singleRingMovementStartFactor +
            (singleRingMovementEndFactor - singleRingMovementStartFactor) * progress;
        const double effectiveSingleRingMovementWeight = singleRingMovementWeight * movementFactor;

        double baseCost = normalizedDisplacement + effectiveSingleRingMovementWeight * normalizedMovement;
        if (singleRingScoreMode == 1) {
            // This alternate mode prefers minimal movement first and uses displacement as a deterministic tie signal
            baseCost = normalizedMovement + effectiveSingleRingMovementWeight * normalizedDisplacement;
        } else if (singleRingScoreMode == 2) {
            // This alternate mode uses movement-first early and displacement-first late to vary sequence structure by progress
            const bool earlyStage = progress < 0.5;
            if (earlyStage) {
                baseCost = normalizedMovement + 1e-3 * normalizedDisplacement;
            } else {
                baseCost = normalizedDisplacement + 1e-3 * normalizedMovement;
            }
        }

        double candidateCost = baseCost;
        const double cornerSharpnessB = ComputeCornerSharpness(A, B, C);
        const double cornerSharpnessC = ComputeCornerSharpness(B, C, D);
        const double cornerSharpness = std::max(cornerSharpnessB, cornerSharpnessC);

        if (stagedPolicyMode != 0) {
            const double earlyStageWeight = (1.0 - progress) * featurePenaltyMultiplier;

            // This keeps high-corner collapses more expensive early and gradually relaxes near target
            candidateCost = candidateCost * (1.0 + earlyStageWeight * cornerSharpness);
        }

        std::size_t minRecentAge = std::numeric_limits<std::size_t>::max();

        const auto considerRecentAge = [&](const std::size_t index) {
            const std::size_t touchedAt = lastTouchStep[index];
            if (touchedAt == std::numeric_limits<std::size_t>::max() || touchedAt > collapseStep) {
                return;
            }

            const std::size_t age = collapseStep - touchedAt;
            if (age < minRecentAge) {
                minRecentAge = age;
            }
        };

        considerRecentAge(indexA);
        considerRecentAge(indexB);
        considerRecentAge(indexC);
        considerRecentAge(indexD);

        outCandidate.recentAge = minRecentAge;
        outCandidate.cornerSharpness = cornerSharpness;

        const bool hasRecentTouchAge = (minRecentAge != std::numeric_limits<std::size_t>::max());
        const bool hasRecentNeighborhoodTouch =
            hasRecentTouchAge && minRecentAge < singleRingDispersionWindow;

        if (singleRingDispersionMode != 0) {
            // This mode enforces early-stage spatial spread by skipping locally clustered collapses when alternatives exist
            if (singleRingDispersionMode == 1) {
                if (hasRecentNeighborhoodTouch && progress < 0.8 && aliveCount > targetVertices + 2U) {
                    return false;
                }
            } else {
                // This mode keeps all candidates but raises recent-neighborhood costs so clustered removals are delayed
                if (hasRecentNeighborhoodTouch) {
                    const double window = static_cast<double>(singleRingDispersionWindow);
                    const double age = static_cast<double>(minRecentAge);
                    const double recentness = (window - age) / window;
                    const double earlyStageWeight = (1.0 - progress);
                    const double multiplier = 1.0 + singleRingDispersionPenaltyMultiplier * recentness * earlyStageWeight;
                    candidateCost *= multiplier;
                }
            }
        }

        outCandidate.cost = candidateCost;
        return true;
    };

    std::priority_queue<SingleRingPairCandidate, std::vector<SingleRingPairCandidate>, SingleRingPairCandidateWorse> heap;
    for (std::size_t i = 1U; i < initialCount; ++i) {
        SingleRingPairCandidate candidate;
        if (buildPairCandidate(i, candidate)) {
            heap.push(candidate);
        }
    }

    const auto candidateIsCurrent = [&](const SingleRingPairCandidate& candidate) {
        const std::size_t indexB = candidate.indexB;
        if (!alive[indexB] || indexB == 0U) {
            return false;
        }

        const std::size_t indexC = next[indexB];
        if (indexC == 0U || !alive[indexC]) {
            return false;
        }

        if (candidate.generationB != generation[indexB] || candidate.generationC != generation[indexC]) {
            return false;
        }

        const std::size_t indexA = prev[indexB];
        const std::size_t indexD = next[indexC];
        if (!alive[indexA] || !alive[indexD]) {
            return false;
        }

        return true;
    };

    while (aliveCount > targetVertices && aliveCount > 3U && !heap.empty()) {
        SingleRingPairCandidate top;
        bool foundCandidate = false;

        if (singleRingSelectionMode == 0) {
            // This path keeps legacy heap-best behavior so baseline results remain reproducible
            while (!heap.empty()) {
                SingleRingPairCandidate candidate = heap.top();
                heap.pop();
                if (!candidateIsCurrent(candidate)) {
                    continue;
                }

                top = candidate;
                foundCandidate = true;
                break;
            }
        } else {
            // This path compares a small deterministic window of valid candidates to avoid over-clustered local choices
            const double remainingRatio = ComputeSingleRingRemainingRatio(initialCount, aliveCount, targetVertices);

            double scheduledPoolExpandFactor = static_cast<double>(singleRingWindowPoolExpandFactor);
            double scheduledWindowPrefilterPercentile = singleRingWindowPrefilterPercentile;
            double scheduledWindowCostPrefilterPercentile = singleRingWindowCostPrefilterPercentile;
            int effectiveWindowPoolMode = singleRingWindowPoolMode;

            if (singleRingProgressScheduleMode != 0) {
                scheduledPoolExpandFactor = SelectSingleRingProgressBandValue(
                    remainingRatio,
                    kSingleRingWindowPoolExpandEarly,
                    kSingleRingWindowPoolExpandMid,
                    kSingleRingWindowPoolExpandLate
                );

                // This pool-only schedule keeps broader spatial diversity early and shifts to strict cost-first near target
                if (singleRingProgressScheduleMode == 1 && remainingRatio <= (1.0 / 3.0)) {
                    effectiveWindowPoolMode = 0;
                }

                // This optional legacy schedule also phases displacement and cost percentiles by progress band
                if (singleRingProgressScheduleMode == 2) {
                    scheduledWindowPrefilterPercentile = SelectSingleRingProgressBandValue(
                        remainingRatio,
                        kSingleRingWindowPrefilterEarly,
                        kSingleRingWindowPrefilterMid,
                        kSingleRingWindowPrefilterLate
                    );
                    scheduledWindowCostPrefilterPercentile = SelectSingleRingProgressBandValue(
                        remainingRatio,
                        kSingleRingWindowCostPrefilterEarly,
                        kSingleRingWindowCostPrefilterMid,
                        kSingleRingWindowCostPrefilterLate
                    );
                }
            }

            std::size_t effectivePoolExpandFactor = static_cast<std::size_t>(std::llround(scheduledPoolExpandFactor));
            if (effectivePoolExpandFactor < 1U) {
                effectivePoolExpandFactor = 1U;
            }
            if (effectivePoolExpandFactor > 16U) {
                effectivePoolExpandFactor = 16U;
            }

            const std::size_t pooledTarget = (effectiveWindowPoolMode == 0)
                ? singleRingSelectionWindow
                : std::min<std::size_t>(
                    256U,
                    std::max(
                        singleRingSelectionWindow,
                        singleRingSelectionWindow * effectivePoolExpandFactor
                    )
                );

            std::vector<SingleRingPairCandidate> candidatePool;
            candidatePool.reserve(pooledTarget);

            while (!heap.empty() && candidatePool.size() < pooledTarget) {
                SingleRingPairCandidate candidate = heap.top();
                heap.pop();
                if (!candidateIsCurrent(candidate)) {
                    continue;
                }

                candidatePool.push_back(candidate);
            }

            std::vector<SingleRingPairCandidate> windowCandidates;
            windowCandidates.reserve(std::min(singleRingSelectionWindow, candidatePool.size()));

            if (!candidatePool.empty()) {
                std::vector<bool> selected(candidatePool.size(), false);

                const bool usePoolComposition =
                    effectiveWindowPoolMode != 0
                    && candidatePool.size() >= singleRingWindowPoolMinCount
                    && singleRingBucketCount > 1U;

                if (usePoolComposition) {
                    std::vector<std::vector<std::size_t>> indicesByBucket(singleRingBucketCount);
                    for (std::size_t i = 0U; i < candidatePool.size(); ++i) {
                        indicesByBucket[candidatePool[i].bucketIndex].push_back(i);
                    }

                    std::vector<std::size_t> bucketCursor(singleRingBucketCount, 0U);
                    while (windowCandidates.size() < singleRingSelectionWindow) {
                        bool addedAny = false;
                        for (std::size_t bucket = 0U; bucket < singleRingBucketCount; ++bucket) {
                            std::vector<std::size_t>& indices = indicesByBucket[bucket];
                            std::size_t& cursor = bucketCursor[bucket];
                            if (cursor >= indices.size()) {
                                continue;
                            }

                            const std::size_t poolIndex = indices[cursor++];
                            if (selected[poolIndex]) {
                                continue;
                            }

                            windowCandidates.push_back(candidatePool[poolIndex]);
                            selected[poolIndex] = true;
                            addedAny = true;

                            if (windowCandidates.size() >= singleRingSelectionWindow) {
                                break;
                            }
                        }

                        if (!addedAny) {
                            break;
                        }
                    }
                }

                for (std::size_t i = 0U; i < candidatePool.size() && windowCandidates.size() < singleRingSelectionWindow; ++i) {
                    if (selected[i]) {
                        continue;
                    }

                    windowCandidates.push_back(candidatePool[i]);
                    selected[i] = true;
                }

                for (std::size_t i = 0U; i < candidatePool.size(); ++i) {
                    if (!selected[i]) {
                        heap.push(candidatePool[i]);
                    }
                }

                if (TraceEnabled() && collapseStep < singleRingTraceSteps) {
                    std::ostringstream line;
                    line << "TRACE_SINGLE_POOL"
                         << " step=" << collapseStep
                            << " mode=" << effectiveWindowPoolMode
                            << " configured_mode=" << singleRingWindowPoolMode
                         << " remaining_ratio=" << remainingRatio
                         << " expand_factor=" << effectivePoolExpandFactor
                         << " pooled_target=" << pooledTarget
                         << " pooled_count=" << candidatePool.size()
                         << " window_count=" << windowCandidates.size();
                    TraceWriteLine(line.str());
                }
            }

            if (TraceEnabled() && collapseStep < singleRingTraceSteps) {
                std::ostringstream line;
                line << "TRACE_SINGLE_WINDOW"
                     << " step=" << collapseStep
                     << " alive=" << aliveCount
                     << " mode=" << singleRingSelectionMode
                     << " count=" << windowCandidates.size();
                TraceWriteLine(line.str());

                for (std::size_t i = 0U; i < windowCandidates.size(); ++i) {
                    const SingleRingPairCandidate& candidate = windowCandidates[i];
                    std::ostringstream candidateLine;
                    candidateLine << "TRACE_SINGLE_WINDOW_CAND"
                                  << " step=" << collapseStep
                                  << " rank=" << i
                                  << " b=" << candidate.indexB
                                  << " bucket=" << candidate.bucketIndex
                                  << " recent_age=" << candidate.recentAge
                                  << " corner=" << candidate.cornerSharpness
                                  << " cost=" << candidate.cost
                                  << " disp=" << candidate.displacement
                                  << " E=(" << FormatPointForTrace(candidate.replacement) << ')';
                    TraceWriteLine(candidateLine.str());
                }
            }

            if (!windowCandidates.empty()) {
                std::vector<std::size_t> selectionIndices;
                selectionIndices.reserve(windowCandidates.size());
                for (std::size_t i = 0U; i < windowCandidates.size(); ++i) {
                    selectionIndices.push_back(i);
                }

                if (
                    singleRingWindowPrefilterMode != 0
                    && windowCandidates.size() >= singleRingWindowPrefilterMinCount
                ) {
                    std::vector<double> sortedDisplacements;
                    sortedDisplacements.reserve(windowCandidates.size());
                    for (const SingleRingPairCandidate& candidate : windowCandidates) {
                        sortedDisplacements.push_back(candidate.displacement);
                    }
                    std::sort(sortedDisplacements.begin(), sortedDisplacements.end());

                    const double scaledRank =
                        scheduledWindowPrefilterPercentile * static_cast<double>(sortedDisplacements.size() - 1U);
                    std::size_t rankIndex = static_cast<std::size_t>(scaledRank);
                    if (rankIndex >= sortedDisplacements.size()) {
                        rankIndex = sortedDisplacements.size() - 1U;
                    }

                    const double displacementThreshold = sortedDisplacements[rankIndex];

                    std::vector<std::size_t> filteredIndices;
                    filteredIndices.reserve(selectionIndices.size());
                    for (const std::size_t index : selectionIndices) {
                        if (windowCandidates[index].displacement <= displacementThreshold + 1e-12) {
                            filteredIndices.push_back(index);
                        }
                    }

                    if (!filteredIndices.empty()) {
                        selectionIndices = std::move(filteredIndices);
                    }

                    if (TraceEnabled() && collapseStep < singleRingTraceSteps) {
                        std::ostringstream line;
                        line << "TRACE_SINGLE_PREFILTER"
                             << " step=" << collapseStep
                             << " mode=" << singleRingWindowPrefilterMode
                                << " percentile=" << scheduledWindowPrefilterPercentile
                             << " threshold=" << displacementThreshold
                             << " before=" << windowCandidates.size()
                             << " after=" << selectionIndices.size();
                        TraceWriteLine(line.str());
                    }
                }

                if (
                    singleRingWindowCostPrefilterMode != 0
                    && selectionIndices.size() >= singleRingWindowCostPrefilterMinCount
                ) {
                    std::vector<double> sortedCosts;
                    sortedCosts.reserve(selectionIndices.size());
                    for (const std::size_t index : selectionIndices) {
                        sortedCosts.push_back(windowCandidates[index].cost);
                    }
                    std::sort(sortedCosts.begin(), sortedCosts.end());

                    const double scaledRank =
                        scheduledWindowCostPrefilterPercentile * static_cast<double>(sortedCosts.size() - 1U);
                    std::size_t rankIndex = static_cast<std::size_t>(scaledRank);
                    if (rankIndex >= sortedCosts.size()) {
                        rankIndex = sortedCosts.size() - 1U;
                    }

                    const double costThreshold = sortedCosts[rankIndex];

                    std::vector<std::size_t> filteredByCost;
                    filteredByCost.reserve(selectionIndices.size());
                    for (const std::size_t index : selectionIndices) {
                        if (windowCandidates[index].cost <= costThreshold + 1e-12) {
                            filteredByCost.push_back(index);
                        }
                    }

                    if (!filteredByCost.empty()) {
                        selectionIndices = std::move(filteredByCost);
                    }

                    if (TraceEnabled() && collapseStep < singleRingTraceSteps) {
                        std::ostringstream line;
                        line << "TRACE_SINGLE_COST_PREFILTER"
                             << " step=" << collapseStep
                             << " mode=" << singleRingWindowCostPrefilterMode
                                << " percentile=" << scheduledWindowCostPrefilterPercentile
                             << " threshold=" << costThreshold
                             << " after=" << selectionIndices.size();
                        TraceWriteLine(line.str());
                    }
                }

                if (
                    singleRingWindowAgePrefilterMode != 0
                    && selectionIndices.size() >= singleRingWindowAgePrefilterMinCount
                ) {
                    double minPoolCost = std::numeric_limits<double>::infinity();
                    double maxPoolCost = 0.0;
                    for (const std::size_t index : selectionIndices) {
                        minPoolCost = std::min(minPoolCost, windowCandidates[index].cost);
                        maxPoolCost = std::max(maxPoolCost, windowCandidates[index].cost);
                    }

                    const bool comparablePool = CostsComparableForSingleRingSelection(minPoolCost, maxPoolCost);
                    const bool applyAgePrefilter =
                        singleRingWindowAgePrefilterMode == 1 ||
                        (singleRingWindowAgePrefilterMode == 2 && comparablePool);

                    if (applyAgePrefilter) {
                        const std::size_t poolSize = selectionIndices.size();
                        std::size_t keepCount = static_cast<std::size_t>(
                            std::ceil(singleRingWindowAgePrefilterPercentile * static_cast<double>(poolSize))
                        );
                        if (keepCount < 1U) {
                            keepCount = 1U;
                        }
                        if (keepCount > poolSize) {
                            keepCount = poolSize;
                        }

                        std::vector<std::size_t> orderedByAge = selectionIndices;
                        std::stable_sort(
                            orderedByAge.begin(),
                            orderedByAge.end(),
                            [&](const std::size_t lhsIndex, const std::size_t rhsIndex) {
                                const SingleRingPairCandidate& lhs = windowCandidates[lhsIndex];
                                const SingleRingPairCandidate& rhs = windowCandidates[rhsIndex];

                                if (lhs.recentAge != rhs.recentAge) {
                                    return lhs.recentAge > rhs.recentAge;
                                }

                                if (LessWithTolerance(lhs.cost, rhs.cost)) {
                                    return true;
                                }
                                if (LessWithTolerance(rhs.cost, lhs.cost)) {
                                    return false;
                                }

                                return lhs.indexB < rhs.indexB;
                            }
                        );

                        std::vector<std::size_t> filteredByAge;
                        filteredByAge.reserve(keepCount);
                        for (std::size_t i = 0U; i < keepCount; ++i) {
                            filteredByAge.push_back(orderedByAge[i]);
                        }

                        selectionIndices = std::move(filteredByAge);

                        if (TraceEnabled() && collapseStep < singleRingTraceSteps) {
                            std::ostringstream line;
                            line << "TRACE_SINGLE_AGE_PREFILTER"
                                 << " step=" << collapseStep
                                 << " mode=" << singleRingWindowAgePrefilterMode
                                 << " percentile=" << singleRingWindowAgePrefilterPercentile
                                 << " keep_count=" << keepCount
                                 << " comparable_pool=" << (comparablePool ? 1 : 0)
                                 << " after=" << selectionIndices.size();
                            TraceWriteLine(line.str());
                        }
                    }
                }

                int effectiveSingleRingSelectionMode = singleRingSelectionMode;
                if (singleRingSelectionMode == 6 && !selectionIndices.empty()) {
                    double minCost = std::numeric_limits<double>::infinity();
                    double maxCost = 0.0;
                    for (const std::size_t index : selectionIndices) {
                        const double cost = windowCandidates[index].cost;
                        minCost = std::min(minCost, cost);
                        maxCost = std::max(maxCost, cost);
                    }

                    const double spread = (maxCost - minCost) / (std::abs(minCost) + 1e-12);
                    effectiveSingleRingSelectionMode =
                        (spread <= singleRingAdaptiveSpreadThreshold) ? 4 : 3;

                    if (TraceEnabled() && collapseStep < singleRingTraceSteps) {
                        std::ostringstream line;
                        line << "TRACE_SINGLE_ADAPTIVE_MODE"
                             << " step=" << collapseStep
                             << " spread=" << spread
                             << " threshold=" << singleRingAdaptiveSpreadThreshold
                             << " effective_mode=" << effectiveSingleRingSelectionMode;
                        TraceWriteLine(line.str());
                    }
                }
                if (singleRingSelectionMode == 7) {
                    effectiveSingleRingSelectionMode =
                        (initialCount <= singleRingAdaptiveSizeThreshold) ? 1 : 3;

                    if (TraceEnabled() && collapseStep < singleRingTraceSteps) {
                        std::ostringstream line;
                        line << "TRACE_SINGLE_ADAPTIVE_SIZE_MODE"
                             << " step=" << collapseStep
                             << " initial_vertices=" << initialCount
                             << " threshold=" << singleRingAdaptiveSizeThreshold
                             << " effective_mode=" << effectiveSingleRingSelectionMode;
                        TraceWriteLine(line.str());
                    }
                }

                std::size_t bestSelectionPos = 0U;
                for (std::size_t pos = 1U; pos < selectionIndices.size(); ++pos) {
                    const std::size_t candidateIndex = selectionIndices[pos];
                    const std::size_t bestIndex = selectionIndices[bestSelectionPos];

                    const SingleRingPairCandidate& candidate = windowCandidates[candidateIndex];
                    const SingleRingPairCandidate& best = windowCandidates[bestIndex];

                    if (
                        effectiveSingleRingSelectionMode == 3
                        || effectiveSingleRingSelectionMode == 4
                        || effectiveSingleRingSelectionMode == 5
                    ) {
                        const auto bucketAge = [&](const SingleRingPairCandidate& windowCandidate) {
                            const std::size_t touchedAt = bucketLastTouch[windowCandidate.bucketIndex];
                            if (touchedAt == std::numeric_limits<std::size_t>::max() || touchedAt > collapseStep) {
                                return std::numeric_limits<std::size_t>::max();
                            }
                            return collapseStep - touchedAt;
                        };

                        const std::size_t candidateBucketAge = bucketAge(candidate);
                        const std::size_t bestBucketAge = bucketAge(best);

                        if (effectiveSingleRingSelectionMode == 3) {
                            if (candidateBucketAge > bestBucketAge) {
                                bestSelectionPos = pos;
                                continue;
                            }
                            if (candidateBucketAge < bestBucketAge) {
                                continue;
                            }

                            if (LessWithTolerance(candidate.cost, best.cost)) {
                                bestSelectionPos = pos;
                                continue;
                            }
                            if (LessWithTolerance(best.cost, candidate.cost)) {
                                continue;
                            }
                        } else if (effectiveSingleRingSelectionMode == 4) {
                            const bool comparable = CostsComparableForSingleRingSelection(candidate.cost, best.cost);
                            if (comparable) {
                                if (candidateBucketAge > bestBucketAge) {
                                    bestSelectionPos = pos;
                                    continue;
                                }
                                if (candidateBucketAge < bestBucketAge) {
                                    continue;
                                }
                            }

                            if (LessWithTolerance(candidate.cost, best.cost)) {
                                bestSelectionPos = pos;
                                continue;
                            }
                            if (LessWithTolerance(best.cost, candidate.cost)) {
                                continue;
                            }
                        } else {
                            const bool comparable = CostsComparableForSingleRingSelection(candidate.cost, best.cost);
                            if (comparable) {
                                if (LessWithTolerance(candidate.cornerSharpness, best.cornerSharpness)) {
                                    bestSelectionPos = pos;
                                    continue;
                                }
                                if (LessWithTolerance(best.cornerSharpness, candidate.cornerSharpness)) {
                                    continue;
                                }

                                if (candidateBucketAge > bestBucketAge) {
                                    bestSelectionPos = pos;
                                    continue;
                                }
                                if (candidateBucketAge < bestBucketAge) {
                                    continue;
                                }
                            }

                            if (LessWithTolerance(candidate.cost, best.cost)) {
                                bestSelectionPos = pos;
                                continue;
                            }
                            if (LessWithTolerance(best.cost, candidate.cost)) {
                                continue;
                            }
                        }

                        if (candidate.recentAge > best.recentAge) {
                            bestSelectionPos = pos;
                            continue;
                        }
                        if (candidate.recentAge < best.recentAge) {
                            continue;
                        }
                    } else if (singleRingSelectionMode == 2) {
                        // This mode prioritizes spatial dispersion directly by preferring candidates far from recent collapses
                        if (candidate.recentAge > best.recentAge) {
                            bestSelectionPos = pos;
                            continue;
                        }
                        if (candidate.recentAge < best.recentAge) {
                            continue;
                        }

                        if (LessWithTolerance(candidate.cost, best.cost)) {
                            bestSelectionPos = pos;
                            continue;
                        }
                        if (LessWithTolerance(best.cost, candidate.cost)) {
                            continue;
                        }
                    } else {
                        if (LessWithTolerance(candidate.cost, best.cost)) {
                            bestSelectionPos = pos;
                            continue;
                        }

                        if (LessWithTolerance(best.cost, candidate.cost)) {
                            continue;
                        }

                        if (CostsComparableForSingleRingSelection(candidate.cost, best.cost)) {
                            if (candidate.recentAge > best.recentAge) {
                                bestSelectionPos = pos;
                                continue;
                            }
                            if (candidate.recentAge < best.recentAge) {
                                continue;
                            }
                        }
                    }

                    if (candidate.indexB < best.indexB) {
                        bestSelectionPos = pos;
                    }
                }

                const std::size_t selectedWindowIndex = selectionIndices[bestSelectionPos];

                top = windowCandidates[selectedWindowIndex];
                foundCandidate = true;

                for (std::size_t i = 0U; i < windowCandidates.size(); ++i) {
                    if (i == selectedWindowIndex) {
                        continue;
                    }
                    heap.push(windowCandidates[i]);
                }
            }
        }

        if (!foundCandidate) {
            break;
        }

        if (TraceEnabled() && collapseStep < singleRingTraceSteps) {
            std::ostringstream line;
            line << "TRACE_SINGLE_PICK"
                 << " step=" << collapseStep
                 << " alive=" << aliveCount
                 << " b=" << top.indexB
                 << " bucket=" << top.bucketIndex
                 << " recent_age=" << top.recentAge
                 << " cost=" << top.cost
                 << " disp=" << top.displacement
                 << " E=(" << FormatPointForTrace(top.replacement) << ')';
            TraceWriteLine(line.str());
        }

        const std::size_t indexB = top.indexB;
        const std::size_t indexC = next[indexB];
        const std::size_t indexA = prev[indexB];
        const std::size_t indexD = next[indexC];

        // This optional guard blocks collapses that would introduce a local self-intersection in single-ring mode
        if (singleRingTopologyGuardEnabled) {
            if (!IsSingleRingPairCollapseTopologicallySafe(ring, prev, next, alive, indexB, top.replacement)) {
                ++topologyGuardRejectCount;
                if (TraceEnabled() && collapseStep < singleRingTraceSteps) {
                    std::ostringstream line;
                    line << "TRACE_SINGLE_REJECT"
                         << " step=" << collapseStep
                         << " reason=topology_guard"
                         << " b=" << top.indexB
                         << " bucket=" << top.bucketIndex;
                    TraceWriteLine(line.str());
                }
                continue;
            }
        }

        // This applies pair-collapse by replacing B with E and removing C from the linked ring
        ring.vertices[indexB] = top.replacement;
        alive[indexC] = false;
        totalDisplacement += top.displacement;
        --aliveCount;

        next[indexB] = indexD;
        prev[indexD] = indexB;

        ++generation[indexA];
        ++generation[indexB];
        ++generation[indexD];

        // This records the local neighborhood touched by the accepted collapse so dispersion policy can spread early edits
        lastTouchStep[indexA] = collapseStep;
        lastTouchStep[indexB] = collapseStep;
        lastTouchStep[indexC] = collapseStep;
        lastTouchStep[indexD] = collapseStep;
        bucketLastTouch[top.bucketIndex] = collapseStep;
        ++collapseStep;

        SingleRingPairCandidate c1;
        if (buildPairCandidate(indexA, c1)) {
            heap.push(c1);
        }

        SingleRingPairCandidate c2;
        if (buildPairCandidate(indexB, c2)) {
            heap.push(c2);
        }

        const std::size_t indexAA = prev[indexA];
        SingleRingPairCandidate c3;
        if (buildPairCandidate(indexAA, c3)) {
            heap.push(c3);
        }

        SingleRingPairCandidate c4;
        if (buildPairCandidate(indexD, c4)) {
            heap.push(c4);
        }
    }

    // This keeps ring start anchored at vertex 0 for deterministic output ordering
    const std::size_t start = 0U;

    if (!alive[start] || aliveCount < 3U) {
        return;
    }

    std::vector<atpps::Point> simplifiedVertices;
    simplifiedVertices.reserve(aliveCount);

    std::size_t cursor = start;
    do {
        simplifiedVertices.push_back(ring.vertices[cursor]);
        cursor = next[cursor];
    } while (cursor != start && simplifiedVertices.size() <= aliveCount + 1U);

    if (simplifiedVertices.size() >= 3U) {
        ring.vertices = std::move(simplifiedVertices);
    }

    if (TraceEnabled()) {
        std::ostringstream line;
        line << "TRACE_SINGLE_END"
             << " final_vertices=" << ring.vertices.size()
             << " displacement=" << totalDisplacement
             << " collapse_steps=" << collapseStep
             << " topo_rejects=" << topologyGuardRejectCount;
        TraceWriteLine(line.str());
        TraceWriteLine("TRACE_SINGLE_OUTPUT_" + FormatRingForTrace(ring));
    }
}

// This picks uniformly spaced vertices in index order for ultra-large rings where iterative collapse is too slow
void SimplifySingleRingByUniformSampling(atpps::Ring& ring, const std::size_t targetVertices) {
    const std::size_t count = ring.vertices.size();
    if (count <= targetVertices || targetVertices < 3U) {
        return;
    }

    std::vector<atpps::Point> sampled;
    sampled.reserve(targetVertices);

    // This preserves vertex 0 as the anchor to keep ring start deterministic
    sampled.push_back(ring.vertices[0]);

    std::size_t lastIndex = std::numeric_limits<std::size_t>::max();
    for (std::size_t i = 1U; i < targetVertices; ++i) {
        std::size_t index = (i * count) / targetVertices;
        if (index >= count) {
            index = count - 1U;
        }

        if (index == 0U) {
            index = 1U;
        }

        if (index == lastIndex && index + 1U < count) {
            ++index;
        }

        sampled.push_back(ring.vertices[index]);
        lastIndex = index;
    }

    if (sampled.size() >= 3U) {
        ring.vertices = std::move(sampled);
    }
}

// This applies one vertex deletion in-place after legality checks already passed
void ApplyRemovalInPlace(atpps::Polygon& polygon, const RemovalCandidate& candidate) {
    auto& ring = polygon.rings[candidate.ringIndex];
    ring.vertices.erase(ring.vertices.begin() + static_cast<std::ptrdiff_t>(candidate.vertexIndex));
}

// This scans all legal removals and returns the best one that preserves topology
bool FindBestLegalRemoval(const atpps::Polygon& polygon, RemovalCandidate& bestCandidate) {
    bool foundAny = false;

    for (std::size_t ringIndex = 0U; ringIndex < polygon.rings.size(); ++ringIndex) {
        const atpps::Ring& ring = polygon.rings[ringIndex];

        // This keeps rings valid by never dropping below a closed triangle
        if (ring.vertices.size() <= 3U) {
            continue;
        }

        for (std::size_t vertexIndex = 0U; vertexIndex < ring.vertices.size(); ++vertexIndex) {
            // This pins ring index 0 as a stable anchor for deterministic output ordering
            if (vertexIndex == 0U) {
                continue;
            }

            const RemovalCandidate candidate = BuildCandidate(ring, ringIndex, vertexIndex);
            if (!IsRemovalLocallyLegal(polygon, candidate)) {
                continue;
            }

            if (!foundAny || IsBetterCandidate(candidate, bestCandidate)) {
                bestCandidate = candidate;
                foundAny = true;
            }
        }
    }

    return foundAny;
}

// This computes arithmetic-mean center used as a stable scaling pivot for one ring
atpps::Point ComputeRingCentroid(const atpps::Ring& ring) {
    atpps::Point centroid;
    if (ring.vertices.empty()) {
        return centroid;
    }

    for (const atpps::Point& point : ring.vertices) {
        centroid.x += point.x;
        centroid.y += point.y;
    }

    const double invCount = 1.0 / static_cast<double>(ring.vertices.size());
    centroid.x *= invCount;
    centroid.y *= invCount;
    return centroid;
}

// This returns a ring scaled around a center so signed area changes while local shape stays similar
atpps::Ring ScaleRingAroundCenter(const atpps::Ring& ring, const atpps::Point& center, const double scale) {
    atpps::Ring scaled = ring;
    for (atpps::Point& point : scaled.vertices) {
        point.x = center.x + (point.x - center.x) * scale;
        point.y = center.y + (point.y - center.y) * scale;
    }
    return scaled;
}

// This checks whether two signed-area values are close enough under mixed relative and absolute tolerance
bool AreasCloseEnough(const double lhs, const double rhs) {
    const double delta = std::abs(lhs - rhs);
    const double scale = std::max(std::abs(lhs), std::max(std::abs(rhs), 1.0));
    return delta <= (kAreaRestoreAbsoluteTolerance + kAreaRestoreRelativeTolerance * scale);
}

// This appends one note fragment while preserving a single-line stderr note style
void AppendNoteMessage(std::string& noteMessage, const std::string& fragment) {
    if (fragment.empty()) {
        return;
    }

    if (!noteMessage.empty()) {
        noteMessage += " | ";
    }
    noteMessage += fragment;
}

// This attempts deterministic safe ring-area restoration by iterative damped scaling with topology checks
bool RestoreRingAreaSafely(
    atpps::Polygon& polygon,
    const std::size_t ringIndex,
    const double targetSignedArea,
    std::string& failureReason
) {
    atpps::Ring& ring = polygon.rings[ringIndex];

    // This keeps restoration logic safe for non-polygonal rings
    if (ring.vertices.size() < 3U) {
        failureReason = "ring has fewer than 3 vertices";
        return false;
    }

    // This uses one stable pivot so repeated scaling stays deterministic
    const atpps::Point center = ComputeRingCentroid(ring);

    // This uses deterministic damping steps from full correction toward conservative updates
    const double dampingFactors[] = {1.0, 0.5, 0.25, 0.125, 0.0625};

    for (std::size_t iteration = 0U; iteration < kMaxAreaRestoreIterations; ++iteration) {
        const double currentArea = atpps::ComputeSignedArea(ring);
        if (AreasCloseEnough(currentArea, targetSignedArea)) {
            failureReason.clear();
            return true;
        }

        // This avoids unstable correction when area sign is incompatible with target sign
        if (currentArea * targetSignedArea <= 0.0) {
            failureReason = "current and target ring signed areas have incompatible signs";
            return false;
        }

        // This avoids division by very small area values that would explode scaling factors
        if (std::abs(currentArea) <= kEpsilon) {
            failureReason = "current ring area is too close to zero for stable scaling";
            return false;
        }

        const double idealScale = std::sqrt(std::abs(targetSignedArea / currentArea));

        bool foundImprovement = false;
        atpps::Ring bestRing = ring;
        double bestError = std::abs(currentArea - targetSignedArea);

        for (const double damping : dampingFactors) {
            const double candidateScale = 1.0 + damping * (idealScale - 1.0);

            // This rejects collapsed or flipped scaling factors before geometry checks
            if (candidateScale <= kEpsilon) {
                continue;
            }

            atpps::Polygon candidatePolygon = polygon;
            candidatePolygon.rings[ringIndex] = ScaleRingAroundCenter(ring, center, candidateScale);

            std::string topologyError;
            if (!atpps::ValidatePolygonTopology(candidatePolygon, topologyError)) {
                continue;
            }

            const double candidateArea = atpps::ComputeSignedArea(candidatePolygon.rings[ringIndex]);
            const double candidateError = std::abs(candidateArea - targetSignedArea);

            // This accepts only strict improvement so iteration progress is monotonic
            if (candidateError + kEpsilon < bestError) {
                bestError = candidateError;
                bestRing = candidatePolygon.rings[ringIndex];
                foundImprovement = true;
            }
        }

        // This reports constrained restoration when no safe improving step exists
        if (!foundImprovement) {
            failureReason = "no topology-safe scaling improvement found";
            return false;
        }

        ring = bestRing;
    }

    const double finalArea = atpps::ComputeSignedArea(ring);
    if (AreasCloseEnough(finalArea, targetSignedArea)) {
        failureReason.clear();
        return true;
    }

    failureReason = "iteration cap reached before meeting area tolerance";
    return false;
}

// This restores ring signed areas toward input targets and records constrained-stop diagnostics
void RestorePolygonRingAreas(const atpps::Polygon& inputPolygon, atpps::Polygon& simplifiedPolygon, std::string& noteMessage) {
    std::map<int, double> targetAreaByRingId;
    for (const atpps::Ring& ring : inputPolygon.rings) {
        targetAreaByRingId[ring.ringId] = atpps::ComputeSignedArea(ring);
    }

    for (std::size_t ringIndex = 0U; ringIndex < simplifiedPolygon.rings.size(); ++ringIndex) {
        const int ringId = simplifiedPolygon.rings[ringIndex].ringId;
        const auto targetAreaIt = targetAreaByRingId.find(ringId);
        if (targetAreaIt == targetAreaByRingId.end()) {
            AppendNoteMessage(noteMessage, "ring " + std::to_string(ringId) + " has no input area target");
            continue;
        }

        const double targetArea = targetAreaIt->second;
        const double currentArea = atpps::ComputeSignedArea(simplifiedPolygon.rings[ringIndex]);
        if (AreasCloseEnough(currentArea, targetArea)) {
            continue;
        }

        std::string failureReason;
        const bool restored = RestoreRingAreaSafely(simplifiedPolygon, ringIndex, targetArea, failureReason);
        if (!restored) {
            std::ostringstream message;
            message << "ring " << ringId << " area restoration constrained: " << failureReason;
            AppendNoteMessage(noteMessage, message.str());
        }
    }
}

// This compares two polygons lexicographically using ring/vertex order and descending coordinates for tie resolution
bool IsPolygonLexicographicallyBetterForBeamTie(const atpps::Polygon& lhs, const atpps::Polygon& rhs) {
    if (lhs.rings.size() != rhs.rings.size()) {
        return lhs.rings.size() < rhs.rings.size();
    }

    for (std::size_t ringIndex = 0U; ringIndex < lhs.rings.size(); ++ringIndex) {
        const atpps::Ring& lhsRing = lhs.rings[ringIndex];
        const atpps::Ring& rhsRing = rhs.rings[ringIndex];

        if (lhsRing.ringId != rhsRing.ringId) {
            return lhsRing.ringId < rhsRing.ringId;
        }

        if (lhsRing.vertices.size() != rhsRing.vertices.size()) {
            return lhsRing.vertices.size() < rhsRing.vertices.size();
        }

        for (std::size_t vertexIndex = 0U; vertexIndex < lhsRing.vertices.size(); ++vertexIndex) {
            const atpps::Point& lhsPoint = lhsRing.vertices[vertexIndex];
            const atpps::Point& rhsPoint = rhsRing.vertices[vertexIndex];

            if (LessWithTolerance(rhsPoint.x, lhsPoint.x)) {
                return true;
            }
            if (LessWithTolerance(lhsPoint.x, rhsPoint.x)) {
                return false;
            }

            if (LessWithTolerance(rhsPoint.y, lhsPoint.y)) {
                return true;
            }
            if (LessWithTolerance(lhsPoint.y, rhsPoint.y)) {
                return false;
            }
        }
    }

    return false;
}

}  // namespace

namespace atpps {

SimplificationResult SimplifyPolygonToTarget(
    const Polygon& inputPolygon,
    const std::size_t targetVertices,
    std::string& noteMessage
) {
    // This starts from input and performs deterministic legal removals toward target
    SimplificationResult result;
    result.polygon = inputPolygon;
    result.requestedTarget = targetVertices;
    noteMessage.clear();

    std::size_t currentVertices = CountTotalVertices(result.polygon);
    if (targetVertices >= currentVertices) {
        result.finalVertexCount = currentVertices;
        result.reachedExactTarget = (result.finalVertexCount == targetVertices);
        if (!result.reachedExactTarget) {
            noteMessage = "target is greater than current vertex count so no removals were needed";
        }
        return result;
    }

    // This uses a faster deterministic path for single-ring cases to avoid repeated full rescans
    if (result.polygon.rings.size() == 1U) {
        const std::size_t clampedTarget = (targetVertices < 3U) ? 3U : targetVertices;
        if (result.polygon.rings[0].vertices.size() >= GetConfiguredHugeSingleRingVertexThreshold()) {
            SimplifySingleRingByUniformSampling(result.polygon.rings[0], clampedTarget);
        } else {
            double singleRingDisplacement = 0.0;
            SimplifySingleRingFast(result.polygon.rings[0], clampedTarget, singleRingDisplacement);
            result.totalArealDisplacement += singleRingDisplacement;
        }
        currentVertices = CountTotalVertices(result.polygon);
        result.finalVertexCount = currentVertices;
        result.reachedExactTarget = (result.finalVertexCount == targetVertices);
        if (!result.reachedExactTarget) {
            noteMessage = "target could not be reached under topology constraints; stopped at feasible vertex count";
        }
        return result;
    }

    // This runs bounded beam search for multi-ring inputs to reduce greedy local-optimum lock-in
    struct BeamState {
        atpps::Polygon polygon;
        std::size_t vertexCount = 0U;
        double cumulativeCost = 0.0;
        double totalDisplacement = 0.0;
    };

    std::vector<BeamState> beam;
    beam.push_back(BeamState{result.polygon, currentVertices, 0.0, 0.0});

    // This snapshots tuning knobs once per run so one process stays internally consistent
    const std::size_t configuredBeamWidth = GetConfiguredMultiRingBeamWidth();
    const std::size_t configuredBranchFactor = GetConfiguredMultiRingBranchFactor();
    const double configuredRingRoleCostMultiplier = GetConfiguredRingRoleComparableCostMultiplier();
    const double configuredMultiRingMovementWeight = GetConfiguredMultiRingMovementWeight();
    const int configuredMultiRingNormalizedScoreMode = GetConfiguredMultiRingNormalizedScoreMode();
    const int configuredApscSameSideMode = GetConfiguredApscSameSideMode();
    const int configuredDepthQuotaMode = GetConfiguredMultiRingDepthQuotaMode();

    if (TraceEnabled()) {
        std::ostringstream header;
        header << "TRACE_BEGIN target=" << targetVertices
               << " initial_vertices=" << currentVertices
               << " rings=" << result.polygon.rings.size()
               << " beam_width=" << configuredBeamWidth
             << " branch_factor=" << configuredBranchFactor
             << " ring_role_mult=" << configuredRingRoleCostMultiplier
             << " movement_weight=" << configuredMultiRingMovementWeight
                         << " normalized_score_mode=" << configuredMultiRingNormalizedScoreMode
               << " apsc_same_side_mode=" << configuredApscSameSideMode
               << " depth_quota_mode=" << configuredDepthQuotaMode;
        TraceWriteLine(header.str());

        for (const atpps::Ring& ring : result.polygon.rings) {
            TraceWriteLine("TRACE_INPUT_" + FormatRingForTrace(ring));
        }
    }

    const std::size_t maxSteps = currentVertices - targetVertices;
    for (std::size_t step = 0U; step < maxSteps; ++step) {
        std::vector<BeamState> nextBeam;

        if (TraceEnabled()) {
            std::ostringstream line;
            line << "TRACE_STEP step=" << step << " beam_states=" << beam.size();
            TraceWriteLine(line.str());
        }

        for (std::size_t stateIndex = 0U; stateIndex < beam.size(); ++stateIndex) {
            const BeamState& state = beam[stateIndex];

            if (TraceEnabled()) {
                std::ostringstream line;
                line << "TRACE_STATE step=" << step
                     << " state=" << stateIndex
                     << " vertices=" << state.vertexCount
                     << " cumulative_cost=" << state.cumulativeCost
                     << " displacement=" << state.totalDisplacement;
                TraceWriteLine(line.str());

                for (const atpps::Ring& ring : state.polygon.rings) {
                    std::ostringstream ringLine;
                    ringLine << "TRACE_STATE_RING step=" << step
                             << " state=" << stateIndex
                             << " " << FormatRingForTrace(ring);
                    TraceWriteLine(ringLine.str());
                }
            }

            if (state.vertexCount <= targetVertices) {
                if (TraceEnabled()) {
                    TraceWriteLine("TRACE_STATE_RETAIN reason=already_at_or_below_target");
                }
                nextBeam.push_back(state);
                continue;
            }

            std::vector<PairCollapseCandidate> candidates;
            CollectPairCollapseCandidates(state.polygon, configuredBranchFactor, candidates);
            if (candidates.empty()) {
                if (TraceEnabled()) {
                    TraceWriteLine("TRACE_STATE_RETAIN reason=no_legal_candidates");
                }
                nextBeam.push_back(state);
                continue;
            }

            if (TraceEnabled()) {
                std::ostringstream line;
                line << "TRACE_CANDIDATES step=" << step
                     << " state=" << stateIndex
                     << " count=" << candidates.size();
                TraceWriteLine(line.str());

                for (std::size_t candidateIndex = 0U; candidateIndex < candidates.size(); ++candidateIndex) {
                    std::ostringstream candidateLine;
                    candidateLine << "TRACE_CANDIDATE step=" << step
                                  << " state=" << stateIndex
                                  << " rank=" << candidateIndex
                                  << ' ' << FormatPairCandidateForTrace(candidates[candidateIndex]);
                    TraceWriteLine(candidateLine.str());
                }
            }

            for (const PairCollapseCandidate& candidate : candidates) {
                BeamState child;
                child.polygon = ApplyPairCollapse(
                    state.polygon,
                    candidate.ringIndex,
                    candidate.indexB,
                    candidate.indexC,
                    candidate.replacement
                );
                child.vertexCount = state.vertexCount - 1U;
                child.cumulativeCost = state.cumulativeCost + candidate.cost;
                child.totalDisplacement = state.totalDisplacement + candidate.displacement;

                if (TraceEnabled()) {
                    std::ostringstream line;
                    line << "TRACE_EXPAND step=" << step
                         << " state=" << stateIndex
                         << ' ' << FormatPairCandidateForTrace(candidate)
                         << " -> child_vertices=" << child.vertexCount
                         << " child_cost=" << child.cumulativeCost
                         << " child_disp=" << child.totalDisplacement;
                    TraceWriteLine(line.str());
                }

                nextBeam.push_back(std::move(child));
            }
        }

        if (nextBeam.empty()) {
            if (TraceEnabled()) {
                TraceWriteLine("TRACE_STOP reason=next_beam_empty");
            }
            break;
        }

        std::stable_sort(
            nextBeam.begin(),
            nextBeam.end(),
            [&](const BeamState& lhs, const BeamState& rhs) {
                const bool lhsHit = lhs.vertexCount == targetVertices;
                const bool rhsHit = rhs.vertexCount == targetVertices;
                if (lhsHit != rhsHit) {
                    return lhsHit;
                }

                if (lhs.vertexCount != rhs.vertexCount) {
                    return lhs.vertexCount < rhs.vertexCount;
                }

                if (LessWithTolerance(lhs.cumulativeCost, rhs.cumulativeCost)) {
                    return true;
                }
                if (LessWithTolerance(rhs.cumulativeCost, lhs.cumulativeCost)) {
                    return false;
                }

                if (LessWithTolerance(lhs.totalDisplacement, rhs.totalDisplacement)) {
                    return true;
                }
                if (LessWithTolerance(rhs.totalDisplacement, lhs.totalDisplacement)) {
                    return false;
                }

                return IsPolygonLexicographicallyBetterForBeamTie(lhs.polygon, rhs.polygon);
            }
        );

        if (nextBeam.size() > configuredBeamWidth) {
            nextBeam.resize(configuredBeamWidth);
        }

        if (TraceEnabled()) {
            for (std::size_t i = 0U; i < nextBeam.size(); ++i) {
                std::ostringstream line;
                line << "TRACE_NEXT_BEAM rank=" << i
                     << " vertices=" << nextBeam[i].vertexCount
                     << " cumulative_cost=" << nextBeam[i].cumulativeCost
                     << " displacement=" << nextBeam[i].totalDisplacement;
                TraceWriteLine(line.str());
            }
        }

        beam = std::move(nextBeam);
    }

    if (!beam.empty()) {
        const BeamState& bestState = beam.front();
        result.polygon = bestState.polygon;
        currentVertices = bestState.vertexCount;
        result.totalArealDisplacement += bestState.totalDisplacement;

        if (TraceEnabled()) {
            std::ostringstream line;
            line << "TRACE_FINAL vertices=" << currentVertices
                 << " cumulative_cost=" << bestState.cumulativeCost
                 << " displacement=" << bestState.totalDisplacement;
            TraceWriteLine(line.str());

            for (const atpps::Ring& ring : result.polygon.rings) {
                TraceWriteLine("TRACE_OUTPUT_" + FormatRingForTrace(ring));
            }
        }
    }

    result.finalVertexCount = currentVertices;
    result.reachedExactTarget = (result.finalVertexCount == targetVertices);
    if (!result.reachedExactTarget) {
        noteMessage = "target could not be reached under topology constraints; stopped at feasible vertex count";
    }
    return result;

    while (currentVertices > targetVertices) {
        RemovalCandidate bestCandidate;
        const bool foundRemoval = FindBestLegalRemoval(result.polygon, bestCandidate);
        if (!foundRemoval) {
            break;
        }

        ApplyRemovalInPlace(result.polygon, bestCandidate);
        --currentVertices;
    }

    // This keeps output on the direct collapse path so large datasets avoid expensive global restoration loops
    // This also preserves the assignment-style areal displacement signal instead of forcing near-zero drift

    result.finalVertexCount = currentVertices;
    result.reachedExactTarget = (result.finalVertexCount == targetVertices);

    // This reports constrained-stop cases through stderr while preserving stdout format
    if (!result.reachedExactTarget) {
        noteMessage = "target could not be reached under topology constraints; stopped at feasible vertex count";
    }

    return result;
}

}  // namespace atpps
