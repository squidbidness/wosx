// This file implements the experimental CPU 2D walk-on-tangent-spheres solver.
// The practical recurrence omits cap and normal-displacement corrections; the
// implementation does not establish convergence or unbiasedness.

#pragma once

#include <wosx/point_estimation/common.h>
#include <wosx/point_estimation/tangent_disk.h>
#include <oneapi/tbb/parallel_for.h>
#include <stdexcept>
#include <string>
#include <cstdint>
#include <sstream>
#include <iomanip>

namespace wosx {

struct TangentSphereSettings {
    double capChordRadius = .01;
    double initialTangentRadius = .4;
    double radiusReduction = .5;
    double capFallbackRatio = 1.9; // d_eff=min(d,1.9r) only if d>=2r
    int maxRadiusTrials = 64;
    bool allowCornerWosFallback = false; // explicit mesh policy, logged below
    double containmentMargin = 0;        // conservative extra center clearance for float BVHs
};

struct TangentSphereDiagnostics {
    uint64_t walksAttempted = 0;
    uint64_t walksCompleted = 0;
    uint64_t failures = 0;
    uint64_t interiorSteps = 0;
    uint64_t tangentSteps = 0;
    uint64_t nearestPointQueries = 0;
    uint64_t cornerFallbacks = 0;
    uint64_t centerCertificates = 0;
    uint64_t radiusRejections = 0;
    uint64_t capAdjustments = 0;
    double minEffectiveCap = std::numeric_limits<double>::infinity();
    double maxEffectiveCap = 0;
    // Each nearestPointQueries call supplies distance, point, normal, type and sign.
    // Total whole-boundary nearest queries = nearestPointQueries + centerCertificates.
};

struct TangentBoundaryPoint2D {
    Vector2 point = Vector2::Zero();
    Vector2 normal = Vector2::Zero();                                  // physical outward unit normal
    double signedClearance = std::numeric_limits<double>::quiet_NaN(); // >=0 inside
    bool isNeumann = false;
    bool normalIsUnique = false;
};
// One search of the COMPLETE boundary. Must also provide a reliable inside/outside
// sign. Do not implement as independent Dirichlet/Neumann searches. A mesh adapter
// must validate its sign and normals; the initial adapter uses analytic geometry.
using TangentNearestBoundary2D = std::function<TangentBoundaryPoint2D(const Vector2&)>;

template <typename T, size_t DIM>
class WalkOnTangentSpheres;

template <typename T>
class WalkOnTangentSpheres<T, 2> {
public:
    // constructor
    WalkOnTangentSpheres(const GeometricQueries<2>& queries_, TangentNearestBoundary2D nearest);

    // Experimental approximate recurrence: drops the handoff's C* and D0.
    // PDE source follows WoSX: -Delta u=s, so s=-f in the handoff.
    // Rejects unsupported modes rather than inheriting WoSt's unrelated options.
    void solve(const PDE<T, 2>& pde,
               const WalkSettings& walkSettings,
               const TangentSphereSettings& tangentSettings,
               int nWalks,
               SamplePoint<T, 2>& samplePt,
               SampleStatistics<T, 2>& statistics,
               TangentSphereDiagnostics& diagnostics) const;

    void solve(const PDE<T, 2>& pde,
               const WalkSettings& walkSettings,
               const TangentSphereSettings& tangentSettings,
               const std::vector<int>& nWalks,
               std::vector<SamplePoint<T, 2>>& points,
               std::vector<SampleStatistics<T, 2>>& statistics,
               std::vector<TangentSphereDiagnostics>& diagnostics,
               bool runSingleThreaded = false,
               std::function<void(int, int)> reportProgress = {}) const;

private:
    // members
    const GeometricQueries<2>& queries;
    TangentNearestBoundary2D nearestBoundary;
    [[noreturn]] static void fail(const char* message);
    void validate(const PDE<T, 2>& pde,
                  const WalkSettings& walkSettings,
                  const TangentSphereSettings& tangentSettings) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation

template <typename T>
WalkOnTangentSpheres<T, 2>::WalkOnTangentSpheres(const GeometricQueries<2>& queries_, TangentNearestBoundary2D nearest)
    : queries(queries_), nearestBoundary(std::move(nearest))
{
    if (!nearestBoundary) {
        throw std::invalid_argument("A certified whole-boundary nearest-point callback is required");
    }
}

template <typename T>
void WalkOnTangentSpheres<T, 2>::solve(const PDE<T, 2>& pde,
                                       const WalkSettings& walkSettings,
                                       const TangentSphereSettings& tangentSettings,
                                       int nWalks,
                                       SamplePoint<T, 2>& samplePt,
                                       SampleStatistics<T, 2>& statistics,
                                       TangentSphereDiagnostics& diagnostics) const
{
    validate(pde, walkSettings, tangentSettings);
    if (samplePt.estimationQuantity == EstimationQuantity::Skip) {
        return;
    }
    if (samplePt.estimationQuantity != EstimationQuantity::Solution) {
        throw std::invalid_argument("Walk on Tangent Spheres initially supports solution estimates only");
    }
    if (nWalks <= 0) {
        throw std::invalid_argument("nWalks must be positive");
    }
    // Commit statistics only if the entire requested batch finishes successfully.
    auto pendingStatistics = statistics;
    for (int walkIndex = 0; walkIndex < nWalks; ++walkIndex) {
        ++diagnostics.walksAttempted;
        try {
            Vector2 currentPt = samplePt.pt;
            T total(0.0f);
            int steps = 0;
            for (;;) {
                // Query once at this position, including the final termination check.
                // Reuse the returned closest point and normal; never re-project.
                ++diagnostics.nearestPointQueries;
                const TangentBoundaryPoint2D nearest = nearestBoundary(currentPt);
                const double distance = nearest.signedClearance;
                if (!std::isfinite(distance) || distance < 0 || !nearest.point.allFinite()) {
                    fail("Invalid whole-boundary query or walk point outside the certified domain");
                }
                if (!nearest.isNeumann && distance <= walkSettings.epsilonShellForAbsorbingBoundary) {
                    total = total + pde.dirichlet(nearest.point, false);
                    break;
                }
                if (steps >= walkSettings.maxWalkLength) {
                    fail("Exceeded maximum walk length; batch has no reported estimate");
                }
                const bool nearNeumann =
                    nearest.isNeumann && distance <= walkSettings.epsilonShellForReflectingBoundary;
                const bool cornerFallback = nearNeumann && !nearest.normalIsUnique;
                if (cornerFallback && !tangentSettings.allowCornerWosFallback) {
                    fail("Nearest Neumann point has no unique outward normal");
                }
                if (cornerFallback) {
                    ++diagnostics.cornerFallbacks;
                }
                if (nearNeumann && !cornerFallback) {
                    const Vector2 reflectingBoundaryPt = nearest.point;
                    Vector2 normal = nearest.normal;
                    if (!nearest.normalIsUnique || !normal.allFinite() || normal.norm() < .5f) {
                        fail("Nearest Neumann point has no unique outward normal");
                    }
                    normal.normalize();
                    if ((currentPt - reflectingBoundaryPt).dot(normal) > 1e-6f ||
                        ((currentPt - reflectingBoundaryPt) + float(distance) * normal).norm() >
                            1e-5f + 1e-3f * float(distance)) {
                        fail("Neumann nearest point is not in a normal tubular neighborhood");
                    }
                    double tangentRadius = tangentSettings.initialTangentRadius;
                    bool accepted = false;
                    double lastRadius = tangentRadius;
                    double lastClearance = std::numeric_limits<double>::quiet_NaN();
                    Vector2 lastCenter = currentPt;
                    for (int trial = 0; trial < tangentSettings.maxRadiusTrials; ++trial) {
                        Vector2 c = currentPt - float(tangentRadius) * normal;
                        ++diagnostics.centerCertificates;
                        const double clearance = nearestBoundary(c).signedClearance;
                        lastRadius = tangentRadius;
                        lastClearance = clearance;
                        lastCenter = c;
                        // Strict conservative certificate; do not add a positive acceptance slack.
                        if (std::isfinite(clearance) &&
                            clearance >= tangentRadius + tangentSettings.containmentMargin) {
                            accepted = true;
                            break;
                        }
                        ++diagnostics.radiusRejections;
                        tangentRadius *= tangentSettings.radiusReduction;
                    }
                    if (!accepted || !(tangentRadius > 0) || float(tangentRadius) == 0) {
                        std::ostringstream detail;
                        detail << std::setprecision(17) << "Could not certify an admissible tangent disk"
                               << "; walk_index=" << walkIndex << " step=" << steps << " y=(" << currentPt.x() << ","
                               << currentPt.y() << ")"
                               << " nearest=(" << reflectingBoundaryPt.x() << "," << reflectingBoundaryPt.y() << ")"
                               << " normal=(" << normal.x() << "," << normal.y() << ")"
                               << " boundary_distance=" << distance << " r0=" << tangentSettings.initialTangentRadius
                               << " last_r=" << lastRadius << " last_center=(" << lastCenter.x() << ","
                               << lastCenter.y() << ")"
                               << " last_clearance=" << lastClearance
                               << " last_clearance_minus_r=" << lastClearance - lastRadius
                               << " required_margin=" << tangentSettings.containmentMargin
                               << " radius_trials=" << tangentSettings.maxRadiusTrials;
                        throw std::runtime_error(detail.str());
                    }
                    double capChordRadius = tangentSettings.capChordRadius;
                    if (capChordRadius >= 2 * tangentRadius) {
                        capChordRadius = tangentSettings.capFallbackRatio * tangentRadius;
                        ++diagnostics.capAdjustments;
                    }
                    diagnostics.minEffectiveCap = std::min(diagnostics.minEffectiveCap, capChordRadius);
                    diagnostics.maxEffectiveCap = std::max(diagnostics.maxEffectiveCap, capChordRadius);
                    const double complementaryMass = tangent_disk::complementaryMass(tangentRadius, capChordRadius);
                    // Known Neumann data evaluated at physical boundary, not at y.
                    T increment = pde.robin(reflectingBoundaryPt, normal, false) +
                                  pde.source(currentPt) * float(tangentRadius / 2);
                    total = total + increment * float(1 / complementaryMass);
                    const double angleSample = samplePt.rng.nextFloat(), signSample = samplePt.rng.nextFloat();
                    const double theta =
                        tangent_disk::sampleAngle(tangentRadius, capChordRadius, angleSample, signSample);
                    Vector2 tangent(-normal.y(), normal.x());
                    const double sinHalfAngle = std::sin(theta / 2);
                    Vector2 next = currentPt + float(tangentRadius * std::sin(theta)) * tangent -
                                   float(2 * tangentRadius * sinHalfAngle * sinHalfAngle) * normal;
                    if ((next - currentPt).squaredNorm() == 0) {
                        fail("Tangent step rounded to zero");
                    }
                    currentPt = next;
                    ++diagnostics.tangentSteps;
                } else {
                    const double sphereRadius = .99 * distance;
                    if (!(sphereRadius > 0)) {
                        fail("Zero interior sphere radius");
                    }
                    total = total + pde.source(currentPt) * float(sphereRadius * sphereRadius / 4);
                    const double theta = 2 * tangent_disk::pi * samplePt.rng.nextFloat();
                    Vector2 next =
                        currentPt + float(sphereRadius) * Vector2(float(std::cos(theta)), float(std::sin(theta)));
                    if ((next - currentPt).squaredNorm() == 0) {
                        fail("Interior step rounded to zero");
                    }
                    currentPt = next;
                    ++diagnostics.interiorSteps;
                }
                ++steps;
            }
            pendingStatistics.addSolutionEstimate(total);
            pendingStatistics.addWalkLength(steps);
            ++diagnostics.walksCompleted;
        } catch (...) {
            ++diagnostics.failures;
            throw;
        }
    }
    statistics = pendingStatistics;
}

template <typename T>
void WalkOnTangentSpheres<T, 2>::solve(const PDE<T, 2>& pde,
                                       const WalkSettings& walkSettings,
                                       const TangentSphereSettings& tangentSettings,
                                       const std::vector<int>& nWalks,
                                       std::vector<SamplePoint<T, 2>>& points,
                                       std::vector<SampleStatistics<T, 2>>& statistics,
                                       std::vector<TangentSphereDiagnostics>& diagnostics,
                                       bool runSingleThreaded,
                                       std::function<void(int, int)> reportProgress) const
{
    validate(pde, walkSettings, tangentSettings);
    if (points.size() != nWalks.size()) {
        throw std::invalid_argument("Point/walk count size mismatch");
    }
    auto pendingStatistics = statistics;
    pendingStatistics.resize(points.size());
    diagnostics.resize(points.size());
    auto run = [&](size_t pointIndex) {
        solve(pde, walkSettings, tangentSettings, nWalks[pointIndex], points[pointIndex], pendingStatistics[pointIndex],
              diagnostics[pointIndex]);
        if (reportProgress) {
            reportProgress(1, 0);
        }
    };
    if (runSingleThreaded) {
        for (size_t pointIndex = 0; pointIndex < points.size(); ++pointIndex) {
            run(pointIndex);
        }
    } else {
        tbb::parallel_for(size_t(0), points.size(), run);
    }
    statistics = std::move(pendingStatistics);
}

template <typename T>
[[noreturn]] void WalkOnTangentSpheres<T, 2>::fail(const char* message)
{
    throw std::runtime_error(message);
}

template <typename T>
void WalkOnTangentSpheres<T, 2>::validate(const PDE<T, 2>& pde,
                                          const WalkSettings& walkSettings,
                                          const TangentSphereSettings& tangentSettings) const
{
    if (!queries.domainIsWatertight || !queries.hasNonEmptyAbsorbingBoundary ||
        !queries.hasNonEmptyReflectingBoundary) {
        throw std::invalid_argument("Walk on Tangent Spheres requires a watertight mixed interior domain");
    }
    if (walkSettings.solveDoubleSided || pde.absorptionCoeff != 0 || !pde.areRobinConditionsPureNeumann ||
        !pde.isSourceConstant || !pde.source || !pde.dirichlet || !pde.robin) {
        throw std::invalid_argument(
            "Walk on Tangent Spheres requires pure Neumann reflection and constant source, without absorption");
    }
    if (walkSettings.ignoreAbsorbingBoundaryContribution || walkSettings.ignoreReflectingBoundaryContribution ||
        walkSettings.ignoreSourceContribution || walkSettings.russianRouletteThreshold != 0 ||
        walkSettings.splittingThreshold != std::numeric_limits<float>::max() ||
        walkSettings.stepsBeforeApplyingTikhonov != 0) {
        throw std::invalid_argument("Contribution suppression, splitting, roulette and Tikhonov are unsupported");
    }
    if (!(std::isfinite(walkSettings.epsilonShellForAbsorbingBoundary) &&
          walkSettings.epsilonShellForAbsorbingBoundary > 0 &&
          std::isfinite(walkSettings.epsilonShellForReflectingBoundary) &&
          walkSettings.epsilonShellForReflectingBoundary > 0 && walkSettings.maxWalkLength > 0 &&
          std::isfinite(tangentSettings.capChordRadius) && tangentSettings.capChordRadius > 0 &&
          std::isfinite(tangentSettings.initialTangentRadius) && tangentSettings.initialTangentRadius > 0 &&
          tangentSettings.radiusReduction > 0 && tangentSettings.radiusReduction < 1 &&
          tangentSettings.capFallbackRatio > 0 && tangentSettings.capFallbackRatio < 2 &&
          tangentSettings.maxRadiusTrials > 0 && std::isfinite(tangentSettings.containmentMargin) &&
          tangentSettings.containmentMargin >= 0)) {
        throw std::invalid_argument("Invalid Walk on Tangent Spheres settings");
    }
}

} // wosx
