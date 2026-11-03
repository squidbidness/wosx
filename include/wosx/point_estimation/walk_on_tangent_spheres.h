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
    double containmentMargin = 0; // optional extra clearance; no mesh-specific margin
};

struct TangentSphereDiagnostics {
    uint64_t walksAttempted = 0;
    uint64_t walksCompleted = 0;
    uint64_t failures = 0;
    uint64_t interiorSteps = 0;
    uint64_t tangentSteps = 0;
    uint64_t nearestPointQueries = 0;
    uint64_t initialMembershipChecks = 0;
    uint64_t centerMembershipChecks = 0;
    uint64_t centerCertificates = 0;
    uint64_t radiusRejections = 0;
    uint64_t capAdjustments = 0;
    double minEffectiveCap = std::numeric_limits<double>::infinity();
    double maxEffectiveCap = 0;
    // Logical projection/center checks are separate from underlying callback counts.
    // Default typed projection searches both boundary structures. Membership checks
    // may issue ray searches through insideDomain; benchmark counters include those.
};

template <typename T, size_t DIM>
class WalkOnTangentSpheres;

template <typename T>
class WalkOnTangentSpheres<T, 2> {
public:
    // constructor
    WalkOnTangentSpheres(const GeometricQueries<2>& queries_);

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
    [[noreturn]] static void fail(const char* message);
    void validate(const PDE<T, 2>& pde,
                  const WalkSettings& walkSettings,
                  const TangentSphereSettings& tangentSettings) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation

template <typename T>
WalkOnTangentSpheres<T, 2>::WalkOnTangentSpheres(const GeometricQueries<2>& queries_) : queries(queries_)
{
    // use the same geometry callbacks and representation as WoSt
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
    ++diagnostics.initialMembershipChecks;
    if (!queries.insideDomain(samplePt.pt)) {
        fail("Initial point is outside the domain");
    }
    for (int walkIndex = 0; walkIndex < nWalks; ++walkIndex) {
        ++diagnostics.walksAttempted;
        try {
            Vector2 currentPt = samplePt.pt;
            T total(0.0f);
            int steps = 0;
            for (;;) {
                // Keep the partition label selected by the existing projection
                // queries, rather than classifying the projected point again.
                ++diagnostics.nearestPointQueries;
                Vector2 boundaryPt = currentPt;
                Vector2 queryNormal;
                float boundaryDistance;
                BoundaryType boundaryType;
                if (!queries.projectToBoundaryWithType(boundaryPt, queryNormal, boundaryDistance, boundaryType,
                                                       false)) {
                    fail("Boundary projection failed");
                }
                const double distance = boundaryDistance;
                if (!std::isfinite(distance) || distance < 0 || !boundaryPt.allFinite() ||
                    boundaryType == BoundaryType::None || queries.outsideBoundingDomain(currentPt)) {
                    fail("Invalid boundary projection or walk escaped the bounding domain");
                }
                if (boundaryType == BoundaryType::Absorbing &&
                    distance <= walkSettings.epsilonShellForAbsorbingBoundary) {
                    total = total + pde.dirichlet(boundaryPt, false);
                    break;
                }
                if (steps >= walkSettings.maxWalkLength) {
                    fail("Exceeded maximum walk length; batch has no reported estimate");
                }
                if (boundaryType == BoundaryType::Reflecting &&
                    distance <= walkSettings.epsilonShellForReflectingBoundary) {
                    const Vector2 reflectingBoundaryPt = boundaryPt;
                    // This approach direction is also parallel to y0-y'. It is
                    // independent of the primitive normal, including at corners.
                    const Vector2 displacement = reflectingBoundaryPt - currentPt;
                    const float displacementLength = displacement.norm();
                    if (!(displacementLength > 0) || !std::isfinite(displacementLength)) {
                        fail("Approach direction is undefined at an exact boundary point");
                    }
                    const Vector2 normal = displacement / displacementLength;
                    double tangentRadius = tangentSettings.initialTangentRadius;
                    bool accepted = false;
                    double lastRadius = tangentRadius;
                    double lastClearance = std::numeric_limits<double>::quiet_NaN();
                    Vector2 lastCenter = currentPt;
                    for (int trial = 0; trial < tangentSettings.maxRadiusTrials; ++trial) {
                        Vector2 c = currentPt - float(tangentRadius) * normal;
                        ++diagnostics.centerCertificates;
                        const double clearance = queries.computeDistToBoundary(c, false);
                        ++diagnostics.centerMembershipChecks;
                        const bool centerInside = queries.insideDomain(c);
                        lastRadius = tangentRadius;
                        lastClearance = clearance;
                        lastCenter = c;
                        // Strict conservative certificate; do not add a positive acceptance slack.
                        if (centerInside && std::isfinite(clearance) &&
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
                    // Neumann data are evaluated with the approach direction.
                    // At corners this is the chosen convention, not necessarily a
                    // face-normal trace; its discrepancy remains an approximation.
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
