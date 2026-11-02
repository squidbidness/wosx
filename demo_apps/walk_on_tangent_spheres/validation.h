// This file checks the tangent-disk formulas, sampling distributions and CPU
// solver interfaces against quadrature and manufactured solutions.

#pragma once
#include "rectangle.h"
#include <wosx/point_estimation/walk_on_stars.h>
#include <iostream>
#include <array>

namespace walk_on_tangent_spheres {
inline void require(bool condition, const char* message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

inline void close(double actual, double expected, double tolerance, const char* message)
{
    require(std::isfinite(actual) && std::abs(actual - expected) <= tolerance, message);
}
template <typename F>
double integrate(F f, double lo, double hi, int n = 200000)
{
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        sum += f(lo + (i + .5) * (hi - lo) / n);
    }
    return sum * (hi - lo) / n;
}

inline void validateDisk()
{
    using namespace wosx::tangent_disk;
    const double r = .4, d = .12, eta = .03, a = 2 * std::asin(d / (2 * r)), delta = .006;
    auto u = [](double x, double z) {
        return 1 + .3 * z + .2 * x * x + .1 * z * z;
    };
    const double uy = u(0, delta), uyp = u(0, delta + eta), g = -.3;
    auto value = [&](double theta) {
        return u(r * std::sin(theta), delta + 2 * r * std::pow(std::sin(theta / 2), 2));
    };
    auto flux = [&](double theta) {
        return r * finiteFlux(r, eta, theta);
    };
    const double mc = integrate(flux, -a, a), mr = integrate(flux, -pi, -a) + integrate(flux, a, pi);
    const double ir = integrate(
                          [&](double x) {
                              return value(x) * flux(x);
                          },
                          -pi, -a) +
                      integrate(
                          [&](double x) {
                              return value(x) * flux(x);
                          },
                          a, pi);
    const double ec = integrate(
        [&](double x) {
            return (value(x) - uy) * flux(x);
        },
        -a, a);
    const double e0 = uy - uyp - eta * g;
    // f=+0.6 here is explicitly the handoff convention; production PDE source s=-f.
    const double volume = finiteConstantVolume(r, eta, .6);
    close(mc + mr, 1, 2e-9, "Finite flux mass");
    close(mc * uy + ir + volume + ec, uyp, 2e-9, "Exact finite-source Green/cap identity");
    close(eta / mr * (mc * (g + e0 / eta) + ir / eta + volume / eta + ec / eta), uyp, 3e-8,
          "Rearranged identity with E0 and EC");
    const double mass = complementaryMass(r, d);
    close(integrate(
              [&](double x) {
                  return r * limitingFlux(r, x);
              },
              -pi, -a) +
              integrate(
                  [&](double x) {
                      return r * limitingFlux(r, x);
                  },
                  a, pi),
          mass, 2e-8, "Complementary limiting mass");
    close(limitingConstantVolume(r, .6), -.12, 1e-15, "Constant volume limit sign");
    const double smallEta = 1e-6;
    close(finiteConstantVolume(r, smallEta, .6) / smallEta, limitingConstantVolume(r, .6), 2e-7,
          "Exact scalar volume convergence");
    close(finiteFlux(r, smallEta, a) / smallEta, limitingFlux(r, a), 2e-4, "Complement flux pointwise limit");
    // Independent numerical area integration of the finite disk Green function.
    // Delta G=delta: log(r*|x-y'|/|r^2-(x-c)conj(y'-c)|)/(2pi).
    double volumeQuadrature = 0;
    const int nr = 700, nt = 1100;
    const double sourceZ = eta - r;
    for (int i = 0; i < nr; ++i) {
        const double rho = r * (i + .5) / nr;
        for (int j = 0; j < nt; ++j) {
            const double th = 2 * pi * (j + .5) / nt;
            const double x = rho * std::cos(th), z = rho * std::sin(th);
            const double numerator = r * std::hypot(x, z - sourceZ);
            const double denominator = std::hypot(r * r - z * sourceZ, x * sourceZ);
            volumeQuadrature += .6 * std::log(numerator / denominator) / (2 * pi) * rho * (r / nr) * (2 * pi / nt);
        }
    }
    close(volumeQuadrature, volume, 3e-7, "Independent finite volume quadrature");
    std::cout << "PASS: finite-source identity includes E0 and EC; disk mass and constant-source volume\n";
}

inline void validateSampler()
{
    using namespace wosx::tangent_disk;
    pcg32 rng(314159, 7);
    const double r = .4, d = .01, cot = cotHalfCap(r, d);
    const int n = 150000;
    double sum = 0;
    double sum2 = 0;
    double meanValue = 0;
    double m2 = 0;
    int positives = 0;
    std::array<int, 10> bins{};
    for (int i = 0; i < n; ++i) {
        const double first = rng.nextFloat(), second = rng.nextFloat();
        const double theta = sampleAngle(r, d, first, second);
        const double chord = 2 * r * std::abs(std::sin(theta / 2));
        require(chord >= d * (1 - 1e-12) && chord <= 2 * r, "Sampler chord cutoff");
        const double cdf = 1 - 1 / (std::tan(std::abs(theta) / 2) * cot);
        sum += cdf;
        sum2 += cdf * cdf;
        positives += theta > 0;
        ++bins[std::min(9, int(std::max(0., cdf) * 10))];
        const double x = r * std::sin(theta), z = .006 + 2 * r * std::pow(std::sin(theta / 2), 2);
        const double value = 1 + .3 * z + .2 * x * x + .1 * z * z;
        const double diff = value - meanValue;
        meanValue += diff / (i + 1);
        m2 += diff * (value - meanValue);
    }
    close(sum / n, .5, .004, "Sampler CDF mean");
    close(sum2 / n, 1. / 3, .004, "Sampler CDF second moment");
    close(double(positives) / n, .5, .005, "Sampler sign symmetry");
    for (int count : bins) {
        require(std::abs(count - n / 10) < 650, "Sampler CDF decile counts");
    }
    const double a = 2 * std::asin(d / (2 * r)), mass = complementaryMass(r, d);
    auto integrand = [&](double theta) {
        const double x = r * std::sin(theta), z = .006 + 2 * r * std::pow(std::sin(theta / 2), 2);
        return (1 + .3 * z + .2 * x * x + .1 * z * z) * r * limitingFlux(r, theta) / mass;
    };
    const double reference = integrate(integrand, -pi, -a) + integrate(integrand, a, pi);
    const double se = std::sqrt(m2 / (n - 1) / n);
    close(meanValue, reference, 6 * se + 1e-7, "Stage B analytic continuation mean");
    // Endpoints, near-diameter caps, and small d/r stability.
    for (double cap : {1e-10, .01, .799999}) {
        for (double uniform : {0., .5, .999999}) {
            const double theta = sampleAngle(r, cap, uniform, .25);
            require(std::isfinite(theta) && 2 * r * std::sin(theta / 2) >= cap * (1 - 1e-12), "Sampler extreme ratio");
        }
    }
    std::cout << "PASS: inverse-CDF distribution, chord cutoff, extreme ratios, Stage B analytic continuation\n";
    std::cout << "Stage B: continuation_mean=" << meanValue << ", quadrature=" << reference << ", SE=" << se
              << " (does not prove approximate recursion unbiased)\n";
}

inline void validateSolver()
{
    using namespace wosx;
    Rectangle rect;
    PDE<float, 2> pde;
    setupPDE(pde, Problem::Constant);
    WalkSettings w(.01f, .01f, 100000, false);
    TangentSphereSettings t;
    WalkOnTangentSpheres<float, 2> solver(rect.queries, rect.nearestQuery());
    auto pt = samplePoint(18);
    SampleStatistics<float, 2> stats;
    TangentSphereDiagnostics diag;
    solver.solve(pde, w, t, 1000, pt, stats, diag);
    close(stats.getEstimatedSolution(), 1, 0, "Constant solution exact preservation");
    require(diag.walksCompleted == 1000 && diag.nearestPointQueries == diag.interiorSteps + diag.tangentSteps + 1000 &&
                diag.tangentSteps > 0 && diag.failures == 0,
            "Recursive constant test completion");
    require(rect.counts.ray == 0 && rect.counts.starRadius == 0, "Walk on Tangent Spheres used ray/star queries");
    require(rect.counts.wholeBoundary.load() == diag.nearestPointQueries + diag.centerCertificates,
            "Whole-boundary query accounting mismatch");
    require(rect.counts.absorbingDistance == 0 && rect.counts.reflectingDistance == 0 &&
                rect.counts.absorbingProjection == 0 && rect.counts.reflectingProjection == 0,
            "Walk on Tangent Spheres used separate distance or projection queries");
    const auto nearJunction = rect.nearestBoundary(Vector2(.995f, .001f));
    require(nearJunction.isNeumann && nearJunction.normalIsUnique, "Nearest-boundary classification near junction");
    require(!rect.nearestBoundary(Vector2(.999f, .02f)).isNeumann, "Nearest Dirichlet classification");
    auto junctionPoint = samplePoint(91, Vector2(.995f, .001f));
    SampleStatistics<float, 2> junctionStats;
    TangentSphereDiagnostics junctionDiag;
    solver.solve(pde, w, t, 1, junctionPoint, junctionStats, junctionDiag);
    require(junctionDiag.tangentSteps > 0, "Dirichlet tolerance incorrectly overrode a nearer Neumann boundary");
    std::cout << "PASS: one whole-boundary query per visited point, no re-projections, nearest-type junction policy\n";
    auto expectReject = [&](auto invoke) {
        bool rejected = false;
        try {
            invoke();
        } catch (const std::exception&) {
            rejected = true;
        }
        require(rejected, "Unsupported input was silently accepted");
    };
    pde.isSourceConstant = false;
    expectReject([&] {
        solver.solve(pde, w, t, 1, pt, stats, diag);
    });
    pde.isSourceConstant = true;
    pde.absorptionCoeff = 1;
    expectReject([&] {
        solver.solve(pde, w, t, 1, pt, stats, diag);
    });
    pde.absorptionCoeff = 0;
    pde.areRobinConditionsPureNeumann = false;
    expectReject([&] {
        solver.solve(pde, w, t, 1, pt, stats, diag);
    });
    pde.areRobinConditionsPureNeumann = true;
    pt.estimationQuantity = EstimationQuantity::SolutionAndGradient;
    expectReject([&] {
        solver.solve(pde, w, t, 1, pt, stats, diag);
    });
    pt.estimationQuantity = EstimationQuantity::Solution;
    // A step-limit failure must not commit a zero or a survivor-only batch mean.
    auto limited = samplePoint(9);
    w.maxWalkLength = 1;
    const int previous = stats.getSolutionEstimateCount();
    expectReject([&] {
        solver.solve(pde, w, t, 2, limited, stats, diag);
    });
    require(stats.getSolutionEstimateCount() == previous && diag.failures > 0, "Failed batch changed statistics");
    w.maxWalkLength = 100000;
    TangentSphereSettings invalid = t;
    invalid.capChordRadius = 0;
    expectReject([&] {
        solver.solve(pde, w, invalid, 1, pt, stats, diag);
    });
    WalkOnTangentSpheres<float, 2> badGeometry(rect.queries, [](const Vector2&) {
        TangentBoundaryPoint2D hit;
        hit.signedClearance = -1;
        return hit;
    });
    expectReject([&] {
        badGeometry.solve(pde, w, t, 1, pt, stats, diag);
    });
    // Near a junction, radius halving must preserve containment and log cap fallback.
    auto nearCorner = samplePoint(77, Vector2(.98f, .005f));
    TangentSphereDiagnostics cornerDiag;
    SampleStatistics<float, 2> cornerStats;
    TangentSphereSettings largeCap = t;
    largeCap.capChordRadius = .1;
    solver.solve(pde, w, largeCap, 25, nearCorner, cornerStats, cornerDiag);
    require(cornerDiag.radiusRejections > 0 && cornerDiag.capAdjustments > 0,
            "Missing radius/cap fallback diagnostics");
    // Exercise the parallel batch overload and deterministic per-point RNG streams.
    std::vector<SamplePoint<float, 2>> points{samplePoint(23), samplePoint(24)};
    std::vector<SampleStatistics<float, 2>> parallelStats;
    std::vector<TangentSphereDiagnostics> parallelDiag;
    solver.solve(pde, w, t, std::vector<int>{20, 20}, points, parallelStats, parallelDiag);
    require(parallelStats.size() == 2 && parallelStats[0].getEstimatedSolution() == 1 &&
                parallelStats[1].getEstimatedSolution() == 1,
            "Parallel solution interface");
    std::cout << "PASS: recursion, failures, mode guards, cap fallback, parallel interface; no ray/star queries\n";
}

inline void validateManufacturedRecursion()
{
    using namespace wosx;
    // Fixed seeds protect against sign/source mistakes. The bias allowance is an
    // explicit practical regression threshold, not a convergence/unbiasedness test.
    constexpr int n = 6000;
    for (Problem problem : {Problem::Polynomial, Problem::Harmonic}) {
        Rectangle rect;
        PDE<float, 2> pde;
        setupPDE(pde, problem);
        WalkSettings w(.01f, .01f, 100000, false);
        TangentSphereSettings t;
        WalkOnTangentSpheres<float, 2> solver(rect.queries, rect.nearestQuery());
        auto pt = samplePoint(20261030);
        SampleStatistics<float, 2> statistics;
        TangentSphereDiagnostics diagnostics;
        solver.solve(pde, w, t, n, pt, statistics, diagnostics);
        const double reference = exact(problem, pt.pt);
        const double se = std::sqrt(double(statistics.getEstimatedSolutionVariance()) / n);
        close(statistics.getEstimatedSolution(), reference, 6 * se + .01,
              "Walk on Tangent Spheres manufactured recursion exceeds sampling+bias regression tolerance");
        require(diagnostics.walksCompleted == n &&
                    diagnostics.nearestPointQueries == diagnostics.interiorSteps + diagnostics.tangentSteps + n &&
                    diagnostics.failures == 0,
                "Manufactured recursion completion/certificate accounting");
        require(rect.counts.ray == 0 && rect.counts.starRadius == 0, "Manufactured recursion used ray/star queries");
        auto baselinePt = samplePoint(20261030);
        SampleStatistics<float, 2> baselineStatistics;
        uint64_t baselineFailures = 0;
        WalkOnStars<float, 2> baseline(rect.queries, {}, [&](WalkCompletionCode, const WalkState<float, 2>&) {
            ++baselineFailures;
            return 0.f;
        });
        baseline.solve(pde, w, n, baselinePt, baselineStatistics);
        require(baselineFailures == 0 && baselineStatistics.getSolutionEstimateCount() == n,
                "WoSt manufactured baseline did not complete");
        const double baselineSe = std::sqrt(double(baselineStatistics.getEstimatedSolutionVariance()) / n);
        close(baselineStatistics.getEstimatedSolution(), reference, 6 * baselineSe + .01,
              "WoSt shared PDE source/normal convention regression");
        std::cout << "PASS: " << (problem == Problem::Polynomial ? "polynomial" : "harmonic")
                  << " recursion; WoTS mean=" << statistics.getEstimatedSolution() << ", SE=" << se
                  << ", WoSt mean=" << baselineStatistics.getEstimatedSolution() << ", SE=" << baselineSe
                  << ", exact=" << reference << " (explicit bias allowance 0.01)\n";
    }
    // Identically seeded point streams must produce identical statistics and
    // query counts whether the per-point solves run serially or in parallel.
    Rectangle serialRect, parallelRect;
    PDE<float, 2> pde;
    setupPDE(pde, Problem::Polynomial);
    WalkSettings w(.01f, .01f, 100000, false);
    TangentSphereSettings t;
    std::vector<SamplePoint<float, 2>> a{samplePoint(202), samplePoint(303)};
    auto b = a;
    std::vector<SampleStatistics<float, 2>> serialStatistics, parallelStatistics;
    std::vector<TangentSphereDiagnostics> serialDiagnostics, parallelDiagnostics;
    WalkOnTangentSpheres<float, 2>(serialRect.queries, serialRect.nearestQuery())
        .solve(pde, w, t, std::vector<int>{200, 200}, a, serialStatistics, serialDiagnostics, true);
    WalkOnTangentSpheres<float, 2>(parallelRect.queries, parallelRect.nearestQuery())
        .solve(pde, w, t, std::vector<int>{200, 200}, b, parallelStatistics, parallelDiagnostics, false);
    for (size_t i = 0; i < a.size(); ++i) {
        close(serialStatistics[i].getEstimatedSolution(), parallelStatistics[i].getEstimatedSolution(), 0,
              "Serial/parallel seeded mean");
        close(serialStatistics[i].getEstimatedSolutionVariance(), parallelStatistics[i].getEstimatedSolutionVariance(),
              0, "Serial/parallel seeded variance");
        require(serialDiagnostics[i].tangentSteps == parallelDiagnostics[i].tangentSteps &&
                    serialDiagnostics[i].centerCertificates == parallelDiagnostics[i].centerCertificates,
                "Serial/parallel diagnostics mismatch");
    }
    const auto& c = serialRect.counts;
    const auto& d = parallelRect.counts;
    require(c.wholeBoundary.load() == d.wholeBoundary.load() &&
                c.absorbingDistance.load() == d.absorbingDistance.load() &&
                c.reflectingDistance.load() == d.reflectingDistance.load() &&
                c.absorbingProjection.load() == d.absorbingProjection.load() &&
                c.reflectingProjection.load() == d.reflectingProjection.load(),
            "Atomic query counts differ between serial and parallel runs");
    std::cout << "PASS: seeded serial/parallel means, variances, diagnostics and atomic query counts agree\n";
}

inline void selfTest()
{
    validateDisk();
    validateSampler();
    validateSolver();
    validateManufacturedRecursion();
}
} // walk_on_tangent_spheres
