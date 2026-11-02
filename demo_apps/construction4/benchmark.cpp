// This file runs CPU rectangle comparisons between Walk on Tangent Spheres and
// Walk on Stars, and exposes the numerical validation checks.

#include "rectangle.h"
#include "validation.h"
#include <wosx/point_estimation/walk_on_stars.h>
#include <chrono>
#include <iomanip>
#include <string>
#include <limits>

int main(int argc, char** argv)
{
    using namespace wosx;
    using namespace construction4;
    try {
        int nWalks = 6000;
        int maxWalkLength = 100000;
        uint64_t seed = 20261030;
        double cap = .01;
        double radius = .4;
        double epsilonShellForAbsorbingBoundary = .01;
        double epsilonShellForReflectingBoundary = .01;
        std::string solverType = "both";
        std::string problemName = "polynomial";
        bool runTests = false;
        for (int i = 1; i < argc; ++i) {
            std::string key = argv[i];
            if (key == "--self-test") {
                runTests = true;
                continue;
            }
            if (key == "--help") {
                std::cout << "construction4_benchmark [--self-test] [--solver both|construction4|wost]\n"
                          << " [--case polynomial|harmonic|constant] [--walks N] [--seed N]\n"
                          << " [--cap d] [--radius r0] [--tol-dirichlet eD] [--tol-neumann eN] [--max-steps N]\n";
                return 0;
            }
            if (i + 1 == argc) {
                throw std::invalid_argument("Missing option value");
            }
            const std::string value = argv[++i];
            if (key == "--walks") {
                nWalks = std::stoi(value);
            } else if (key == "--seed") {
                seed = std::stoull(value);
            } else if (key == "--cap") {
                cap = std::stod(value);
            } else if (key == "--radius") {
                radius = std::stod(value);
            } else if (key == "--tol-dirichlet") {
                epsilonShellForAbsorbingBoundary = std::stod(value);
            } else if (key == "--tol-neumann") {
                epsilonShellForReflectingBoundary = std::stod(value);
            } else if (key == "--max-steps") {
                maxWalkLength = std::stoi(value);
            } else if (key == "--solver") {
                solverType = value;
            } else if (key == "--case") {
                problemName = value;
            } else {
                throw std::invalid_argument("Unknown option: " + key);
            }
        }
        if (runTests) {
            selfTest();
            return 0;
        }
        if (nWalks < 2 || maxWalkLength < 1 ||
            !(epsilonShellForAbsorbingBoundary > 0 && epsilonShellForReflectingBoundary > 0) ||
            !std::isfinite(epsilonShellForAbsorbingBoundary) || !std::isfinite(epsilonShellForReflectingBoundary)) {
            throw std::invalid_argument("Need >=2 walks, positive step limit and finite positive tolerances");
        }
        if (solverType != "both" && solverType != "construction4" && solverType != "wost") {
            throw std::invalid_argument("Unknown solver");
        }
        if (problemName != "polynomial" && problemName != "harmonic" && problemName != "constant") {
            throw std::invalid_argument("Unknown manufactured case");
        }
        Problem problem = problemName == "polynomial"
                              ? Problem::Polynomial
                              : (problemName == "harmonic" ? Problem::Harmonic : Problem::Constant);
        std::cout << std::setprecision(10);
        std::cout << "# Experimental Construction 4; cap/displacement errors omitted; no convergence claim\n"
                  << "# -Delta u=s; polynomial s=-0.6, outward bottom gN=-0.3\n"
                  << "# case=" << problemName << " seed=" << seed << " cap=" << cap << " radius=" << radius
                  << " tolD=" << epsilonShellForAbsorbingBoundary << " tolN=" << epsilonShellForReflectingBoundary
                  << " maxSteps=" << maxWalkLength << "\n";
        std::cout
            << "solver,walks,reference,mean,signed_error,mc_se,seconds,mean_steps,neumann_encounters,"
            << "distance_D,distance_N,project_D,project_N,ray_queries,boundary_samples,star_queries,"
            << "point_nearest_queries,center_nearest_queries,total_nearest_queries,radius_rejections,cap_adjustments,min_effective_cap,max_effective_cap,failures\n";
        for (const std::string algorithm : {std::string("construction4"), std::string("wost")}) {
            if (solverType != "both" && solverType != algorithm) {
                continue;
            }
            Rectangle rectangle;
            PDE<float, 2> pde;
            setupPDE(pde, problem);
            auto samplePt = samplePoint(seed);
            SampleStatistics<float, 2> statistics;
            WalkSettings walkSettings(float(epsilonShellForAbsorbingBoundary), float(epsilonShellForReflectingBoundary),
                                      maxWalkLength, false);
            TangentSphereDiagnostics diagnostics;
            uint64_t failures = 0;
            uint64_t encounters = 0;
            const auto start = std::chrono::steady_clock::now();
            if (algorithm == "construction4") {
                TangentSphereSettings tangentSettings;
                tangentSettings.capChordRadius = cap;
                tangentSettings.initialTangentRadius = radius;
                WalkOnTangentSpheres<float, 2> solver(rectangle.queries, rectangle.nearestQuery());
                try {
                    solver.solve(pde, walkSettings, tangentSettings, nWalks, samplePt, statistics, diagnostics);
                } catch (const std::exception& error) {
                    std::cerr << "Construction 4 failed: " << error.what()
                              << "; completed=" << diagnostics.walksCompleted
                              << " attempted=" << diagnostics.walksAttempted << " failures=" << diagnostics.failures
                              << ". No estimate reported for the failed batch.\n";
                    return 2;
                }
                encounters = diagnostics.tangentSteps;
            } else {
                WalkOnStars<float, 2> solver(
                    rectangle.queries,
                    [&](const WalkState<float, 2>& state) {
                        encounters += state.onReflectingBoundary;
                    },
                    [&](WalkCompletionCode, const WalkState<float, 2>&) {
                        ++failures;
                        return 0.f;
                    });
                solver.solve(pde, walkSettings, nWalks, samplePt, statistics);
                if (failures || statistics.getSolutionEstimateCount() != nWalks) {
                    std::cerr << "WoSt had " << failures
                              << " terminal failures; estimates=" << statistics.getSolutionEstimateCount()
                              << ". No comparison mean reported.\n";
                    return 2;
                }
            }
            const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            const double reference = exact(problem, samplePt.pt), mean = statistics.getEstimatedSolution();
            const double se =
                std::sqrt(double(statistics.getEstimatedSolutionVariance()) / statistics.getSolutionEstimateCount());
            const auto& queryCounts = rectangle.counts;
            std::cout << algorithm << ',' << statistics.getSolutionEstimateCount() << ',' << reference << ',' << mean
                      << ',' << mean - reference << ',' << se << ',' << seconds << ',' << statistics.getMeanWalkLength()
                      << ',' << encounters << ',' << queryCounts.absorbingDistance << ','
                      << queryCounts.reflectingDistance << ',' << queryCounts.absorbingProjection << ','
                      << queryCounts.reflectingProjection << ',' << queryCounts.ray << ',' << queryCounts.boundarySample
                      << ',' << queryCounts.starRadius << ',' << diagnostics.nearestPointQueries << ','
                      << diagnostics.centerCertificates << ','
                      << queryCounts.wholeBoundary.load() + queryCounts.absorbingDistance.load() +
                             queryCounts.reflectingDistance.load() + queryCounts.absorbingProjection.load() +
                             queryCounts.reflectingProjection.load()
                      << ',' << diagnostics.radiusRejections << ',' << diagnostics.capAdjustments << ',';
            if (std::isfinite(diagnostics.minEffectiveCap)) {
                std::cout << diagnostics.minEffectiveCap;
            } else {
                std::cout << "nan";
            }
            std::cout << ',' << diagnostics.maxEffectiveCap << ',' << failures + diagnostics.failures << '\n';
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
