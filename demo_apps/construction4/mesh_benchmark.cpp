// This file compares CPU solvers on the supplied 2D mesh problems using the
// existing WoSX geometry and PDE callbacks.

#include "mesh_adapter.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>

namespace {
using namespace wosx;
using namespace construction4_mesh;
struct BenchmarkResult {
    double mean = 0;
    double standardError = 0;
    double seconds = 0;
    uint64_t walks = 0;
    uint64_t steps = 0;
    uint64_t queries = 0;
    uint64_t centers = 0;
    uint64_t cap = 0;
    uint64_t corners = 0;
    uint64_t rejections = 0;
    uint64_t failures = 0;
};

struct BenchmarkOptions {
    std::string root = ".";
    std::string config = "demo_apps/basic_2d/model_problems/engine/wost.json";
    std::string problem = "manufactured";
    int points = 16;
    int walks = 6000;
    int referenceWalks = 65536;
    int maxSteps = 100000;
    uint64_t seed = 20261030;
    double cap = .01;
    double radius = .4;
    double tolD = .01;
    double tolN = .01;
    double referenceTolerance = 1e-4;
    bool validateOnly = false;
};
inline double evaluateExactSolution(const Vector2& pt)
{
    return 1 + .3 * pt.y() + .2 * double(pt.x()) * pt.x() + .1 * double(pt.y()) * pt.y();
}
void setupPDE(PDE<float, 2>& pde, MeshModel& model, bool manufactured)
{
    pde.isSourceConstant = true;
    pde.absorptionCoeff = 0;
    pde.areRobinConditionsPureNeumann = true;
    if (manufactured) {
        pde.source = [](const Vector2&) {
            return -.6f;
        };
        pde.dirichlet = [](const Vector2& pt, bool) {
            return float(evaluateExactSolution(pt));
        };
        pde.robin = [](const Vector2& pt, const Vector2& n, bool) {
            return .4f * pt.x() * n.x() + (.3f + .2f * pt.y()) * n.y();
        };
    } else {
        // Match disabled source and reflecting contributions in the supplied config.
        pde.source = [](const Vector2&) {
            return 0.f;
        };
        pde.dirichlet = model.getPDE().dirichlet;
        pde.robin = [](const Vector2&, const Vector2&, bool) {
            return 0.f;
        };
    }
    pde.robinCoeff = [](const Vector2&, const Vector2&, bool) {
        return 0.f;
    };
    pde.hasReflectingBoundaryConditions = model.getPDE().hasReflectingBoundaryConditions;
}
BenchmarkResult runSolver(MeshModel& model,
                          MeshAdapter& mesh,
                          const PDE<float, 2>& pde,
                          const Vector2& pt,
                          const BenchmarkOptions& options,
                          bool useTangentSpheres,
                          bool computeReference,
                          uint64_t seed)
{
    const int walks = computeReference ? options.referenceWalks : options.walks;
    const auto& original = model.getGeometricQueries();
    SamplePoint<float, 2> sample(pt, Vector2::Zero(), SampleType::InDomain, EstimationQuantity::Solution, 1,
                                 original.computeDistToAbsorbingBoundary(pt, false),
                                 original.computeDistToReflectingBoundary(pt, false));
    sample.rng = pcg32(seed, 17);
    SampleStatistics<float, 2> statistics;
    WalkSettings walkSettings(float(computeReference ? options.referenceTolerance : options.tolD),
                              float(computeReference ? options.referenceTolerance : options.tolN), options.maxSteps,
                              false);
    TangentSphereDiagnostics diagnostics;
    WoStCounts counts;
    auto queries = countedQueries(original, counts);
    BenchmarkResult result;
    result.walks = walks;
    mesh.queryCount = 0;
    const auto start = std::chrono::steady_clock::now();
    if (useTangentSpheres) {
        TangentSphereSettings tangentSettings;
        tangentSettings.capChordRadius = options.cap;
        tangentSettings.initialTangentRadius = options.radius;
        tangentSettings.allowCornerWosFallback = true;
        tangentSettings.containmentMargin = 8 * std::numeric_limits<float>::epsilon() * mesh.scale;
        WalkOnTangentSpheres<float, 2> solver(queries, mesh.callback());
        try {
            solver.solve(pde, walkSettings, tangentSettings, walks, sample, statistics, diagnostics);
        } catch (const std::exception& error) {
            throw std::runtime_error(
                std::string("Construction 4 failed: ") + error.what() + "; evaluation_point=(" +
                std::to_string(pt.x()) + "," + std::to_string(pt.y()) + ")" + "; seed=" + std::to_string(seed) +
                "; completed=" + std::to_string(diagnostics.walksCompleted) +
                " attempted=" + std::to_string(diagnostics.walksAttempted) +
                " corners=" + std::to_string(diagnostics.cornerFallbacks) + "; no failed-batch estimate reported");
        }
        result.queries = mesh.queryCount.load();
        result.steps = diagnostics.interiorSteps + diagnostics.tangentSteps;
        result.centers = diagnostics.centerCertificates;
        result.cap = diagnostics.capAdjustments;
        result.corners = diagnostics.cornerFallbacks;
        result.rejections = diagnostics.radiusRejections;
        if (result.queries != diagnostics.nearestPointQueries + diagnostics.centerCertificates ||
            diagnostics.nearestPointQueries != result.steps + result.walks) {
            throw std::runtime_error("Mesh nearest-query accounting mismatch");
        }
    } else {
        walkSettings.ignoreSourceContribution = options.problem == "engine";
        walkSettings.ignoreReflectingBoundaryContribution = options.problem == "engine";
        WalkOnStars<float, 2> solver(queries, {}, [&](WalkCompletionCode, const WalkState<float, 2>&) {
            ++result.failures;
            return 0.f;
        });
        solver.solve(pde, walkSettings, walks, sample, statistics);
        if (result.failures || statistics.getSolutionEstimateCount() != walks) {
            throw std::runtime_error("WoSt/reference had failed or missing walks; comparison aborted");
        }
        result.queries = counts.total();
        result.steps = uint64_t(std::llround(double(statistics.getMeanWalkLength()) * walks));
    }
    result.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    result.mean = statistics.getEstimatedSolution();
    result.standardError =
        std::sqrt(double(statistics.getEstimatedSolutionVariance()) / statistics.getSolutionEstimateCount());
    return result;
}
void writeResult(const char* kind,
                 int index,
                 const Vector2& pt,
                 const char* solver,
                 const BenchmarkResult& result,
                 double referenceValue,
                 double referenceStandardError)
{
    std::cout << kind << ',' << index << ',' << pt.x() << ',' << pt.y() << ',' << solver << ',' << referenceValue << ','
              << referenceStandardError << ',' << result.mean << ',' << result.mean - referenceValue << ','
              << result.standardError << ',' << std::hypot(result.standardError, referenceStandardError) << ','
              << result.seconds << ',' << double(result.steps) / result.walks << ','
              << double(result.queries) / result.walks << ',' << double(result.centers) / result.walks << ','
              << double(result.cap) / result.walks << ',' << double(result.corners) / result.walks << ','
              << double(result.rejections) / result.walks << ',' << result.failures << '\n';
}
struct BenchmarkSummary {
    double squaredErrorSum = 0;
    double varianceSum = 0;
    double referenceVarianceSum = 0;
    double maxAbsoluteError = 0;
    double seconds = 0;
    uint64_t walks = 0;
    uint64_t steps = 0;
    uint64_t queries = 0;
    uint64_t centers = 0;
    uint64_t cap = 0;
    uint64_t corners = 0;
    uint64_t rejections = 0;
    uint64_t failures = 0;
    void add(const BenchmarkResult& result, double referenceValue, double referenceStandardError)
    {
        const double error = result.mean - referenceValue;
        squaredErrorSum += error * error;
        varianceSum += result.standardError * result.standardError;
        referenceVarianceSum += referenceStandardError * referenceStandardError;
        maxAbsoluteError = std::max(maxAbsoluteError, std::abs(error));
        seconds += result.seconds;
        walks += result.walks;
        steps += result.steps;
        queries += result.queries;
        centers += result.centers;
        cap += result.cap;
        corners += result.corners;
        rejections += result.rejections;
        failures += result.failures;
    }
    void print(const char* solver, int n) const
    {
        std::cout << "SUMMARY," << solver << ",points=" << n << ",walks=" << walks
                  << ",rmse=" << std::sqrt(squaredErrorSum / n) << ",max_abs_error=" << maxAbsoluteError
                  << ",rms_mc_se=" << std::sqrt(varianceSum / n)
                  << ",rms_reference_se=" << std::sqrt(referenceVarianceSum / n) << ",seconds=" << seconds
                  << ",avg_steps=" << double(steps) / walks << ",avg_queries=" << double(queries) / walks
                  << ",avg_center_checks=" << double(centers) / walks << ",avg_cap_fallback=" << double(cap) / walks
                  << ",avg_corner_fallback=" << double(corners) / walks
                  << ",avg_radius_rejections=" << double(rejections) / walks << ",failures=" << failures << '\n';
    }
};
}
int main(int argc, char** argv)
{
    try {
        BenchmarkOptions options;
        for (int i = 1; i < argc; ++i) {
            std::string key = argv[i];
            if (key == "--validate-only") {
                options.validateOnly = true;
                continue;
            }
            if (key == "--help") {
                std::cout << "construction4_mesh_benchmark --case manufactured|engine [--root repo] [--config path]\n"
                          << " [--points 16] [--walks 6000] [--seed 20261030] [--cap .01] [--radius .4]\n"
                          << " [--tol-dirichlet .01] [--tol-neumann .01] [--max-steps 100000]\n"
                          << " [--reference-walks 65536] [--reference-tolerance .0001] [--validate-only]\n";
                return 0;
            }
            if (i + 1 == argc) {
                throw std::invalid_argument("Missing option value");
            }
            std::string value = argv[++i];
            if (key == "--case") {
                options.problem = value;
            } else if (key == "--root") {
                options.root = value;
            } else if (key == "--config") {
                options.config = value;
            } else if (key == "--points") {
                options.points = std::stoi(value);
            } else if (key == "--walks") {
                options.walks = std::stoi(value);
            } else if (key == "--reference-walks") {
                options.referenceWalks = std::stoi(value);
            } else if (key == "--seed") {
                options.seed = std::stoull(value);
            } else if (key == "--cap") {
                options.cap = std::stod(value);
            } else if (key == "--radius") {
                options.radius = std::stod(value);
            } else if (key == "--tol-dirichlet") {
                options.tolD = std::stod(value);
            } else if (key == "--tol-neumann") {
                options.tolN = std::stod(value);
            } else if (key == "--reference-tolerance") {
                options.referenceTolerance = std::stod(value);
            } else if (key == "--max-steps") {
                options.maxSteps = std::stoi(value);
            } else {
                throw std::invalid_argument("Unknown option: " + key);
            }
        }
        if (options.problem != "engine" && options.problem != "manufactured") {
            throw std::invalid_argument("Unknown problem");
        }
        if (options.points < 1 || options.walks < 2 || options.referenceWalks < 2 || options.maxSteps < 1) {
            throw std::invalid_argument("Invalid counts");
        }
        for (double x : {options.cap, options.radius, options.tolD, options.tolN, options.referenceTolerance}) {
            if (!std::isfinite(x) || x <= 0) {
                throw std::invalid_argument("Parameters must be finite and positive");
            }
        }
        const auto root = std::filesystem::absolute(options.root);
        std::ifstream input(root / options.config);
        if (!input) {
            throw std::runtime_error("Cannot read mesh configuration; run from repository root or set --root");
        }
        json config;
        input >> config;
        const auto& modelProblemConfig = config.at("modelProblem");
        if (getOptional<bool>(modelProblemConfig, "solveDoubleSided", false) ||
            getOptional<bool>(modelProblemConfig, "solveExterior", false) ||
            !getOptional<bool>(modelProblemConfig, "domainIsWatertight", false) ||
            getOptional<float>(modelProblemConfig, "robinCoeff", 0) != 0 ||
            getOptional<float>(modelProblemConfig, "absorptionCoeff", 0) != 0) {
            throw std::invalid_argument("Only watertight single-sided mixed meshes without Robin/screening supported");
        }
        if (options.problem == "engine" &&
            (!getOptional<bool>(config.at("solver"), "ignoreSourceContribution", false) ||
             !getOptional<bool>(config.at("solver"), "ignoreReflectingBoundaryContribution", false))) {
            throw std::invalid_argument("Supplied-data mode requires disabled source and reflecting contributions");
        }
        MeshModel model(modelProblemConfig, root.string());
        MeshAdapter mesh(model);
        const size_t checks = mesh.validateQueries();
        std::cout << std::setprecision(10) << "# mesh=" << modelProblemConfig.at("geometry")
                  << " nonzero_segments=" << mesh.edges.size() << " omitted_zero_edges=" << mesh.omittedZeroEdges
                  << " geometry_validation_points=" << checks << '\n';
        if (options.validateOnly) {
            std::cout << "PASS: closed topology, crossings, orientation, winding signs and brute-force distances\n";
            return 0;
        }
        const auto box = model.getBoundingBox();
        std::vector<Vector2> candidates;
        for (int j = 0; j < 64; ++j) {
            for (int i = 0; i < 64; ++i) {
                Vector2 pt = box.first + Vector2((i + .5f) / 64 * (box.second.x() - box.first.x()),
                                                 (j + .5f) / 64 * (box.second.y() - box.first.y()));
                if (mesh.nearest(pt).signedClearance > .02 * mesh.scale) {
                    candidates.push_back(pt);
                }
            }
        }
        if (candidates.size() < size_t(options.points)) {
            throw std::runtime_error("Not enough interior evaluation candidates");
        }
        std::vector<Vector2> points;
        for (int i = 0; i < options.points; ++i) {
            points.push_back(candidates[(size_t(2 * i + 1) * candidates.size()) / (2 * options.points)]);
        }
        std::cout
            << "# case=" << options.problem << " points=" << options.points << " walks_per_point=" << options.walks
            << " seed=" << options.seed << " cap=" << options.cap << " radius=" << options.radius
            << " epsilonD=" << options.tolD << " epsilonN=" << options.tolN
            << " reference_walks=" << options.referenceWalks << " reference_tolerance=" << options.referenceTolerance
            << '\n'
            << "# C4 queries=whole-boundary position+center searches; WoSt queries=distance,projection,star,ray,sample,visibility.\n"
            << "# Single-threaded timings exclude loading, BVH builds, validation, reference and initial point metadata.\n"
            << "# Engine differences use an independently seeded finer-tolerance WoSt reference, not an exact solution.\n";
        std::cout
            << "kind,point,x,y,solver,reference,reference_se,mean,signed_difference,mc_se,combined_se,seconds,avg_steps,avg_queries,avg_center_checks,avg_cap_fallback,avg_corner_fallback,avg_radius_rejections,failures\n";
        PDE<float, 2> pde;
        setupPDE(pde, model, options.problem == "manufactured");
        BenchmarkSummary summaries[2];
        for (size_t i = 0; i < points.size(); ++i) {
            double referenceValue = evaluateExactSolution(points[i]);
            double referenceStandardError = 0;
            if (options.problem == "engine") {
                const auto result = runSolver(model, mesh, pde, points[i], options, false, true,
                                              options.seed + 1000000007ULL + i * 104729ULL);
                referenceValue = result.mean;
                referenceStandardError = result.standardError;
                writeResult("reference", int(i), points[i], "wost_reference", result, referenceValue,
                            referenceStandardError);
            }
            for (int method = 0; method < 2; ++method) {
                const auto result =
                    runSolver(model, mesh, pde, points[i], options, method == 0, false, options.seed + i * 104729ULL);
                writeResult(options.problem == "manufactured" ? "exact_error" : "reference_discrepancy", int(i),
                            points[i], method == 0 ? "construction4" : "wost", result, referenceValue,
                            referenceStandardError);
                summaries[method].add(result, referenceValue, referenceStandardError);
            }
        }
        summaries[0].print("construction4", int(points.size()));
        summaries[1].print("wost", int(points.size()));
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
