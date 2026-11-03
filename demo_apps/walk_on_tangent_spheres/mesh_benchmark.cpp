// This file compares CPU solvers on the supplied 2D mesh problems using the
// existing WoSX geometry and PDE callbacks.

#include "../basic_2d/model_problem.h"
#include "query_statistics.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <chrono>

namespace {
using namespace wosx;
using namespace walk_on_tangent_spheres;
struct BenchmarkResult {
    // members
    double mean = 0;
    double standardError = 0;
    double seconds = 0;
    uint64_t walks = 0;
    uint64_t steps = 0;
    uint64_t queries = 0;
    uint64_t centers = 0;
    uint64_t cap = 0;
    uint64_t memberships = 0;
    uint64_t rejections = 0;
    uint64_t failures = 0;
};

struct BenchmarkOptions {
    // members
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
};
inline double evaluateExactSolution(const Vector2& pt)
{
    return 1 + .3 * pt.y() + .2 * double(pt.x()) * pt.x() + .1 * double(pt.y()) * pt.y();
}
void setupPDE(PDE<float, 2>& pde, ModelProblem<float>& model, bool manufactured)
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
BenchmarkResult runSolver(ModelProblem<float>& model,
                          QueryStatistics& counts,
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
    const auto& queries = original;
    BenchmarkResult result;
    result.walks = walks;
    counts.reset();
    const auto start = std::chrono::steady_clock::now();
    if (useTangentSpheres) {
        TangentSphereSettings tangentSettings;
        tangentSettings.capChordRadius = options.cap;
        tangentSettings.initialTangentRadius = options.radius;
        WalkOnTangentSpheres<float, 2> solver(queries);
        try {
            solver.solve(pde, walkSettings, tangentSettings, walks, sample, statistics, diagnostics);
        } catch (const std::exception& error) {
            throw std::runtime_error(
                std::string("Walk on Tangent Spheres failed: ") + error.what() + "; evaluation_point=(" +
                std::to_string(pt.x()) + "," + std::to_string(pt.y()) + ")" + "; seed=" + std::to_string(seed) +
                "; completed=" + std::to_string(diagnostics.walksCompleted) +
                " attempted=" + std::to_string(diagnostics.walksAttempted) + "; no failed-batch estimate reported");
        }
        result.queries = counts.total();
        result.steps = diagnostics.interiorSteps + diagnostics.tangentSteps;
        result.centers = diagnostics.centerCertificates;
        result.cap = diagnostics.capAdjustments;
        result.memberships = counts.membership;
        result.rejections = diagnostics.radiusRejections;
        if (diagnostics.nearestPointQueries != result.steps + result.walks ||
            counts.projection != 2 * diagnostics.nearestPointQueries ||
            counts.distance != 2 * diagnostics.centerCertificates ||
            counts.membership != diagnostics.initialMembershipChecks + diagnostics.centerMembershipChecks) {
            throw std::runtime_error("Shared geometry query accounting mismatch");
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
        result.memberships = counts.membership;
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
              << double(result.cap) / result.walks << ',' << double(result.memberships) / result.walks << ','
              << double(result.rejections) / result.walks << ',' << result.failures << '\n';
}

struct BenchmarkSummary {
    // accumulated statistics
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
    uint64_t memberships = 0;
    uint64_t rejections = 0;
    uint64_t failures = 0;

    // adds one evaluation point to the summary
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
        memberships += result.memberships;
        rejections += result.rejections;
        failures += result.failures;
    }
    // writes the aggregate comparison statistics
    void print(const char* solver, int n) const
    {
        std::cout << "SUMMARY," << solver << ",points=" << n << ",walks=" << walks
                  << ",rmse=" << std::sqrt(squaredErrorSum / n) << ",max_abs_error=" << maxAbsoluteError
                  << ",rms_mc_se=" << std::sqrt(varianceSum / n)
                  << ",rms_reference_se=" << std::sqrt(referenceVarianceSum / n) << ",seconds=" << seconds
                  << ",avg_steps=" << double(steps) / walks << ",avg_queries=" << double(queries) / walks
                  << ",avg_center_checks=" << double(centers) / walks << ",avg_cap_fallback=" << double(cap) / walks
                  << ",avg_membership_checks=" << double(memberships) / walks
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
            if (key == "--help") {
                std::cout
                    << "walk_on_tangent_spheres_mesh_benchmark --case manufactured|engine [--root repo] [--config path]\n"
                    << " [--points 16] [--walks 6000] [--seed 20261030] [--cap .01] [--radius .4]\n"
                    << " [--tol-dirichlet .01] [--tol-neumann .01] [--max-steps 100000]\n"
                    << " [--reference-walks 65536] [--reference-tolerance .0001]\n";
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
        ModelProblem<float> model(modelProblemConfig, root.string());
        const auto& queries = model.getGeometricQueries();
        const auto box = model.getBoundingBox();
        const double scale = (box.second - box.first).norm();
        std::vector<Vector2> candidates;
        for (int j = 0; j < 64; ++j) {
            for (int i = 0; i < 64; ++i) {
                Vector2 pt = box.first + Vector2((i + .5f) / 64 * (box.second.x() - box.first.x()),
                                                 (j + .5f) / 64 * (box.second.y() - box.first.y()));
                if (queries.insideDomain(pt) && queries.computeDistToBoundary(pt, false) > .02 * scale) {
                    candidates.push_back(pt);
                }
            }
        }
        QueryStatistics counts;
        instrumentQueries(model.getGeometricQueries(), counts);
        std::cout << std::setprecision(10) << "# mesh=" << modelProblemConfig.at("geometry")
                  << " shared_WoSX_geometry=true approach_direction_normal=true\n";
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
            << "# Both totals count underlying distance,projection,star,ray,sample,visibility,all-hit searches; membership calls reported separately.\n"
            << "# Single-threaded timings exclude loading, shared BVH builds, reference and initial point metadata.\n"
            << "# Engine differences use an independently seeded finer-tolerance WoSt reference, not an exact solution.\n";
        std::cout
            << "kind,point,x,y,solver,reference,reference_se,mean,signed_difference,mc_se,combined_se,seconds,avg_steps,avg_queries,avg_center_checks,avg_cap_fallback,avg_membership_checks,avg_radius_rejections,failures\n";
        PDE<float, 2> pde;
        setupPDE(pde, model, options.problem == "manufactured");
        BenchmarkSummary summaries[2];
        for (size_t i = 0; i < points.size(); ++i) {
            double referenceValue = evaluateExactSolution(points[i]);
            double referenceStandardError = 0;
            if (options.problem == "engine") {
                const auto result = runSolver(model, counts, pde, points[i], options, false, true,
                                              options.seed + 1000000007ULL + i * 104729ULL);
                referenceValue = result.mean;
                referenceStandardError = result.standardError;
                writeResult("reference", int(i), points[i], "wost_reference", result, referenceValue,
                            referenceStandardError);
            }
            for (int method = 0; method < 2; ++method) {
                const auto result =
                    runSolver(model, counts, pde, points[i], options, method == 0, false, options.seed + i * 104729ULL);
                writeResult(options.problem == "manufactured" ? "exact_error" : "reference_discrepancy", int(i),
                            points[i], method == 0 ? "walk_on_tangent_spheres" : "wost", result, referenceValue,
                            referenceStandardError);
                summaries[method].add(result, referenceValue, referenceStandardError);
            }
        }
        summaries[0].print("walk_on_tangent_spheres", int(points.size()));
        summaries[1].print("wost", int(points.size()));
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "ERROR: " << error.what() << '\n';
        return 1;
    }
}
