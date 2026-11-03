// This file defines analytic rectangle geometry and manufactured PDE data for
// CPU solver comparisons.

#pragma once

#include <wosx/point_estimation/walk_on_tangent_spheres.h>
#include <algorithm>
#include <atomic>

namespace walk_on_tangent_spheres {
struct QueryCounts {
    std::atomic<uint64_t> absorbingDistance{0};
    std::atomic<uint64_t> reflectingDistance{0};
    std::atomic<uint64_t> absorbingProjection{0};
    std::atomic<uint64_t> reflectingProjection{0};
    std::atomic<uint64_t> ray{0};
    std::atomic<uint64_t> boundarySample{0};
    std::atomic<uint64_t> starRadius{0};
    std::atomic<uint64_t> membership{0};
};

// All callbacks capture this object: do not copy/move it after construction.
// Rectangle (-1,1)x(0,1), bottom Neumann, other sides Dirichlet.
struct Rectangle {
    // constructors
    Rectangle();
    Rectangle(const Rectangle&) = delete;
    Rectangle& operator=(const Rectangle&) = delete;

    // geometric queries
    static double clearance(const wosx::Vector2& pt);
    static bool contains(const wosx::Vector2& pt);
    static wosx::Vector2 nearestAbsorbingBoundaryPoint(const wosx::Vector2& pt);

    // members
    wosx::GeometricQueries<2> queries;
    QueryCounts counts;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation

inline double Rectangle::clearance(const wosx::Vector2& pt)
{
    return std::min({double(pt.x()) + 1, 1 - double(pt.x()), double(pt.y()), 1 - double(pt.y())});
}

inline bool Rectangle::contains(const wosx::Vector2& pt)
{
    return clearance(pt) >= 0;
}

inline wosx::Vector2 Rectangle::nearestAbsorbingBoundaryPoint(const wosx::Vector2& pt)
{
    const float x = std::clamp(pt.x(), -1.f, 1.f), z = std::clamp(pt.y(), 0.f, 1.f);
    wosx::Vector2 a(-1, z), b(1, z), c(x, 1);
    auto hit = a;
    if ((b - pt).squaredNorm() < (hit - pt).squaredNorm()) {
        hit = b;
    }
    if ((c - pt).squaredNorm() < (hit - pt).squaredNorm()) {
        hit = c;
    }
    return hit;
}

inline Rectangle::Rectangle()
{
    using namespace wosx;
    queries.domainIsWatertight = true;
    queries.hasNonEmptyAbsorbingBoundary = true;
    queries.hasNonEmptyReflectingBoundary = true;
    queries.domainMin = Vector2(-1, 0);
    queries.domainMax = Vector2(1, 1);
    queries.computeDistToAbsorbingBoundary = [this](const Vector2& pt, bool signedDistance) {
        ++counts.absorbingDistance;
        const float d = (pt - nearestAbsorbingBoundaryPoint(pt)).norm();
        return signedDistance && contains(pt) ? -d : d;
    };
    queries.computeDistToReflectingBoundary = [this](const Vector2& pt, bool signedDistance) {
        ++counts.reflectingDistance;
        const Vector2 hit(std::clamp(pt.x(), -1.f, 1.f), 0);
        const float d = (pt - hit).norm();
        return signedDistance && pt.y() >= 0 ? -d : d;
    };
    queries.projectToAbsorbingBoundary = [this](Vector2& pt, Vector2& normal, float& d, bool computeSignedDistance) {
        ++counts.absorbingProjection;
        const Vector2 queryPt = pt;
        pt = nearestAbsorbingBoundaryPoint(queryPt);
        normal = pt.y() == 1 ? Vector2(0, 1) : (pt.x() < 0 ? Vector2(-1, 0) : Vector2(1, 0));
        d = (queryPt - pt).norm();
        if (computeSignedDistance && contains(queryPt)) {
            d = -d;
        }
        return true;
    };
    queries.projectToReflectingBoundary = [this](Vector2& pt, Vector2& normal, float& d, bool computeSignedDistance) {
        ++counts.reflectingProjection;
        const Vector2 queryPt = pt;
        pt = Vector2(std::clamp(pt.x(), -1.f, 1.f), 0);
        normal = Vector2(0, -1);
        d = (queryPt - pt).norm();
        if (computeSignedDistance && queryPt.y() >= 0) {
            d = -d;
        }
        return true;
    };
    queries.insideDomain = [this](const Vector2& pt) {
        ++counts.membership;
        return contains(pt);
    };
    // The rectangle intersected with any disk is convex: no silhouette bound needed.
    queries.computeStarRadiusForReflectingBoundary = [this](const Vector2&, float, float maxRadius, float, bool) {
        ++counts.starRadius;
        return maxRadius;
    };
    queries.intersectReflectingBoundary = [this](const Vector2& origin, const Vector2& normal, const Vector2& direction,
                                                 float maxRadius, bool onBoundary, IntersectionPoint<2>& hit) {
        ++counts.ray;
        Vector2 pt = onBoundary ? queries.offsetPointAlongDirection(origin, -normal) : origin;
        if (direction.y() >= 0) {
            return false;
        }
        const float t = -pt.y() / direction.y();
        const float x = pt.x() + t * direction.x();
        if (t < 0 || t > maxRadius || x < -1 || x > 1) {
            return false;
        }
        hit = IntersectionPoint<2>(Vector2(x, 0), Vector2(0, -1), t);
        return true;
    };
    queries.sampleReflectingBoundary = [this](const Vector2& pt, float r, const Vector2& u, BoundarySample<2>& hit) {
        ++counts.boundarySample;
        if (std::abs(pt.y()) >= r) {
            return false;
        }
        const double h = std::sqrt(double(r) * r - double(pt.y()) * pt.y());
        const double lo = std::max(-1., double(pt.x()) - h), hi = std::min(1., double(pt.x()) + h);
        if (!(hi > lo)) {
            return false;
        }
        hit = BoundarySample<2>(Vector2(float(lo + (hi - lo) * u.x()), 0), Vector2(0, -1), float(1 / (hi - lo)));
        return true;
    };
    // All points used here are in the rectangle and it is convex.
    queries.intersectsWithReflectingBoundary = [](const Vector2&, const Vector2&, const Vector2&, const Vector2&, bool,
                                                  bool) {
        return false;
    };
}

enum class Problem { Polynomial, Harmonic, Constant };
inline double exact(Problem problem, const wosx::Vector2& pt)
{
    if (problem == Problem::Constant) {
        return 1;
    }
    if (problem == Problem::Harmonic) {
        return 1 + .3 * pt.y() + .2 * (double(pt.x()) * pt.x() - double(pt.y()) * pt.y());
    }
    return 1 + .3 * pt.y() + .2 * double(pt.x()) * pt.x() + .1 * double(pt.y()) * pt.y();
}

inline void setupPDE(wosx::PDE<float, 2>& pde, Problem problem)
{
    pde.absorptionCoeff = 0;
    pde.isSourceConstant = true;
    pde.areRobinConditionsPureNeumann = true;
    // WoSX -Delta u=s. Handoff polynomial has Delta u=+0.6, hence s=-0.6.
    pde.source = [problem](const wosx::Vector2&) {
        return problem == Problem::Polynomial ? -.6f : 0.f;
    };
    pde.dirichlet = [problem](const wosx::Vector2& pt, bool) {
        return float(exact(problem, pt));
    };
    pde.robin = [problem](const wosx::Vector2&, const wosx::Vector2&, bool) {
        return problem == Problem::Constant ? 0.f : -.3f;
    };
    pde.robinCoeff = [](const wosx::Vector2&, const wosx::Vector2&, bool) {
        return 0.f;
    };
    pde.hasReflectingBoundaryConditions = [](const wosx::Vector2& pt) {
        return std::abs(pt.y()) < 1e-6f;
    };
}

inline wosx::SamplePoint<float, 2> samplePoint(uint64_t seed, const wosx::Vector2& pt = wosx::Vector2(0, .1f))
{
    wosx::SamplePoint<float, 2> sample(pt, wosx::Vector2::Zero(), wosx::SampleType::InDomain,
                                       wosx::EstimationQuantity::Solution, 1,
                                       (pt - Rectangle::nearestAbsorbingBoundaryPoint(pt)).norm(), std::abs(pt.y()));
    sample.rng = pcg32(seed, 1);
    return sample;
}
} // walk_on_tangent_spheres
