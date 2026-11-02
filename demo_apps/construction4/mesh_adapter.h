// This file supplies the prototype whole-boundary mesh queries and geometry
// checks used by the original tangent-sphere mesh benchmark.

#pragma once
#include "../basic_2d/model_problem.h"
#include <map>
#include <set>
#include <atomic>
#include <stdexcept>

namespace construction4_mesh {
using wosx::Vector2;
using wosx::Vector2i;
// Expose geometry only in this benchmark; leave the supplied loader unchanged.
class MeshModel : public ModelProblem<float> {
public:
    using ModelProblem<float>::ModelProblem;
    const std::vector<Vector2>& positions() const { return mPositions; }
    const std::vector<Vector2i>& indices() const { return mIndices; }
};

class MeshAdapter {
public:
    std::vector<Vector2> positions;
    std::vector<Vector2i> edges;
    std::vector<bool> neumann;
    std::vector<Vector2> normals;
    std::vector<int> incoming, outgoing;
    fcpw::Scene<2> scene;
    std::atomic<uint64_t> queryCount{0};
    int omittedZeroEdges = 0;
    double scale = 1;
    MeshAdapter(const MeshAdapter&) = delete;
    MeshAdapter& operator=(const MeshAdapter&) = delete;

    explicit MeshAdapter(MeshModel& model);
    static double cross(const Vector2& a, const Vector2& b);
    bool insideWinding(const Vector2& p) const;
    wosx::TangentBoundaryPoint2D nearest(const Vector2& p);
    wosx::TangentNearestBoundary2D callback();
    // Setup-only checks. Not charged to either solver's walk counts or timings.
    size_t validateQueries();

private:
    void validateNoCrossings() const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation

inline MeshAdapter::MeshAdapter(MeshModel& model)
{
    std::map<std::pair<float, float>, int> unique;
    std::vector<int> remap;
    for (const auto& p : model.positions()) {
        const auto key = std::make_pair(p.x(), p.y());
        auto found = unique.find(key);
        if (found == unique.end()) {
            const int id = int(positions.size());
            unique[key] = id;
            positions.push_back(p);
            remap.push_back(id);
        } else {
            remap.push_back(found->second);
        }
    }
    const auto& pde = model.getPDE();
    for (const auto& e : model.indices()) {
        const int a = remap.at(e.x()), b = remap.at(e.y());
        if (a == b) {
            ++omittedZeroEdges;
            continue;
        }
        if ((positions[a] - positions[b]).norm() < 1e-8f) {
            throw std::runtime_error("Unresolved tiny mesh edge");
        }
        edges.emplace_back(a, b);
        neumann.push_back(pde.hasReflectingBoundaryConditions((positions[a] + positions[b]) * .5f));
        Vector2 v = positions[b] - positions[a];
        normals.emplace_back(Vector2(v.y(), -v.x()).normalized());
    }
    if (edges.empty()) {
        throw std::runtime_error("Empty boundary mesh");
    }
    incoming.assign(positions.size(), -1);
    outgoing = incoming;
    for (size_t j = 0; j < edges.size(); ++j) {
        const auto& e = edges[j];
        if (outgoing[e.x()] != -1 || incoming[e.y()] != -1) {
            throw std::runtime_error("Mesh is not an oriented manifold after welding");
        }
        outgoing[e.x()] = int(j);
        incoming[e.y()] = int(j);
    }
    for (size_t i = 0; i < positions.size(); ++i) {
        if (incoming[i] < 0 || outgoing[i] < 0) {
            throw std::runtime_error("Boundary is not closed");
        }
    }
    Vector2 lo = positions[0], hi = lo;
    for (const auto& p : positions) {
        lo = lo.cwiseMin(p);
        hi = hi.cwiseMax(p);
    }
    scale = (hi - lo).norm();
    validateNoCrossings();
    // Validate every oriented edge using an independent winding-number predicate.
    for (size_t j = 0; j < edges.size(); ++j) {
        const auto& e = edges[j];
        Vector2 mid = (positions[e.x()] + positions[e.y()]) * .5f;
        const float offset = float(std::min(1e-5 * scale, .001 * double((positions[e.x()] - positions[e.y()]).norm())));
        if (!insideWinding(mid - offset * normals[j]) || insideWinding(mid + offset * normals[j])) {
            throw std::runtime_error("Boundary orientation/clearance could not be certified by winding tests");
        }
    }
    scene.setObjectCount(1);
    scene.setObjectVertices(positions, 0);
    scene.setObjectLineSegments(edges, 0);
    scene.computeObjectNormals(0, false); // unweighted pseudonormals for signed distance at vertices
    scene.build(fcpw::AggregateType::Bvh_SurfaceArea, false, false, false);
}

inline double MeshAdapter::cross(const Vector2& a, const Vector2& b)
{
    return double(a.x()) * b.y() - double(a.y()) * b.x();
}

inline bool MeshAdapter::insideWinding(const Vector2& p) const
{
    int winding = 0;
    for (const auto& e : edges) {
        const Vector2& a = positions[e.x()];
        const Vector2& b = positions[e.y()];
        const double side = cross(b - a, p - a);
        if (a.y() <= p.y()) {
            if (b.y() > p.y() && side > 0) {
                ++winding;
            }
        } else if (b.y() <= p.y() && side < 0) {
            --winding;
        }
    }
    return winding != 0;
}

inline wosx::TangentBoundaryPoint2D MeshAdapter::nearest(const Vector2& p)
{
    ++queryCount;
    fcpw::Interaction<2> i;
    if (!scene.findClosestPoint(p, i, fcpw::maxFloat, true) || i.primitiveIndex < 0 ||
        size_t(i.primitiveIndex) >= edges.size() || !std::isfinite(i.d)) {
        throw std::runtime_error("Whole-mesh nearest-point query failed");
    }
    const int j = i.primitiveIndex;
    wosx::TangentBoundaryPoint2D hit;
    hit.point = i.p;
    hit.normal = normals[j];
    hit.isNeumann = neumann[j];
    hit.signedClearance = -double(i.signedDistance(p));
    const float uv = i.uv.x();
    const int vertex = uv <= 1e-6f ? edges[j].x() : (uv >= 1 - 1e-6f ? edges[j].y() : -1);
    hit.normalIsUnique = true;
    if (vertex >= 0) {
        const int a = incoming[vertex], b = outgoing[vertex];
        hit.normalIsUnique = normals[a].dot(normals[b]) > 1 - 1e-6f && neumann[a] == neumann[b];
    }
    return hit;
}

inline wosx::TangentNearestBoundary2D MeshAdapter::callback()
{
    return [this](const Vector2& p) {
        return nearest(p);
    };
}

inline size_t MeshAdapter::validateQueries()
{
    Vector2 lo = positions[0], hi = lo;
    for (const auto& p : positions) {
        lo = lo.cwiseMin(p);
        hi = hi.cwiseMax(p);
    }
    lo -= Vector2::Constant(float(.03 * scale));
    hi += Vector2::Constant(float(.03 * scale));
    size_t tested = 0;
    auto check = [&](const Vector2& p, bool bruteForce) {
        const auto hit = nearest(p);
        if (std::abs(hit.signedClearance) > 2e-6 * scale && (hit.signedClearance > 0) != insideWinding(p)) {
            throw std::runtime_error("FCPW sign disagrees with winding predicate");
        }
        if (bruteForce) {
            double best = std::numeric_limits<double>::infinity();
            for (const auto& e : edges) {
                Eigen::Vector2d a = positions[e.x()].cast<double>(), b = positions[e.y()].cast<double>();
                const auto ab = (b - a).eval();
                const double t = std::clamp((p.cast<double>() - a).dot(ab) / ab.squaredNorm(), 0., 1.);
                best = std::min(best, (p.cast<double>() - a - t * ab).norm());
            }
            if (std::abs(std::abs(hit.signedClearance) - best) > 2e-6 * scale) {
                throw std::runtime_error("FCPW nearest distance disagrees with brute-force segments");
            }
        }
        ++tested;
    };
    for (int j = 0; j < 41; ++j) {
        for (int i = 0; i < 41; ++i) {
            check(Vector2(lo.x() + (i + .37f) / 41 * (hi.x() - lo.x()), lo.y() + (j + .61f) / 41 * (hi.y() - lo.y())),
                  true);
        }
    }
    for (size_t j = 0; j < positions.size(); ++j) {
        const float eps = float(2e-5 * scale);
        for (int k = 0; k < 8; ++k) {
            const double th = 2 * wosx::tangent_disk::pi * (k + .25) / 8;
            check(positions[j] + eps * Vector2(float(std::cos(th)), float(std::sin(th))), false);
        }
    }
    queryCount = 0;
    return tested;
}

inline void MeshAdapter::validateNoCrossings() const
{
    for (size_t i = 0; i < edges.size(); ++i) {
        for (size_t j = i + 1; j < edges.size(); ++j) {
            const auto a = edges[i], b = edges[j];
            if (a.x() == b.x() || a.x() == b.y() || a.y() == b.x() || a.y() == b.y()) {
                continue;
            }
            const auto p = positions[a.x()], q = positions[a.y()], r = positions[b.x()], s = positions[b.y()];
            const double c1 = cross(q - p, r - p), c2 = cross(q - p, s - p), c3 = cross(s - r, p - r),
                         c4 = cross(s - r, q - r);
            if (c1 * c2 < 0 && c3 * c4 < 0) {
                throw std::runtime_error("Self-intersecting mesh boundary");
            }
            // Also reject non-adjacent contact/overlap within a tight geometry tolerance.
            const double eps = 1e-12 * scale * scale;
            auto on = [&](const Vector2& x, const Vector2& u, const Vector2& v, double c) {
                return std::abs(c) < eps && x.x() >= std::min(u.x(), v.x()) && x.x() <= std::max(u.x(), v.x()) &&
                       x.y() >= std::min(u.y(), v.y()) && x.y() <= std::max(u.y(), v.y());
            };
            if (on(r, p, q, c1) || on(s, p, q, c2) || on(p, r, s, c3) || on(q, r, s, c4)) {
                throw std::runtime_error("Non-manifold boundary contact/overlap");
            }
        }
    }
}

struct WoStCounts {
    uint64_t distance = 0;
    uint64_t projection = 0;
    uint64_t star = 0;
    uint64_t ray = 0;
    uint64_t sample = 0;
    uint64_t visibility = 0;
    uint64_t total() const { return distance + projection + star + ray + sample + visibility; }
};
// Wrap only in the new benchmark; existing WoSt and FCPW source is unchanged.
inline wosx::GeometricQueries<2> countedQueries(const wosx::GeometricQueries<2>& original, WoStCounts& c)
{
    auto q = original;
    auto d = original.computeDistToAbsorbingBoundary;
    q.computeDistToAbsorbingBoundary = [d, &c](const Vector2& p, bool s) {
        ++c.distance;
        return d(p, s);
    };
    auto dn = original.computeDistToReflectingBoundary;
    q.computeDistToReflectingBoundary = [dn, &c](const Vector2& p, bool s) {
        ++c.distance;
        return dn(p, s);
    };
    auto proj = original.projectToAbsorbingBoundary;
    q.projectToAbsorbingBoundary = [proj, &c](Vector2& p, Vector2& n, float& d, bool s) {
        ++c.projection;
        return proj(p, n, d, s);
    };
    auto star = original.computeStarRadiusForReflectingBoundary;
    q.computeStarRadiusForReflectingBoundary = [star, &c](const Vector2& p, float a, float b, float e, bool f) {
        ++c.star;
        return star(p, a, b, e, f);
    };
    auto ray = original.intersectReflectingBoundary;
    q.intersectReflectingBoundary = [ray, &c](const Vector2& p, const Vector2& n, const Vector2& v, float r, bool on,
                                              wosx::IntersectionPoint<2>& hit) {
        ++c.ray;
        return ray(p, n, v, r, on, hit);
    };
    auto sample = original.sampleReflectingBoundary;
    q.sampleReflectingBoundary = [sample, &c](const Vector2& p, float r, const Vector2& u,
                                              wosx::BoundarySample<2>& hit) {
        ++c.sample;
        return sample(p, r, u, hit);
    };
    auto vis = original.intersectsWithReflectingBoundary;
    q.intersectsWithReflectingBoundary = [vis, &c](const Vector2& p, const Vector2& v, const Vector2& n,
                                                   const Vector2& m, bool a, bool b) {
        ++c.visibility;
        return vis(p, v, n, m, a, b);
    };
    return q;
}
} // construction4_mesh
