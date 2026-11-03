// This file instruments existing geometry callbacks for benchmark query counts.

#pragma once

#include <wosx/core/geometric_queries.h>
#include <atomic>
#include <cstdint>
#include <utility>

namespace walk_on_tangent_spheres {

// Instrument existing callbacks in place so the shared composite callbacks keep
// referring to the same object. No geometry is copied, cleaned or rebuilt here.
struct QueryStatistics {
    std::atomic<uint64_t> distance{0};
    std::atomic<uint64_t> projection{0};
    std::atomic<uint64_t> star{0};
    std::atomic<uint64_t> ray{0};
    std::atomic<uint64_t> sample{0};
    std::atomic<uint64_t> visibility{0};
    std::atomic<uint64_t> allHits{0};
    std::atomic<uint64_t> membership{0};

    void reset();

    uint64_t total() const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Implementation

inline void QueryStatistics::reset()
{
    distance = 0;
    projection = 0;
    star = 0;
    ray = 0;
    sample = 0;
    visibility = 0;
    allHits = 0;
    membership = 0;
}

inline uint64_t QueryStatistics::total() const
{
    // A membership callback's underlying all-hit ray searches are already
    // counted. Do not count the composite invocation again in this total.
    return distance + projection + star + ray + sample + visibility + allHits;
}

template <typename Callback>
auto countQuery(Callback callback, std::atomic<uint64_t>& count)
{
    return [callback, &count](auto&&... args) -> decltype(auto) {
        ++count;
        return callback(std::forward<decltype(args)>(args)...);
    };
}

inline void instrumentQueries(wosx::GeometricQueries<2>& queries, QueryStatistics& counts)
{
    queries.computeDistToAbsorbingBoundary = countQuery(queries.computeDistToAbsorbingBoundary, counts.distance);
    queries.computeDistToReflectingBoundary = countQuery(queries.computeDistToReflectingBoundary, counts.distance);
    queries.projectToAbsorbingBoundary = countQuery(queries.projectToAbsorbingBoundary, counts.projection);
    queries.projectToReflectingBoundary = countQuery(queries.projectToReflectingBoundary, counts.projection);
    queries.computeStarRadiusForReflectingBoundary =
        countQuery(queries.computeStarRadiusForReflectingBoundary, counts.star);
    queries.intersectAbsorbingBoundary = countQuery(queries.intersectAbsorbingBoundary, counts.ray);
    queries.intersectReflectingBoundary = countQuery(queries.intersectReflectingBoundary, counts.ray);
    queries.sampleReflectingBoundary = countQuery(queries.sampleReflectingBoundary, counts.sample);
    queries.intersectsWithReflectingBoundary = countQuery(queries.intersectsWithReflectingBoundary, counts.visibility);
    queries.intersectAbsorbingBoundaryAllHits = countQuery(queries.intersectAbsorbingBoundaryAllHits, counts.allHits);
    queries.intersectReflectingBoundaryAllHits = countQuery(queries.intersectReflectingBoundaryAllHits, counts.allHits);
    queries.insideDomain = countQuery(queries.insideDomain, counts.membership);
}

} // walk_on_tangent_spheres
