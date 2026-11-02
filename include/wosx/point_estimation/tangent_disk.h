// This file provides the explicit 2D tangent-disk flux, mass and sampling formulas.

#pragma once

#include <cmath>
#include <stdexcept>

namespace wosx {
namespace tangent_disk {

// These are explicit 2D disk calculations, not claims about global convergence.
constexpr double pi = 3.141592653589793238462643383279502884;

inline void validate(double r, double d)
{
    if (!std::isfinite(r) || !std::isfinite(d) || !(r > 0 && d > 0 && d < 2 * r)) {
        throw std::invalid_argument("Tangent disk requires finite r > 0 and 0 < d < 2r");
    }
}

inline double cotHalfCap(double r, double d)
{
    validate(r, d);
    const double h = d / (2 * r);
    return std::sqrt((1 - h) * (1 + h)) / h;
}

inline double complementaryMass(double r, double d)
{
    return cotHalfCap(r, d) / (pi * r);
}
// Density per unit arc length. Stable near theta=0 for positive eta.
inline double finiteFlux(double r, double eta, double theta)
{
    if (!(r > 0 && eta > 0 && eta < 2 * r)) {
        throw std::invalid_argument("Finite flux requires 0 < eta < 2r");
    }
    const double s = std::sin(theta / 2);
    return eta * (2 * r - eta) / (2 * pi * r * (eta * eta + 4 * r * (r - eta) * s * s));
}

inline double limitingFlux(double r, double theta)
{
    const double chord = 2 * r * std::sin(theta / 2);
    return 1 / (pi * chord * chord);
}
// u and v are independent uniforms in [0,1). No evaluation of singular flux.
inline double sampleAngle(double r, double d, double u, double v)
{
    if (!(u >= 0 && u < 1 && v >= 0 && v < 1)) {
        throw std::invalid_argument("Sampler requires uniforms in [0,1)");
    }
    const double theta = 2 * std::atan2(1.0, (1 - u) * cotHalfCap(r, d));
    return v < .5 ? theta : -theta;
}
// Handoff convention Delta u=f, G<=0. Exact scalar volume formula and limit.
inline double finiteConstantVolume(double r, double eta, double f)
{
    return f * eta * (eta - 2 * r) / 4;
}

inline double limitingConstantVolume(double r, double f)
{
    return -f * r / 2;
}

} // tangent_disk
} // wosx
