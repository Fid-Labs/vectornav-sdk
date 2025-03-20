// The MIT License (MIT)
// 
//  VectorNav Software Development Kit (v0.14.2)
// Copyright (c) 2024 VectorNav Technologies, LLC
// 
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#pragma once

#ifndef CORE_CONVERSIONS_HPP
#define CORE_CONVERSIONS_HPP

#define _USE_MATH_DEFINES

#include <vector>
#include <math.h>
#include <cmath>
#include "TemplateLibrary/Matrix.hpp"
#include "../../../tests/GTestUtils.hpp"
#include "Implementation/MeasurementDatatypes.hpp"

namespace VN
{
namespace Conversions
{

// ******
// Angles
// ******

inline float rad2deg(float angleInRads) noexcept { return angleInRads * 180.0f / static_cast<float>(M_PI); }

inline double rad2deg(double angleInRads) noexcept { return angleInRads * 180.0 / static_cast<double>(M_PI); }

inline float deg2rad(float angleInDegs) noexcept { return angleInDegs * static_cast<float>(M_PI) / 180.0f; }

inline double deg2rad(double angleInDegs) noexcept { return angleInDegs * static_cast<double>(M_PI) / 180.0; }

inline Ypr deg2rad(Ypr ypr) noexcept
{
    ypr.yaw *= (static_cast<float>(M_PI) / 180.0f);
    ypr.pitch *= (static_cast<float>(M_PI) / 180.0f);
    ypr.roll *= (static_cast<float>(M_PI) / 180.0f);
    return ypr;
}

inline Ypr rad2deg(Ypr ypr) noexcept
{
    ypr.yaw *= (180.0f / static_cast<float>(M_PI));
    ypr.pitch *= (180.0f / static_cast<float>(M_PI));
    ypr.roll *= (180.0f / static_cast<float>(M_PI));
    return ypr;
}

// ***********
// Temperature
// ***********

inline float celsius2fahren(float tempInCelsius) noexcept { return (tempInCelsius * 9.0f) / 5.0f + 32.0f; }

inline double celsius2fahren(double tempInCelsius) noexcept { return (tempInCelsius * 9.0) / 5.0 + 32.0; }

inline float fahren2celsius(float tempInFahren) noexcept { return (tempInFahren - 32.0f) * 5.0f / 9.0f; }

inline double fahren2celsius(double tempInFahren) noexcept { return (tempInFahren - 32.0) * 5.0 / 9.0; }

inline float celsius2kelvin(float tempInCelsius) noexcept { return tempInCelsius + 273.15f; }

inline double celsius2kelvin(double tempInCelsius) noexcept { return tempInCelsius + 273.15; }

inline float kelvin2celsius(float tempInKelvin) noexcept { return tempInKelvin - 273.15f; }

inline double kelvin2celsius(double tempInKelvin) noexcept { return tempInKelvin - 273.15; }

inline float fahren2kelvin(float tempInFahren) noexcept { return celsius2kelvin(fahren2celsius(tempInFahren)); }

inline double fahren2kelvin(double tempInFahren) noexcept { return celsius2kelvin(fahren2celsius(tempInFahren)); }

inline float kelvin2fahren(float tempInKelvin) noexcept { return celsius2fahren(tempInKelvin - 273.15f); }

inline double kelvin2fahren(double tempInKelvin) noexcept { return celsius2fahren(tempInKelvin - 273.15); }

// ********
// Attitude
// ********

inline Quat yprInRads2Quat(const Ypr& ypr) noexcept
{
    float c1 = std::cos(ypr.yaw / 2.0f);
    float s1 = std::sin(ypr.yaw / 2.0f);
    float c2 = std::cos(ypr.pitch / 2.0f);
    float s2 = std::sin(ypr.pitch / 2.0f);
    float c3 = std::cos(ypr.roll / 2.0f);
    float s3 = std::sin(ypr.roll / 2.0f);

    return Quat({c1 * c2 * s3 - s1 * s2 * c3, c1 * s2 * c3 + s1 * c2 * s3, s1 * c2 * c3 - c1 * s2 * s3}, c1 * c2 * c3 + s1 * s2 * s3);
}

inline Quat yprInDegs2Quat(const Ypr& ypr) noexcept { return yprInRads2Quat(deg2rad(ypr)); }

inline Mat3f yprInRads2Dcm(const Ypr& ypr) noexcept
{
    float st1 = std::sin(ypr.yaw);
    float ct1 = std::cos(ypr.yaw);
    float st2 = std::sin(ypr.pitch);
    float ct2 = std::cos(ypr.pitch);
    float st3 = std::sin(ypr.roll);
    float ct3 = std::cos(ypr.roll);

#if PROPOSED_NEW_FORMUALA
    // PL-175 - YPR to DCM conversion error
    return Mat3f({ct2 * ct1, -ct2 * st1, st2, st3 * st2 * ct1 + ct3 * st1, -st3 * st2 * st1 + ct3 * ct1, -st3 * ct2, -ct3 * st2 * ct1 + st3 * st1,
                  ct3 * st2 * st1 + st3 * ct1, ct3 * ct2});
#else
    return Mat3f({ct2 * ct1, st3 * st2 * ct1 - ct3 * st1, ct3 * st2 * ct1 + st3 * st1, ct2 * st1, st3 * st2 * st1 + ct3 * ct1, ct3 * st2 * st1 - st3 * ct1,
                  -st2, st3 * ct2, ct3 * ct2});
#endif
}

inline Mat3f yprInDegs2Dcm(const Ypr& ypr) noexcept { return yprInRads2Dcm(deg2rad(ypr)); }

inline Ypr quat2YprInRads(const Quat& quat) noexcept
{
    float q1 = quat.vector[0];
    float q2 = quat.vector[1];
    float q3 = quat.vector[2];
    float q0 = quat.scalar;

    return Ypr(std::atan2(2.0f * (q1 * q2 + q0 * q3), q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3), std::asin(-2.0f * (q1 * q3 - q0 * q2)),
               std::atan2(2.0f * (q2 * q3 + q0 * q1), q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3));
}

inline Ypr quat2YprInDegs(const Quat& quat) noexcept { return rad2deg(quat2YprInRads(quat)); }

inline Mat3f quat2dcm(const Quat& quat) noexcept
{
    float q1 = quat.vector[0];
    float q2 = quat.vector[1];
    float q3 = quat.vector[2];
    float q0 = quat.scalar;

    return Mat3f({q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3, 2.f * (q1 * q2 - q0 * q3), 2.f * (q1 * q3 + q0 * q2), 2.f * (q1 * q2 + q0 * q3),
                  q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3, 2.f * (q2 * q3 - q0 * q1), 2.f * (q1 * q3 - q0 * q2), 2.f * (q2 * q3 + q0 * q1),
                  q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3});
}

inline Ypr dcm2YprInRads(const Mat3f& dcm) noexcept { return Ypr(std::atan2(dcm(1, 0), dcm(0, 0)), -(std::asin(dcm(2, 0))), std::atan2(dcm(2, 1), dcm(2, 2))); }

inline Ypr dcm2YprInDegs(const Mat3f& dcm) noexcept { return rad2deg(dcm2YprInRads(dcm)); }

inline Quat dcm2quat(const Mat3f& dcm) noexcept
{
    float tr = dcm(0, 0) + dcm(1, 1) + dcm(2, 2);
    Quat q;

    if (tr > 0)
    {
        float S = std::sqrt(tr + 1.0f) * 2.0f;
        q.scalar = 0.25f * S;
        q.vector[0] = (dcm(2, 1) - dcm(1, 2)) / S;
        q.vector[1] = (dcm(0, 2) - dcm(2, 0)) / S;
        q.vector[2] = (dcm(1, 0) - dcm(0, 1)) / S;
    }
    else if ((dcm(0, 0) > dcm(1, 1)) && (dcm(0, 0) > dcm(2, 2)))
    {
        float S = std::sqrt(1.0f + dcm(0, 0) - dcm(1, 1) - dcm(2, 2) * 2.0f);
        q.scalar = (dcm(2, 1) - dcm(1, 2)) / S;
        q.vector[0] = 0.25f * S;
        q.vector[1] = (dcm(0, 1) + dcm(1, 0)) / S;
        q.vector[2] = (dcm(0, 2) + dcm(2, 0)) / S;
    }
    else if (dcm(1, 1) > dcm(2, 2))
    {
        float S = std::sqrt(1.0f + dcm(1, 1) - dcm(0, 0) - dcm(2, 2) * 2.0f);
        q.scalar = (dcm(0, 2) - dcm(2, 0)) / S;
        q.vector[0] = (dcm(0, 1) + dcm(1, 0)) / S;
        q.vector[1] = 0.25f * S;
        q.vector[2] = (dcm(1, 2) + dcm(2, 1)) / S;
    }
    else
    {
        float S = std::sqrt(1.0f + dcm(2, 2) - dcm(0, 0) - dcm(1, 1) * 2.0f);
        q.scalar = (dcm(1, 0) - dcm(0, 1)) / S;
        q.vector[0] = (dcm(0, 2) + dcm(2, 0)) / S;
        q.vector[1] = (dcm(1, 2) + dcm(2, 1)) / S;
        q.vector[2] = 0.25f * S;
    }

    float norm = std::sqrt(q.vector[0] * q.vector[0] + q.vector[1] * q.vector[1] + q.vector[2] * q.vector[2] + q.scalar * q.scalar);
    q.scalar /= norm;
    q.vector[0] /= norm;
    q.vector[1] /= norm;
    q.vector[2] /= norm;

    return q;
}

inline float course_over_ground(float velNedX, float velNedY) noexcept { return std::atan2(velNedY, velNedX); }

inline float course_over_ground(const Vec3f& velNed) noexcept { return course_over_ground(velNed[0], velNed[1]); }

inline float speed_over_ground(const float velNedX, const float velNedY) noexcept { return std::sqrt(velNedX * velNedX + velNedY * velNedY); }

inline float speed_over_ground(const Vec3f& velNed) noexcept { return speed_over_ground(velNed[0], velNed[1]); }

constexpr double C_E2 = 0.006694379990141;
constexpr double C_EPSILON = 0.996647189335253;
constexpr double C_ABAR = 42.69767270717997;
constexpr double C_BBAR = 42.84131151331357;
constexpr double C_A = 6378.137;

inline Lla ecef_to_lla_v3d(const Vec3d& ecef) noexcept
{
    double x, y, z, r, c_phi, c_phi0, s_phi, s_phi0, tau, lat, lon, /*alt,*/ eta, h;

    const double Rthresh = 0.001;  // Limit on distance from pole in km to switch calculation.

    x = ecef[0];
    y = ecef[1];
    z = ecef[2];

    r = std::sqrt(x * x + y * y);

    if (r < Rthresh)
    {
        c_phi0 = 0;
        s_phi0 = (z > 0) - (z < 0);  // computes sign
    }
    else
    {
        double tau0 = z / (C_EPSILON * r);
        c_phi0 = 1 / std::sqrt(1 + tau0 * tau0);
        s_phi0 = tau0 * c_phi0;
    }

    tau = (z + C_BBAR * s_phi0 * s_phi0 * s_phi0) / (r - C_ABAR * c_phi0 * c_phi0 * c_phi0);
    lat = std::atan(tau);

    if (r < Rthresh)
    {
        c_phi = 0;
        s_phi = (z > 0) - (z < 0);  // computes sign
    }
    else
    {
        c_phi = 1 / std::sqrt(1 + tau * tau);
        s_phi = tau * c_phi;
    }

    eta = std::sqrt(1 - C_E2 * s_phi * s_phi);
    h = r * c_phi + z * s_phi - C_A * eta;

    lon = std::atan2(y, x);

    // Convert latitude and longitude to degrees and altitude to meters
    return Lla(lat * 180 / static_cast<double>(M_PI), lon * 180 / static_cast<double>(M_PI), h * 1000);
}

inline Vec3d lla_to_ecef_v3d(Lla lla) noexcept
{
    double n, x, y, z;
    double t1; /* TEMPS */

    lla.lat *= static_cast<double>(M_PI) / 180.0f;  // Geodetic latitude in radians.
    lla.lon *= static_cast<double>(M_PI) / 180.0f;  // Longitude in radians.
    lla.alt /= 1000;                                // Altitude above WGS84 in km.

    t1 = std::sin(lla.lat);
    n = C_A / std::sqrt(1 - C_E2 * t1 * t1);

    t1 = lla.alt + n;
    x = t1 * std::cos(lla.lat) * std::cos(lla.lon);
    y = t1 * std::cos(lla.lat) * std::sin(lla.lon);
    z = (t1 - C_E2 * n) * std::sin(lla.lat);

    return Vec3d({x, y, z});
}

}  // namespace Conversions
}  // namespace VN

#endif  // CORE_CONVERSIONS_HPP
