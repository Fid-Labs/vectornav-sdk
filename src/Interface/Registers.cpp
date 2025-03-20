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


#include "Interface/Registers.hpp"
#include <cstdint>
#include "TemplateLibrary/String.hpp"

namespace VN
{
namespace Registers
{
Vector<const char*, Config::PacketFinders::asciiMaxFieldCount> findIndexOfFieldSeparators(const AsciiMessage& input);
namespace Attitude
{
bool YawPitchRoll::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 3 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto yaw_t = StringUtils::fromString<float>(start, end);
    if (!yaw_t.has_value()) return true;
    yaw = yaw_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto pitch_t = StringUtils::fromString<float>(start, end);
    if (!pitch_t.has_value()) return true;
    pitch = pitch_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto roll_t = StringUtils::fromString<float>(start, end);
    if (!roll_t.has_value()) return true;
    roll = roll_t.value();

    return false;  // no errors occurred
}

bool operator==(const YawPitchRoll& lhs, const YawPitchRoll& rhs) { return ((lhs.yaw == rhs.yaw) && (lhs.pitch == rhs.pitch) && (lhs.roll == rhs.roll)); }

bool Quaternion::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 4 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto quatX_t = StringUtils::fromString<float>(start, end);
    if (!quatX_t.has_value()) return true;
    quatX = quatX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto quatY_t = StringUtils::fromString<float>(start, end);
    if (!quatY_t.has_value()) return true;
    quatY = quatY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto quatZ_t = StringUtils::fromString<float>(start, end);
    if (!quatZ_t.has_value()) return true;
    quatZ = quatZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto quatS_t = StringUtils::fromString<float>(start, end);
    if (!quatS_t.has_value()) return true;
    quatS = quatS_t.value();

    return false;  // no errors occurred
}

bool operator==(const Quaternion& lhs, const Quaternion& rhs)
{
    return ((lhs.quatX == rhs.quatX) && (lhs.quatY == rhs.quatY) && (lhs.quatZ == rhs.quatZ) && (lhs.quatS == rhs.quatS));
}

bool QuatMagAccelRate::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 13 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto quatX_t = StringUtils::fromString<float>(start, end);
    if (!quatX_t.has_value()) return true;
    quatX = quatX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto quatY_t = StringUtils::fromString<float>(start, end);
    if (!quatY_t.has_value()) return true;
    quatY = quatY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto quatZ_t = StringUtils::fromString<float>(start, end);
    if (!quatZ_t.has_value()) return true;
    quatZ = quatZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto quatS_t = StringUtils::fromString<float>(start, end);
    if (!quatS_t.has_value()) return true;
    quatS = quatS_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magX_t = StringUtils::fromString<float>(start, end);
    if (!magX_t.has_value()) return true;
    magX = magX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magY_t = StringUtils::fromString<float>(start, end);
    if (!magY_t.has_value()) return true;
    magY = magY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magZ_t = StringUtils::fromString<float>(start, end);
    if (!magZ_t.has_value()) return true;
    magZ = magZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelX_t = StringUtils::fromString<float>(start, end);
    if (!accelX_t.has_value()) return true;
    accelX = accelX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelY_t = StringUtils::fromString<float>(start, end);
    if (!accelY_t.has_value()) return true;
    accelY = accelY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelZ_t = StringUtils::fromString<float>(start, end);
    if (!accelZ_t.has_value()) return true;
    accelZ = accelZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroX_t = StringUtils::fromString<float>(start, end);
    if (!gyroX_t.has_value()) return true;
    gyroX = gyroX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroY_t = StringUtils::fromString<float>(start, end);
    if (!gyroY_t.has_value()) return true;
    gyroY = gyroY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroZ_t = StringUtils::fromString<float>(start, end);
    if (!gyroZ_t.has_value()) return true;
    gyroZ = gyroZ_t.value();

    return false;  // no errors occurred
}

bool operator==(const QuatMagAccelRate& lhs, const QuatMagAccelRate& rhs)
{
    return ((lhs.quatX == rhs.quatX) && (lhs.quatY == rhs.quatY) && (lhs.quatZ == rhs.quatZ) && (lhs.quatS == rhs.quatS) && (lhs.magX == rhs.magX) &&
            (lhs.magY == rhs.magY) && (lhs.magZ == rhs.magZ) && (lhs.accelX == rhs.accelX) && (lhs.accelY == rhs.accelY) && (lhs.accelZ == rhs.accelZ) &&
            (lhs.gyroX == rhs.gyroX) && (lhs.gyroY == rhs.gyroY) && (lhs.gyroZ == rhs.gyroZ));
}

bool MagGravRefVec::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 6 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magRefN_t = StringUtils::fromString<float>(start, end);
    if (!magRefN_t.has_value()) return true;
    magRefN = magRefN_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magRefE_t = StringUtils::fromString<float>(start, end);
    if (!magRefE_t.has_value()) return true;
    magRefE = magRefE_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magRefD_t = StringUtils::fromString<float>(start, end);
    if (!magRefD_t.has_value()) return true;
    magRefD = magRefD_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gravRefN_t = StringUtils::fromString<float>(start, end);
    if (!gravRefN_t.has_value()) return true;
    gravRefN = gravRefN_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gravRefE_t = StringUtils::fromString<float>(start, end);
    if (!gravRefE_t.has_value()) return true;
    gravRefE = gravRefE_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gravRefD_t = StringUtils::fromString<float>(start, end);
    if (!gravRefD_t.has_value()) return true;
    gravRefD = gravRefD_t.value();

    return false;  // no errors occurred
}

AsciiMessage MagGravRefVec::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%f,%f,%f,%f,%f,%f", magRefN, magRefE, magRefD, gravRefN, gravRefE, gravRefD);
    return result;
}

bool operator==(const MagGravRefVec& lhs, const MagGravRefVec& rhs)
{
    return ((lhs.magRefN == rhs.magRefN) && (lhs.magRefE == rhs.magRefE) && (lhs.magRefD == rhs.magRefD) && (lhs.gravRefN == rhs.gravRefN) &&
            (lhs.gravRefE == rhs.gravRefE) && (lhs.gravRefD == rhs.gravRefD));
}

bool YprMagAccelAngularRates::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 12 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto yaw_t = StringUtils::fromString<float>(start, end);
    if (!yaw_t.has_value()) return true;
    yaw = yaw_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto pitch_t = StringUtils::fromString<float>(start, end);
    if (!pitch_t.has_value()) return true;
    pitch = pitch_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto roll_t = StringUtils::fromString<float>(start, end);
    if (!roll_t.has_value()) return true;
    roll = roll_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magX_t = StringUtils::fromString<float>(start, end);
    if (!magX_t.has_value()) return true;
    magX = magX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magY_t = StringUtils::fromString<float>(start, end);
    if (!magY_t.has_value()) return true;
    magY = magY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magZ_t = StringUtils::fromString<float>(start, end);
    if (!magZ_t.has_value()) return true;
    magZ = magZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelX_t = StringUtils::fromString<float>(start, end);
    if (!accelX_t.has_value()) return true;
    accelX = accelX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelY_t = StringUtils::fromString<float>(start, end);
    if (!accelY_t.has_value()) return true;
    accelY = accelY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelZ_t = StringUtils::fromString<float>(start, end);
    if (!accelZ_t.has_value()) return true;
    accelZ = accelZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroX_t = StringUtils::fromString<float>(start, end);
    if (!gyroX_t.has_value()) return true;
    gyroX = gyroX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroY_t = StringUtils::fromString<float>(start, end);
    if (!gyroY_t.has_value()) return true;
    gyroY = gyroY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroZ_t = StringUtils::fromString<float>(start, end);
    if (!gyroZ_t.has_value()) return true;
    gyroZ = gyroZ_t.value();

    return false;  // no errors occurred
}

bool operator==(const YprMagAccelAngularRates& lhs, const YprMagAccelAngularRates& rhs)
{
    return ((lhs.yaw == rhs.yaw) && (lhs.pitch == rhs.pitch) && (lhs.roll == rhs.roll) && (lhs.magX == rhs.magX) && (lhs.magY == rhs.magY) &&
            (lhs.magZ == rhs.magZ) && (lhs.accelX == rhs.accelX) && (lhs.accelY == rhs.accelY) && (lhs.accelZ == rhs.accelZ) && (lhs.gyroX == rhs.gyroX) &&
            (lhs.gyroY == rhs.gyroY) && (lhs.gyroZ == rhs.gyroZ));
}

bool VpeBasicControl::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 4 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto resv_t = StringUtils::fromString<uint8_t>(start, end);
    if (!resv_t.has_value()) return true;
    resv = resv_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto headingMode_t = StringUtils::fromString<uint8_t>(start, end);
    if (!headingMode_t.has_value()) return true;
    headingMode = static_cast<VpeBasicControl::HeadingMode>(headingMode_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto filteringMode_t = StringUtils::fromString<uint8_t>(start, end);
    if (!filteringMode_t.has_value()) return true;
    filteringMode = static_cast<VpeBasicControl::FilteringMode>(filteringMode_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto tuningMode_t = StringUtils::fromString<uint8_t>(start, end);
    if (!tuningMode_t.has_value()) return true;
    tuningMode = static_cast<VpeBasicControl::TuningMode>(tuningMode_t.value());

    return false;  // no errors occurred
}

AsciiMessage VpeBasicControl::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%u,%u,%u,%u", resv, static_cast<uint8_t>(headingMode), static_cast<uint8_t>(filteringMode),
                  static_cast<uint8_t>(tuningMode));
    return result;
}

bool operator==(const VpeBasicControl& lhs, const VpeBasicControl& rhs)
{
    return ((lhs.resv == rhs.resv) && (lhs.headingMode == rhs.headingMode) && (lhs.filteringMode == rhs.filteringMode) && (lhs.tuningMode == rhs.tuningMode));
}

bool VpeMagBasicTuning::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 9 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto baseTuningX_t = StringUtils::fromString<float>(start, end);
    if (!baseTuningX_t.has_value()) return true;
    baseTuningX = baseTuningX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto baseTuningY_t = StringUtils::fromString<float>(start, end);
    if (!baseTuningY_t.has_value()) return true;
    baseTuningY = baseTuningY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto baseTuningZ_t = StringUtils::fromString<float>(start, end);
    if (!baseTuningZ_t.has_value()) return true;
    baseTuningZ = baseTuningZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto adaptiveTuningX_t = StringUtils::fromString<float>(start, end);
    if (!adaptiveTuningX_t.has_value()) return true;
    adaptiveTuningX = adaptiveTuningX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto adaptiveTuningY_t = StringUtils::fromString<float>(start, end);
    if (!adaptiveTuningY_t.has_value()) return true;
    adaptiveTuningY = adaptiveTuningY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto adaptiveTuningZ_t = StringUtils::fromString<float>(start, end);
    if (!adaptiveTuningZ_t.has_value()) return true;
    adaptiveTuningZ = adaptiveTuningZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto adaptiveFilteringX_t = StringUtils::fromString<float>(start, end);
    if (!adaptiveFilteringX_t.has_value()) return true;
    adaptiveFilteringX = adaptiveFilteringX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto adaptiveFilteringY_t = StringUtils::fromString<float>(start, end);
    if (!adaptiveFilteringY_t.has_value()) return true;
    adaptiveFilteringY = adaptiveFilteringY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto adaptiveFilteringZ_t = StringUtils::fromString<float>(start, end);
    if (!adaptiveFilteringZ_t.has_value()) return true;
    adaptiveFilteringZ = adaptiveFilteringZ_t.value();

    return false;  // no errors occurred
}

AsciiMessage VpeMagBasicTuning::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%f,%f,%f,%f,%f,%f,%f,%f,%f", baseTuningX, baseTuningY, baseTuningZ, adaptiveTuningX, adaptiveTuningY,
                  adaptiveTuningZ, adaptiveFilteringX, adaptiveFilteringY, adaptiveFilteringZ);
    return result;
}

bool operator==(const VpeMagBasicTuning& lhs, const VpeMagBasicTuning& rhs)
{
    return ((lhs.baseTuningX == rhs.baseTuningX) && (lhs.baseTuningY == rhs.baseTuningY) && (lhs.baseTuningZ == rhs.baseTuningZ) &&
            (lhs.adaptiveTuningX == rhs.adaptiveTuningX) && (lhs.adaptiveTuningY == rhs.adaptiveTuningY) && (lhs.adaptiveTuningZ == rhs.adaptiveTuningZ) &&
            (lhs.adaptiveFilteringX == rhs.adaptiveFilteringX) && (lhs.adaptiveFilteringY == rhs.adaptiveFilteringY) &&
            (lhs.adaptiveFilteringZ == rhs.adaptiveFilteringZ));
}

bool VpeAccelBasicTuning::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 9 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto baseTuningX_t = StringUtils::fromString<float>(start, end);
    if (!baseTuningX_t.has_value()) return true;
    baseTuningX = baseTuningX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto baseTuningY_t = StringUtils::fromString<float>(start, end);
    if (!baseTuningY_t.has_value()) return true;
    baseTuningY = baseTuningY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto baseTuningZ_t = StringUtils::fromString<float>(start, end);
    if (!baseTuningZ_t.has_value()) return true;
    baseTuningZ = baseTuningZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto adaptiveTuningX_t = StringUtils::fromString<float>(start, end);
    if (!adaptiveTuningX_t.has_value()) return true;
    adaptiveTuningX = adaptiveTuningX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto adaptiveTuningY_t = StringUtils::fromString<float>(start, end);
    if (!adaptiveTuningY_t.has_value()) return true;
    adaptiveTuningY = adaptiveTuningY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto adaptiveTuningZ_t = StringUtils::fromString<float>(start, end);
    if (!adaptiveTuningZ_t.has_value()) return true;
    adaptiveTuningZ = adaptiveTuningZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto adaptiveFilteringX_t = StringUtils::fromString<float>(start, end);
    if (!adaptiveFilteringX_t.has_value()) return true;
    adaptiveFilteringX = adaptiveFilteringX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto adaptiveFilteringY_t = StringUtils::fromString<float>(start, end);
    if (!adaptiveFilteringY_t.has_value()) return true;
    adaptiveFilteringY = adaptiveFilteringY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto adaptiveFilteringZ_t = StringUtils::fromString<float>(start, end);
    if (!adaptiveFilteringZ_t.has_value()) return true;
    adaptiveFilteringZ = adaptiveFilteringZ_t.value();

    return false;  // no errors occurred
}

AsciiMessage VpeAccelBasicTuning::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%f,%f,%f,%f,%f,%f,%f,%f,%f", baseTuningX, baseTuningY, baseTuningZ, adaptiveTuningX, adaptiveTuningY,
                  adaptiveTuningZ, adaptiveFilteringX, adaptiveFilteringY, adaptiveFilteringZ);
    return result;
}

bool operator==(const VpeAccelBasicTuning& lhs, const VpeAccelBasicTuning& rhs)
{
    return ((lhs.baseTuningX == rhs.baseTuningX) && (lhs.baseTuningY == rhs.baseTuningY) && (lhs.baseTuningZ == rhs.baseTuningZ) &&
            (lhs.adaptiveTuningX == rhs.adaptiveTuningX) && (lhs.adaptiveTuningY == rhs.adaptiveTuningY) && (lhs.adaptiveTuningZ == rhs.adaptiveTuningZ) &&
            (lhs.adaptiveFilteringX == rhs.adaptiveFilteringX) && (lhs.adaptiveFilteringY == rhs.adaptiveFilteringY) &&
            (lhs.adaptiveFilteringZ == rhs.adaptiveFilteringZ));
}

bool YprLinearBodyAccelAngularRates::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 9 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto yaw_t = StringUtils::fromString<float>(start, end);
    if (!yaw_t.has_value()) return true;
    yaw = yaw_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto pitch_t = StringUtils::fromString<float>(start, end);
    if (!pitch_t.has_value()) return true;
    pitch = pitch_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto roll_t = StringUtils::fromString<float>(start, end);
    if (!roll_t.has_value()) return true;
    roll = roll_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto linAccelX_t = StringUtils::fromString<float>(start, end);
    if (!linAccelX_t.has_value()) return true;
    linAccelX = linAccelX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto linAccelY_t = StringUtils::fromString<float>(start, end);
    if (!linAccelY_t.has_value()) return true;
    linAccelY = linAccelY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto linAccelZ_t = StringUtils::fromString<float>(start, end);
    if (!linAccelZ_t.has_value()) return true;
    linAccelZ = linAccelZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroX_t = StringUtils::fromString<float>(start, end);
    if (!gyroX_t.has_value()) return true;
    gyroX = gyroX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroY_t = StringUtils::fromString<float>(start, end);
    if (!gyroY_t.has_value()) return true;
    gyroY = gyroY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroZ_t = StringUtils::fromString<float>(start, end);
    if (!gyroZ_t.has_value()) return true;
    gyroZ = gyroZ_t.value();

    return false;  // no errors occurred
}

bool operator==(const YprLinearBodyAccelAngularRates& lhs, const YprLinearBodyAccelAngularRates& rhs)
{
    return ((lhs.yaw == rhs.yaw) && (lhs.pitch == rhs.pitch) && (lhs.roll == rhs.roll) && (lhs.linAccelX == rhs.linAccelX) &&
            (lhs.linAccelY == rhs.linAccelY) && (lhs.linAccelZ == rhs.linAccelZ) && (lhs.gyroX == rhs.gyroX) && (lhs.gyroY == rhs.gyroY) &&
            (lhs.gyroZ == rhs.gyroZ));
}

bool YprLinearInertialAccelAngularRates::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 9 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto yaw_t = StringUtils::fromString<float>(start, end);
    if (!yaw_t.has_value()) return true;
    yaw = yaw_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto pitch_t = StringUtils::fromString<float>(start, end);
    if (!pitch_t.has_value()) return true;
    pitch = pitch_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto roll_t = StringUtils::fromString<float>(start, end);
    if (!roll_t.has_value()) return true;
    roll = roll_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto linAccelN_t = StringUtils::fromString<float>(start, end);
    if (!linAccelN_t.has_value()) return true;
    linAccelN = linAccelN_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto linAccelE_t = StringUtils::fromString<float>(start, end);
    if (!linAccelE_t.has_value()) return true;
    linAccelE = linAccelE_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto linAccelD_t = StringUtils::fromString<float>(start, end);
    if (!linAccelD_t.has_value()) return true;
    linAccelD = linAccelD_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroX_t = StringUtils::fromString<float>(start, end);
    if (!gyroX_t.has_value()) return true;
    gyroX = gyroX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroY_t = StringUtils::fromString<float>(start, end);
    if (!gyroY_t.has_value()) return true;
    gyroY = gyroY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroZ_t = StringUtils::fromString<float>(start, end);
    if (!gyroZ_t.has_value()) return true;
    gyroZ = gyroZ_t.value();

    return false;  // no errors occurred
}

bool operator==(const YprLinearInertialAccelAngularRates& lhs, const YprLinearInertialAccelAngularRates& rhs)
{
    return ((lhs.yaw == rhs.yaw) && (lhs.pitch == rhs.pitch) && (lhs.roll == rhs.roll) && (lhs.linAccelN == rhs.linAccelN) &&
            (lhs.linAccelE == rhs.linAccelE) && (lhs.linAccelD == rhs.linAccelD) && (lhs.gyroX == rhs.gyroX) && (lhs.gyroY == rhs.gyroY) &&
            (lhs.gyroZ == rhs.gyroZ));
}

}  // namespace Attitude

namespace GNSS
{
bool GnssBasicConfig::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 5 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto receiverEnable_t = StringUtils::fromString<uint8_t>(start, end);
    if (!receiverEnable_t.has_value()) return true;
    receiverEnable = static_cast<GnssBasicConfig::ReceiverEnable>(receiverEnable_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto ppsSource_t = StringUtils::fromString<uint8_t>(start, end);
    if (!ppsSource_t.has_value()) return true;
    ppsSource = static_cast<GnssBasicConfig::PpsSource>(ppsSource_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto rate_t = StringUtils::fromString<uint8_t>(start, end);
    if (!rate_t.has_value()) return true;
    rate = static_cast<GnssBasicConfig::Rate>(rate_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto resv4_t = StringUtils::fromString<uint8_t>(start, end);
    if (!resv4_t.has_value()) return true;
    resv4 = resv4_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto antPower_t = StringUtils::fromString<uint8_t>(start, end);
    if (!antPower_t.has_value()) return true;
    antPower = static_cast<GnssBasicConfig::AntPower>(antPower_t.value());

    return false;  // no errors occurred
}

AsciiMessage GnssBasicConfig::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%u,%u,%u,%u,%u", static_cast<uint8_t>(receiverEnable), static_cast<uint8_t>(ppsSource),
                  static_cast<uint8_t>(rate), resv4, static_cast<uint8_t>(antPower));
    return result;
}

bool operator==(const GnssBasicConfig& lhs, const GnssBasicConfig& rhs)
{
    return ((lhs.receiverEnable == rhs.receiverEnable) && (lhs.ppsSource == rhs.ppsSource) && (lhs.rate == rhs.rate) && (lhs.resv4 == rhs.resv4) &&
            (lhs.antPower == rhs.antPower));
}

bool GnssAOffset::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 3 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto positionX_t = StringUtils::fromString<float>(start, end);
    if (!positionX_t.has_value()) return true;
    positionX = positionX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto positionY_t = StringUtils::fromString<float>(start, end);
    if (!positionY_t.has_value()) return true;
    positionY = positionY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto positionZ_t = StringUtils::fromString<float>(start, end);
    if (!positionZ_t.has_value()) return true;
    positionZ = positionZ_t.value();

    return false;  // no errors occurred
}

AsciiMessage GnssAOffset::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%f,%f,%f", positionX, positionY, positionZ);
    return result;
}

bool operator==(const GnssAOffset& lhs, const GnssAOffset& rhs)
{
    return ((lhs.positionX == rhs.positionX) && (lhs.positionY == rhs.positionY) && (lhs.positionZ == rhs.positionZ));
}

bool GnssSolLla::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 15 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gpsTow_t = StringUtils::fromString<double>(start, end);
    if (!gpsTow_t.has_value()) return true;
    gpsTow = gpsTow_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gps1Week_t = StringUtils::fromString<uint16_t>(start, end);
    if (!gps1Week_t.has_value()) return true;
    gps1Week = gps1Week_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnss1Fix_t = StringUtils::fromString<uint8_t>(start, end);
    if (!gnss1Fix_t.has_value()) return true;
    gnss1Fix = static_cast<GnssSolLla::Gnss1Fix>(gnss1Fix_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnss1NumSats_t = StringUtils::fromString<uint8_t>(start, end);
    if (!gnss1NumSats_t.has_value()) return true;
    gnss1NumSats = gnss1NumSats_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto lat_t = StringUtils::fromString<double>(start, end);
    if (!lat_t.has_value()) return true;
    lat = lat_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto lon_t = StringUtils::fromString<double>(start, end);
    if (!lon_t.has_value()) return true;
    lon = lon_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto alt_t = StringUtils::fromString<double>(start, end);
    if (!alt_t.has_value()) return true;
    alt = alt_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velN_t = StringUtils::fromString<float>(start, end);
    if (!velN_t.has_value()) return true;
    velN = velN_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velE_t = StringUtils::fromString<float>(start, end);
    if (!velE_t.has_value()) return true;
    velE = velE_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velD_t = StringUtils::fromString<float>(start, end);
    if (!velD_t.has_value()) return true;
    velD = velD_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posUncertaintyN_t = StringUtils::fromString<float>(start, end);
    if (!posUncertaintyN_t.has_value()) return true;
    posUncertaintyN = posUncertaintyN_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posUncertaintyE_t = StringUtils::fromString<float>(start, end);
    if (!posUncertaintyE_t.has_value()) return true;
    posUncertaintyE = posUncertaintyE_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posUncertaintyD_t = StringUtils::fromString<float>(start, end);
    if (!posUncertaintyD_t.has_value()) return true;
    posUncertaintyD = posUncertaintyD_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnss1VelUncertainty_t = StringUtils::fromString<float>(start, end);
    if (!gnss1VelUncertainty_t.has_value()) return true;
    gnss1VelUncertainty = gnss1VelUncertainty_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnss1TimeUncertainty_t = StringUtils::fromString<float>(start, end);
    if (!gnss1TimeUncertainty_t.has_value()) return true;
    gnss1TimeUncertainty = gnss1TimeUncertainty_t.value();

    return false;  // no errors occurred
}

bool operator==(const GnssSolLla& lhs, const GnssSolLla& rhs)
{
    return ((lhs.gpsTow == rhs.gpsTow) && (lhs.gps1Week == rhs.gps1Week) && (lhs.gnss1Fix == rhs.gnss1Fix) && (lhs.gnss1NumSats == rhs.gnss1NumSats) &&
            (lhs.lat == rhs.lat) && (lhs.lon == rhs.lon) && (lhs.alt == rhs.alt) && (lhs.velN == rhs.velN) && (lhs.velE == rhs.velE) &&
            (lhs.velD == rhs.velD) && (lhs.posUncertaintyN == rhs.posUncertaintyN) && (lhs.posUncertaintyE == rhs.posUncertaintyE) &&
            (lhs.posUncertaintyD == rhs.posUncertaintyD) && (lhs.gnss1VelUncertainty == rhs.gnss1VelUncertainty) &&
            (lhs.gnss1TimeUncertainty == rhs.gnss1TimeUncertainty));
}

bool GnssSolEcef::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 15 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gpsTow_t = StringUtils::fromString<double>(start, end);
    if (!gpsTow_t.has_value()) return true;
    gpsTow = gpsTow_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gps1Week_t = StringUtils::fromString<uint16_t>(start, end);
    if (!gps1Week_t.has_value()) return true;
    gps1Week = gps1Week_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnss1Fix_t = StringUtils::fromString<uint8_t>(start, end);
    if (!gnss1Fix_t.has_value()) return true;
    gnss1Fix = static_cast<GnssSolEcef::Gnss1Fix>(gnss1Fix_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnss1NumSats_t = StringUtils::fromString<uint8_t>(start, end);
    if (!gnss1NumSats_t.has_value()) return true;
    gnss1NumSats = gnss1NumSats_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posX_t = StringUtils::fromString<double>(start, end);
    if (!posX_t.has_value()) return true;
    posX = posX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posY_t = StringUtils::fromString<double>(start, end);
    if (!posY_t.has_value()) return true;
    posY = posY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posZ_t = StringUtils::fromString<double>(start, end);
    if (!posZ_t.has_value()) return true;
    posZ = posZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velX_t = StringUtils::fromString<float>(start, end);
    if (!velX_t.has_value()) return true;
    velX = velX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velY_t = StringUtils::fromString<float>(start, end);
    if (!velY_t.has_value()) return true;
    velY = velY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velZ_t = StringUtils::fromString<float>(start, end);
    if (!velZ_t.has_value()) return true;
    velZ = velZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posUncertaintyX_t = StringUtils::fromString<float>(start, end);
    if (!posUncertaintyX_t.has_value()) return true;
    posUncertaintyX = posUncertaintyX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posUncertaintyY_t = StringUtils::fromString<float>(start, end);
    if (!posUncertaintyY_t.has_value()) return true;
    posUncertaintyY = posUncertaintyY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posUncertaintyZ_t = StringUtils::fromString<float>(start, end);
    if (!posUncertaintyZ_t.has_value()) return true;
    posUncertaintyZ = posUncertaintyZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnss1VelUncertainty_t = StringUtils::fromString<float>(start, end);
    if (!gnss1VelUncertainty_t.has_value()) return true;
    gnss1VelUncertainty = gnss1VelUncertainty_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnss1TimeUncertainty_t = StringUtils::fromString<float>(start, end);
    if (!gnss1TimeUncertainty_t.has_value()) return true;
    gnss1TimeUncertainty = gnss1TimeUncertainty_t.value();

    return false;  // no errors occurred
}

bool operator==(const GnssSolEcef& lhs, const GnssSolEcef& rhs)
{
    return ((lhs.gpsTow == rhs.gpsTow) && (lhs.gps1Week == rhs.gps1Week) && (lhs.gnss1Fix == rhs.gnss1Fix) && (lhs.gnss1NumSats == rhs.gnss1NumSats) &&
            (lhs.posX == rhs.posX) && (lhs.posY == rhs.posY) && (lhs.posZ == rhs.posZ) && (lhs.velX == rhs.velX) && (lhs.velY == rhs.velY) &&
            (lhs.velZ == rhs.velZ) && (lhs.posUncertaintyX == rhs.posUncertaintyX) && (lhs.posUncertaintyY == rhs.posUncertaintyY) &&
            (lhs.posUncertaintyZ == rhs.posUncertaintyZ) && (lhs.gnss1VelUncertainty == rhs.gnss1VelUncertainty) &&
            (lhs.gnss1TimeUncertainty == rhs.gnss1TimeUncertainty));
}

bool GnssSystemConfig::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 9 + 1)
    {
        // This register may not have an optional field
        if (tokens.size() != 8 + 1) { return true; }
    }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto systems_t = StringUtils::fromStringHex<uint16_t>(start, end);
    if (!systems_t.has_value()) return true;
    systems._value = systems_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto minCno_t = StringUtils::fromString<uint8_t>(start, end);
    if (!minCno_t.has_value()) return true;
    minCno = minCno_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto minElev_t = StringUtils::fromString<uint8_t>(start, end);
    if (!minElev_t.has_value()) return true;
    minElev = minElev_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto maxSats_t = StringUtils::fromString<uint8_t>(start, end);
    if (!maxSats_t.has_value()) return true;
    maxSats = maxSats_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto sbasMode_t = StringUtils::fromStringHex<uint8_t>(start, end);
    if (!sbasMode_t.has_value()) return true;
    sbasMode._value = sbasMode_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto sbasSelect1_t = StringUtils::fromStringHex<uint16_t>(start, end);
    if (!sbasSelect1_t.has_value()) return true;
    sbasSelect1._value = sbasSelect1_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto sbasSelect2_t = StringUtils::fromStringHex<uint16_t>(start, end);
    if (!sbasSelect2_t.has_value()) return true;
    sbasSelect2._value = sbasSelect2_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto sbasSelect3_t = StringUtils::fromStringHex<uint16_t>(start, end);
    if (!sbasSelect3_t.has_value()) return true;
    sbasSelect3._value = sbasSelect3_t.value();

    if (tokens.size() == (9 + 1))
    {
        // This is an optional parameter.
        start = tokens[index++] + 1;
        end = tokens[index];
        const auto receiverSelect_t = StringUtils::fromString<uint8_t>(start, end);
        if (!receiverSelect_t.has_value()) return true;
        receiverSelect = static_cast<GnssSystemConfig::ReceiverSelect>(receiverSelect_t.value());
    }

    return false;  // no errors occurred
}

Command GnssSystemConfig::toReadCommand()
{
    AsciiMessage commandString;
    if (receiverSelect == GNSS::GnssSystemConfig::ReceiverSelect::GnssAB) { std::snprintf(commandString.begin(), commandString.capacity(), "RRG,%02d", _id); }
    else { std::snprintf(commandString.begin(), commandString.capacity(), "RRG,%02d,%1d", _id, static_cast<uint8_t>(receiverSelect)); }
    Command readCommand(commandString, 6);
    return readCommand;
}

AsciiMessage GnssSystemConfig::toString() const
{
    AsciiMessage result = "";
    if (receiverSelect == GNSS::GnssSystemConfig::ReceiverSelect::GnssAB)
    {
        std::snprintf(result.begin(), result.capacity(), "%04X,%d,%d,%d,%02X,%04X,%04X,%04X", systems._value, minCno, minElev, maxSats, sbasMode._value,
                      sbasSelect1._value, sbasSelect2._value, sbasSelect3._value);
    }
    else
    {
        std::snprintf(result.begin(), result.capacity(), "%04X,%d,%d,%d,%02X,%04X,%04X,%04X,%d", systems._value, minCno, minElev, maxSats, sbasMode._value,
                      sbasSelect1._value, sbasSelect2._value, sbasSelect3._value, static_cast<uint8_t>(receiverSelect));
    }
    return result;
}

bool operator==(const GnssSystemConfig& lhs, const GnssSystemConfig& rhs)
{
    return ((lhs.systems._value == rhs.systems._value) && (lhs.minCno == rhs.minCno) && (lhs.minElev == rhs.minElev) && (lhs.maxSats == rhs.maxSats) &&
            (lhs.sbasMode._value == rhs.sbasMode._value) && (lhs.sbasSelect1._value == rhs.sbasSelect1._value) &&
            (lhs.sbasSelect2._value == rhs.sbasSelect2._value) && (lhs.sbasSelect3._value == rhs.sbasSelect3._value) &&
            (lhs.receiverSelect == rhs.receiverSelect));
}

bool GnssSyncConfig::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 7 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnssSyncEnable_t = StringUtils::fromString<uint8_t>(start, end);
    if (!gnssSyncEnable_t.has_value()) return true;
    gnssSyncEnable = static_cast<GnssSyncConfig::GnssSyncEnable>(gnssSyncEnable_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto polarity_t = StringUtils::fromString<uint8_t>(start, end);
    if (!polarity_t.has_value()) return true;
    polarity = static_cast<GnssSyncConfig::Polarity>(polarity_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto specType_t = StringUtils::fromString<uint8_t>(start, end);
    if (!specType_t.has_value()) return true;
    specType = static_cast<GnssSyncConfig::SpecType>(specType_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto resv_t = StringUtils::fromString<uint8_t>(start, end);
    if (!resv_t.has_value()) return true;
    resv = resv_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto period_t = StringUtils::fromString<uint32_t>(start, end);
    if (!period_t.has_value()) return true;
    period = period_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto pulseWidth_t = StringUtils::fromString<uint32_t>(start, end);
    if (!pulseWidth_t.has_value()) return true;
    pulseWidth = pulseWidth_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto offset_t = StringUtils::fromString<int32_t>(start, end);
    if (!offset_t.has_value()) return true;
    offset = offset_t.value();

    return false;  // no errors occurred
}

AsciiMessage GnssSyncConfig::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%u,%u,%u,%u,%u,%u,%d", static_cast<uint8_t>(gnssSyncEnable), static_cast<uint8_t>(polarity),
                  static_cast<uint8_t>(specType), resv, period, pulseWidth, offset);
    return result;
}

bool operator==(const GnssSyncConfig& lhs, const GnssSyncConfig& rhs)
{
    return ((lhs.gnssSyncEnable == rhs.gnssSyncEnable) && (lhs.polarity == rhs.polarity) && (lhs.specType == rhs.specType) && (lhs.resv == rhs.resv) &&
            (lhs.period == rhs.period) && (lhs.pulseWidth == rhs.pulseWidth) && (lhs.offset == rhs.offset));
}

bool Gnss2SolLla::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 15 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gpsTow_t = StringUtils::fromString<double>(start, end);
    if (!gpsTow_t.has_value()) return true;
    gpsTow = gpsTow_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gps2Week_t = StringUtils::fromString<uint16_t>(start, end);
    if (!gps2Week_t.has_value()) return true;
    gps2Week = gps2Week_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnss2Fix_t = StringUtils::fromString<uint8_t>(start, end);
    if (!gnss2Fix_t.has_value()) return true;
    gnss2Fix = static_cast<Gnss2SolLla::Gnss2Fix>(gnss2Fix_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnss2NumSats_t = StringUtils::fromString<uint8_t>(start, end);
    if (!gnss2NumSats_t.has_value()) return true;
    gnss2NumSats = gnss2NumSats_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto lat_t = StringUtils::fromString<double>(start, end);
    if (!lat_t.has_value()) return true;
    lat = lat_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto lon_t = StringUtils::fromString<double>(start, end);
    if (!lon_t.has_value()) return true;
    lon = lon_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto alt_t = StringUtils::fromString<double>(start, end);
    if (!alt_t.has_value()) return true;
    alt = alt_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velN_t = StringUtils::fromString<float>(start, end);
    if (!velN_t.has_value()) return true;
    velN = velN_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velE_t = StringUtils::fromString<float>(start, end);
    if (!velE_t.has_value()) return true;
    velE = velE_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velD_t = StringUtils::fromString<float>(start, end);
    if (!velD_t.has_value()) return true;
    velD = velD_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posUncertaintyN_t = StringUtils::fromString<float>(start, end);
    if (!posUncertaintyN_t.has_value()) return true;
    posUncertaintyN = posUncertaintyN_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posUncertaintyE_t = StringUtils::fromString<float>(start, end);
    if (!posUncertaintyE_t.has_value()) return true;
    posUncertaintyE = posUncertaintyE_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posUncertaintyD_t = StringUtils::fromString<float>(start, end);
    if (!posUncertaintyD_t.has_value()) return true;
    posUncertaintyD = posUncertaintyD_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnss2VelUncertainty_t = StringUtils::fromString<float>(start, end);
    if (!gnss2VelUncertainty_t.has_value()) return true;
    gnss2VelUncertainty = gnss2VelUncertainty_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnss2TimeUncertainty_t = StringUtils::fromString<float>(start, end);
    if (!gnss2TimeUncertainty_t.has_value()) return true;
    gnss2TimeUncertainty = gnss2TimeUncertainty_t.value();

    return false;  // no errors occurred
}

bool operator==(const Gnss2SolLla& lhs, const Gnss2SolLla& rhs)
{
    return ((lhs.gpsTow == rhs.gpsTow) && (lhs.gps2Week == rhs.gps2Week) && (lhs.gnss2Fix == rhs.gnss2Fix) && (lhs.gnss2NumSats == rhs.gnss2NumSats) &&
            (lhs.lat == rhs.lat) && (lhs.lon == rhs.lon) && (lhs.alt == rhs.alt) && (lhs.velN == rhs.velN) && (lhs.velE == rhs.velE) &&
            (lhs.velD == rhs.velD) && (lhs.posUncertaintyN == rhs.posUncertaintyN) && (lhs.posUncertaintyE == rhs.posUncertaintyE) &&
            (lhs.posUncertaintyD == rhs.posUncertaintyD) && (lhs.gnss2VelUncertainty == rhs.gnss2VelUncertainty) &&
            (lhs.gnss2TimeUncertainty == rhs.gnss2TimeUncertainty));
}

bool Gnss2SolEcef::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 15 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gpsTow_t = StringUtils::fromString<double>(start, end);
    if (!gpsTow_t.has_value()) return true;
    gpsTow = gpsTow_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gps2Week_t = StringUtils::fromString<uint16_t>(start, end);
    if (!gps2Week_t.has_value()) return true;
    gps2Week = gps2Week_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnss2Fix_t = StringUtils::fromString<uint8_t>(start, end);
    if (!gnss2Fix_t.has_value()) return true;
    gnss2Fix = static_cast<Gnss2SolEcef::Gnss2Fix>(gnss2Fix_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnss2NumSats_t = StringUtils::fromString<uint8_t>(start, end);
    if (!gnss2NumSats_t.has_value()) return true;
    gnss2NumSats = gnss2NumSats_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posX_t = StringUtils::fromString<double>(start, end);
    if (!posX_t.has_value()) return true;
    posX = posX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posY_t = StringUtils::fromString<double>(start, end);
    if (!posY_t.has_value()) return true;
    posY = posY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posZ_t = StringUtils::fromString<double>(start, end);
    if (!posZ_t.has_value()) return true;
    posZ = posZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velX_t = StringUtils::fromString<float>(start, end);
    if (!velX_t.has_value()) return true;
    velX = velX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velY_t = StringUtils::fromString<float>(start, end);
    if (!velY_t.has_value()) return true;
    velY = velY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velZ_t = StringUtils::fromString<float>(start, end);
    if (!velZ_t.has_value()) return true;
    velZ = velZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posUncertaintyX_t = StringUtils::fromString<float>(start, end);
    if (!posUncertaintyX_t.has_value()) return true;
    posUncertaintyX = posUncertaintyX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posUncertaintyY_t = StringUtils::fromString<float>(start, end);
    if (!posUncertaintyY_t.has_value()) return true;
    posUncertaintyY = posUncertaintyY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posUncertaintyZ_t = StringUtils::fromString<float>(start, end);
    if (!posUncertaintyZ_t.has_value()) return true;
    posUncertaintyZ = posUncertaintyZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnss2VelUncertainty_t = StringUtils::fromString<float>(start, end);
    if (!gnss2VelUncertainty_t.has_value()) return true;
    gnss2VelUncertainty = gnss2VelUncertainty_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnss2TimeUncertainty_t = StringUtils::fromString<float>(start, end);
    if (!gnss2TimeUncertainty_t.has_value()) return true;
    gnss2TimeUncertainty = gnss2TimeUncertainty_t.value();

    return false;  // no errors occurred
}

bool operator==(const Gnss2SolEcef& lhs, const Gnss2SolEcef& rhs)
{
    return ((lhs.gpsTow == rhs.gpsTow) && (lhs.gps2Week == rhs.gps2Week) && (lhs.gnss2Fix == rhs.gnss2Fix) && (lhs.gnss2NumSats == rhs.gnss2NumSats) &&
            (lhs.posX == rhs.posX) && (lhs.posY == rhs.posY) && (lhs.posZ == rhs.posZ) && (lhs.velX == rhs.velX) && (lhs.velY == rhs.velY) &&
            (lhs.velZ == rhs.velZ) && (lhs.posUncertaintyX == rhs.posUncertaintyX) && (lhs.posUncertaintyY == rhs.posUncertaintyY) &&
            (lhs.posUncertaintyZ == rhs.posUncertaintyZ) && (lhs.gnss2VelUncertainty == rhs.gnss2VelUncertainty) &&
            (lhs.gnss2TimeUncertainty == rhs.gnss2TimeUncertainty));
}

bool ExtGnssOffset::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 3 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto positionX_t = StringUtils::fromString<float>(start, end);
    if (!positionX_t.has_value()) return true;
    positionX = positionX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto positionY_t = StringUtils::fromString<float>(start, end);
    if (!positionY_t.has_value()) return true;
    positionY = positionY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto positionZ_t = StringUtils::fromString<float>(start, end);
    if (!positionZ_t.has_value()) return true;
    positionZ = positionZ_t.value();

    return false;  // no errors occurred
}

AsciiMessage ExtGnssOffset::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%f,%f,%f", positionX, positionY, positionZ);
    return result;
}

bool operator==(const ExtGnssOffset& lhs, const ExtGnssOffset& rhs)
{
    return ((lhs.positionX == rhs.positionX) && (lhs.positionY == rhs.positionY) && (lhs.positionZ == rhs.positionZ));
}

}  // namespace GNSS

namespace GNSSCompass
{
bool GnssCompassSignalHealthStatus::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 8 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto numSatsPvtA_t = StringUtils::fromString<float>(start, end);
    if (!numSatsPvtA_t.has_value()) return true;
    numSatsPvtA = numSatsPvtA_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto numSatsRtkA_t = StringUtils::fromString<float>(start, end);
    if (!numSatsRtkA_t.has_value()) return true;
    numSatsRtkA = numSatsRtkA_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto highestCn0A_t = StringUtils::fromString<float>(start, end);
    if (!highestCn0A_t.has_value()) return true;
    highestCn0A = highestCn0A_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto numSatsPvtB_t = StringUtils::fromString<float>(start, end);
    if (!numSatsPvtB_t.has_value()) return true;
    numSatsPvtB = numSatsPvtB_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto numSatsRtkB_t = StringUtils::fromString<float>(start, end);
    if (!numSatsRtkB_t.has_value()) return true;
    numSatsRtkB = numSatsRtkB_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto highestCn0B_t = StringUtils::fromString<float>(start, end);
    if (!highestCn0B_t.has_value()) return true;
    highestCn0B = highestCn0B_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto numComSatsPVT_t = StringUtils::fromString<float>(start, end);
    if (!numComSatsPVT_t.has_value()) return true;
    numComSatsPVT = numComSatsPVT_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto numComSatsRTK_t = StringUtils::fromString<float>(start, end);
    if (!numComSatsRTK_t.has_value()) return true;
    numComSatsRTK = numComSatsRTK_t.value();

    return false;  // no errors occurred
}

bool operator==(const GnssCompassSignalHealthStatus& lhs, const GnssCompassSignalHealthStatus& rhs)
{
    return ((lhs.numSatsPvtA == rhs.numSatsPvtA) && (lhs.numSatsRtkA == rhs.numSatsRtkA) && (lhs.highestCn0A == rhs.highestCn0A) &&
            (lhs.numSatsPvtB == rhs.numSatsPvtB) && (lhs.numSatsRtkB == rhs.numSatsRtkB) && (lhs.highestCn0B == rhs.highestCn0B) &&
            (lhs.numComSatsPVT == rhs.numComSatsPVT) && (lhs.numComSatsRTK == rhs.numComSatsRTK));
}

bool GnssCompassBaseline::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 6 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto positionX_t = StringUtils::fromString<float>(start, end);
    if (!positionX_t.has_value()) return true;
    positionX = positionX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto positionY_t = StringUtils::fromString<float>(start, end);
    if (!positionY_t.has_value()) return true;
    positionY = positionY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto positionZ_t = StringUtils::fromString<float>(start, end);
    if (!positionZ_t.has_value()) return true;
    positionZ = positionZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto uncertaintyX_t = StringUtils::fromString<float>(start, end);
    if (!uncertaintyX_t.has_value()) return true;
    uncertaintyX = uncertaintyX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto uncertaintyY_t = StringUtils::fromString<float>(start, end);
    if (!uncertaintyY_t.has_value()) return true;
    uncertaintyY = uncertaintyY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto uncertaintyZ_t = StringUtils::fromString<float>(start, end);
    if (!uncertaintyZ_t.has_value()) return true;
    uncertaintyZ = uncertaintyZ_t.value();

    return false;  // no errors occurred
}

AsciiMessage GnssCompassBaseline::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%f,%f,%f,%f,%f,%f", positionX, positionY, positionZ, uncertaintyX, uncertaintyY, uncertaintyZ);
    return result;
}

bool operator==(const GnssCompassBaseline& lhs, const GnssCompassBaseline& rhs)
{
    return ((lhs.positionX == rhs.positionX) && (lhs.positionY == rhs.positionY) && (lhs.positionZ == rhs.positionZ) &&
            (lhs.uncertaintyX == rhs.uncertaintyX) && (lhs.uncertaintyY == rhs.uncertaintyY) && (lhs.uncertaintyZ == rhs.uncertaintyZ));
}

bool GnssCompassEstBaseline::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 9 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto estBaselineComplete_t = StringUtils::fromString<uint8_t>(start, end);
    if (!estBaselineComplete_t.has_value()) return true;
    estBaselineComplete = estBaselineComplete_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto resv_t = StringUtils::fromString<uint8_t>(start, end);
    if (!resv_t.has_value()) return true;
    resv = resv_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto numMeas_t = StringUtils::fromString<uint16_t>(start, end);
    if (!numMeas_t.has_value()) return true;
    numMeas = numMeas_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto positionX_t = StringUtils::fromString<float>(start, end);
    if (!positionX_t.has_value()) return true;
    positionX = positionX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto positionY_t = StringUtils::fromString<float>(start, end);
    if (!positionY_t.has_value()) return true;
    positionY = positionY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto positionZ_t = StringUtils::fromString<float>(start, end);
    if (!positionZ_t.has_value()) return true;
    positionZ = positionZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto uncertaintyX_t = StringUtils::fromString<float>(start, end);
    if (!uncertaintyX_t.has_value()) return true;
    uncertaintyX = uncertaintyX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto uncertaintyY_t = StringUtils::fromString<float>(start, end);
    if (!uncertaintyY_t.has_value()) return true;
    uncertaintyY = uncertaintyY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto uncertaintyZ_t = StringUtils::fromString<float>(start, end);
    if (!uncertaintyZ_t.has_value()) return true;
    uncertaintyZ = uncertaintyZ_t.value();

    return false;  // no errors occurred
}

bool operator==(const GnssCompassEstBaseline& lhs, const GnssCompassEstBaseline& rhs)
{
    return ((lhs.estBaselineComplete == rhs.estBaselineComplete) && (lhs.resv == rhs.resv) && (lhs.numMeas == rhs.numMeas) &&
            (lhs.positionX == rhs.positionX) && (lhs.positionY == rhs.positionY) && (lhs.positionZ == rhs.positionZ) &&
            (lhs.uncertaintyX == rhs.uncertaintyX) && (lhs.uncertaintyY == rhs.uncertaintyY) && (lhs.uncertaintyZ == rhs.uncertaintyZ));
}

bool GnssCompassStartupStatus::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 2 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto percentComplete_t = StringUtils::fromString<uint8_t>(start, end);
    if (!percentComplete_t.has_value()) return true;
    percentComplete = percentComplete_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto currentHeading_t = StringUtils::fromString<float>(start, end);
    if (!currentHeading_t.has_value()) return true;
    currentHeading = currentHeading_t.value();

    return false;  // no errors occurred
}

bool operator==(const GnssCompassStartupStatus& lhs, const GnssCompassStartupStatus& rhs)
{
    return ((lhs.percentComplete == rhs.percentComplete) && (lhs.currentHeading == rhs.currentHeading));
}

}  // namespace GNSSCompass

namespace HardSoftIronEstimator
{
bool RealTimeHsiControl::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 3 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto mode_t = StringUtils::fromString<uint8_t>(start, end);
    if (!mode_t.has_value()) return true;
    mode = static_cast<RealTimeHsiControl::Mode>(mode_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto applyCompensation_t = StringUtils::fromString<uint8_t>(start, end);
    if (!applyCompensation_t.has_value()) return true;
    applyCompensation = static_cast<RealTimeHsiControl::ApplyCompensation>(applyCompensation_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto convergeRate_t = StringUtils::fromString<uint8_t>(start, end);
    if (!convergeRate_t.has_value()) return true;
    convergeRate = convergeRate_t.value();

    return false;  // no errors occurred
}

AsciiMessage RealTimeHsiControl::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%u,%u,%u", static_cast<uint8_t>(mode), static_cast<uint8_t>(applyCompensation), convergeRate);
    return result;
}

bool operator==(const RealTimeHsiControl& lhs, const RealTimeHsiControl& rhs)
{
    return ((lhs.mode == rhs.mode) && (lhs.applyCompensation == rhs.applyCompensation) && (lhs.convergeRate == rhs.convergeRate));
}

bool EstMagCal::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 12 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magGain00_t = StringUtils::fromString<float>(start, end);
    if (!magGain00_t.has_value()) return true;
    magGain00 = magGain00_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magGain01_t = StringUtils::fromString<float>(start, end);
    if (!magGain01_t.has_value()) return true;
    magGain01 = magGain01_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magGain02_t = StringUtils::fromString<float>(start, end);
    if (!magGain02_t.has_value()) return true;
    magGain02 = magGain02_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magGain10_t = StringUtils::fromString<float>(start, end);
    if (!magGain10_t.has_value()) return true;
    magGain10 = magGain10_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magGain11_t = StringUtils::fromString<float>(start, end);
    if (!magGain11_t.has_value()) return true;
    magGain11 = magGain11_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magGain12_t = StringUtils::fromString<float>(start, end);
    if (!magGain12_t.has_value()) return true;
    magGain12 = magGain12_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magGain20_t = StringUtils::fromString<float>(start, end);
    if (!magGain20_t.has_value()) return true;
    magGain20 = magGain20_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magGain21_t = StringUtils::fromString<float>(start, end);
    if (!magGain21_t.has_value()) return true;
    magGain21 = magGain21_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magGain22_t = StringUtils::fromString<float>(start, end);
    if (!magGain22_t.has_value()) return true;
    magGain22 = magGain22_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magBiasX_t = StringUtils::fromString<float>(start, end);
    if (!magBiasX_t.has_value()) return true;
    magBiasX = magBiasX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magBiasY_t = StringUtils::fromString<float>(start, end);
    if (!magBiasY_t.has_value()) return true;
    magBiasY = magBiasY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magBiasZ_t = StringUtils::fromString<float>(start, end);
    if (!magBiasZ_t.has_value()) return true;
    magBiasZ = magBiasZ_t.value();

    return false;  // no errors occurred
}

bool operator==(const EstMagCal& lhs, const EstMagCal& rhs)
{
    return ((lhs.magGain00 == rhs.magGain00) && (lhs.magGain01 == rhs.magGain01) && (lhs.magGain02 == rhs.magGain02) && (lhs.magGain10 == rhs.magGain10) &&
            (lhs.magGain11 == rhs.magGain11) && (lhs.magGain12 == rhs.magGain12) && (lhs.magGain20 == rhs.magGain20) && (lhs.magGain21 == rhs.magGain21) &&
            (lhs.magGain22 == rhs.magGain22) && (lhs.magBiasX == rhs.magBiasX) && (lhs.magBiasY == rhs.magBiasY) && (lhs.magBiasZ == rhs.magBiasZ));
}

}  // namespace HardSoftIronEstimator

namespace Heave
{
bool HeaveOutputs::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 3 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto heave_t = StringUtils::fromString<float>(start, end);
    if (!heave_t.has_value()) return true;
    heave = heave_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto heaveRate_t = StringUtils::fromString<float>(start, end);
    if (!heaveRate_t.has_value()) return true;
    heaveRate = heaveRate_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto delayedHeave_t = StringUtils::fromString<float>(start, end);
    if (!delayedHeave_t.has_value()) return true;
    delayedHeave = delayedHeave_t.value();

    return false;  // no errors occurred
}

bool operator==(const HeaveOutputs& lhs, const HeaveOutputs& rhs)
{
    return ((lhs.heave == rhs.heave) && (lhs.heaveRate == rhs.heaveRate) && (lhs.delayedHeave == rhs.delayedHeave));
}

bool HeaveBasicConfig::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 7 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto initialWavePeriod_t = StringUtils::fromString<float>(start, end);
    if (!initialWavePeriod_t.has_value()) return true;
    initialWavePeriod = initialWavePeriod_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto initialWaveAmplitude_t = StringUtils::fromString<float>(start, end);
    if (!initialWaveAmplitude_t.has_value()) return true;
    initialWaveAmplitude = initialWaveAmplitude_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto maxWavePeriod_t = StringUtils::fromString<float>(start, end);
    if (!maxWavePeriod_t.has_value()) return true;
    maxWavePeriod = maxWavePeriod_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto minWaveAmplitude_t = StringUtils::fromString<float>(start, end);
    if (!minWaveAmplitude_t.has_value()) return true;
    minWaveAmplitude = minWaveAmplitude_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto delayedHeaveCutoffFreq_t = StringUtils::fromString<float>(start, end);
    if (!delayedHeaveCutoffFreq_t.has_value()) return true;
    delayedHeaveCutoffFreq = delayedHeaveCutoffFreq_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto heaveCutoffFreq_t = StringUtils::fromString<float>(start, end);
    if (!heaveCutoffFreq_t.has_value()) return true;
    heaveCutoffFreq = heaveCutoffFreq_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto heaveRateCutoffFreq_t = StringUtils::fromString<float>(start, end);
    if (!heaveRateCutoffFreq_t.has_value()) return true;
    heaveRateCutoffFreq = heaveRateCutoffFreq_t.value();

    return false;  // no errors occurred
}

AsciiMessage HeaveBasicConfig::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%f,%f,%f,%f,%f,%f,%f", initialWavePeriod, initialWaveAmplitude, maxWavePeriod, minWaveAmplitude,
                  delayedHeaveCutoffFreq, heaveCutoffFreq, heaveRateCutoffFreq);
    return result;
}

bool operator==(const HeaveBasicConfig& lhs, const HeaveBasicConfig& rhs)
{
    return ((lhs.initialWavePeriod == rhs.initialWavePeriod) && (lhs.initialWaveAmplitude == rhs.initialWaveAmplitude) &&
            (lhs.maxWavePeriod == rhs.maxWavePeriod) && (lhs.minWaveAmplitude == rhs.minWaveAmplitude) &&
            (lhs.delayedHeaveCutoffFreq == rhs.delayedHeaveCutoffFreq) && (lhs.heaveCutoffFreq == rhs.heaveCutoffFreq) &&
            (lhs.heaveRateCutoffFreq == rhs.heaveRateCutoffFreq));
}

}  // namespace Heave

namespace IMU
{
bool Mag::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 3 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magX_t = StringUtils::fromString<float>(start, end);
    if (!magX_t.has_value()) return true;
    magX = magX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magY_t = StringUtils::fromString<float>(start, end);
    if (!magY_t.has_value()) return true;
    magY = magY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magZ_t = StringUtils::fromString<float>(start, end);
    if (!magZ_t.has_value()) return true;
    magZ = magZ_t.value();

    return false;  // no errors occurred
}

bool operator==(const Mag& lhs, const Mag& rhs) { return ((lhs.magX == rhs.magX) && (lhs.magY == rhs.magY) && (lhs.magZ == rhs.magZ)); }

bool Accel::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 3 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelX_t = StringUtils::fromString<float>(start, end);
    if (!accelX_t.has_value()) return true;
    accelX = accelX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelY_t = StringUtils::fromString<float>(start, end);
    if (!accelY_t.has_value()) return true;
    accelY = accelY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelZ_t = StringUtils::fromString<float>(start, end);
    if (!accelZ_t.has_value()) return true;
    accelZ = accelZ_t.value();

    return false;  // no errors occurred
}

bool operator==(const Accel& lhs, const Accel& rhs) { return ((lhs.accelX == rhs.accelX) && (lhs.accelY == rhs.accelY) && (lhs.accelZ == rhs.accelZ)); }

bool Gyro::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 3 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroX_t = StringUtils::fromString<float>(start, end);
    if (!gyroX_t.has_value()) return true;
    gyroX = gyroX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroY_t = StringUtils::fromString<float>(start, end);
    if (!gyroY_t.has_value()) return true;
    gyroY = gyroY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroZ_t = StringUtils::fromString<float>(start, end);
    if (!gyroZ_t.has_value()) return true;
    gyroZ = gyroZ_t.value();

    return false;  // no errors occurred
}

bool operator==(const Gyro& lhs, const Gyro& rhs) { return ((lhs.gyroX == rhs.gyroX) && (lhs.gyroY == rhs.gyroY) && (lhs.gyroZ == rhs.gyroZ)); }

bool MagAccelGyro::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 9 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magX_t = StringUtils::fromString<float>(start, end);
    if (!magX_t.has_value()) return true;
    magX = magX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magY_t = StringUtils::fromString<float>(start, end);
    if (!magY_t.has_value()) return true;
    magY = magY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magZ_t = StringUtils::fromString<float>(start, end);
    if (!magZ_t.has_value()) return true;
    magZ = magZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelX_t = StringUtils::fromString<float>(start, end);
    if (!accelX_t.has_value()) return true;
    accelX = accelX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelY_t = StringUtils::fromString<float>(start, end);
    if (!accelY_t.has_value()) return true;
    accelY = accelY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelZ_t = StringUtils::fromString<float>(start, end);
    if (!accelZ_t.has_value()) return true;
    accelZ = accelZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroX_t = StringUtils::fromString<float>(start, end);
    if (!gyroX_t.has_value()) return true;
    gyroX = gyroX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroY_t = StringUtils::fromString<float>(start, end);
    if (!gyroY_t.has_value()) return true;
    gyroY = gyroY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroZ_t = StringUtils::fromString<float>(start, end);
    if (!gyroZ_t.has_value()) return true;
    gyroZ = gyroZ_t.value();

    return false;  // no errors occurred
}

bool operator==(const MagAccelGyro& lhs, const MagAccelGyro& rhs)
{
    return ((lhs.magX == rhs.magX) && (lhs.magY == rhs.magY) && (lhs.magZ == rhs.magZ) && (lhs.accelX == rhs.accelX) && (lhs.accelY == rhs.accelY) &&
            (lhs.accelZ == rhs.accelZ) && (lhs.gyroX == rhs.gyroX) && (lhs.gyroY == rhs.gyroY) && (lhs.gyroZ == rhs.gyroZ));
}

bool MagCal::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 12 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magGain00_t = StringUtils::fromString<float>(start, end);
    if (!magGain00_t.has_value()) return true;
    magGain00 = magGain00_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magGain01_t = StringUtils::fromString<float>(start, end);
    if (!magGain01_t.has_value()) return true;
    magGain01 = magGain01_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magGain02_t = StringUtils::fromString<float>(start, end);
    if (!magGain02_t.has_value()) return true;
    magGain02 = magGain02_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magGain10_t = StringUtils::fromString<float>(start, end);
    if (!magGain10_t.has_value()) return true;
    magGain10 = magGain10_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magGain11_t = StringUtils::fromString<float>(start, end);
    if (!magGain11_t.has_value()) return true;
    magGain11 = magGain11_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magGain12_t = StringUtils::fromString<float>(start, end);
    if (!magGain12_t.has_value()) return true;
    magGain12 = magGain12_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magGain20_t = StringUtils::fromString<float>(start, end);
    if (!magGain20_t.has_value()) return true;
    magGain20 = magGain20_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magGain21_t = StringUtils::fromString<float>(start, end);
    if (!magGain21_t.has_value()) return true;
    magGain21 = magGain21_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magGain22_t = StringUtils::fromString<float>(start, end);
    if (!magGain22_t.has_value()) return true;
    magGain22 = magGain22_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magBiasX_t = StringUtils::fromString<float>(start, end);
    if (!magBiasX_t.has_value()) return true;
    magBiasX = magBiasX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magBiasY_t = StringUtils::fromString<float>(start, end);
    if (!magBiasY_t.has_value()) return true;
    magBiasY = magBiasY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magBiasZ_t = StringUtils::fromString<float>(start, end);
    if (!magBiasZ_t.has_value()) return true;
    magBiasZ = magBiasZ_t.value();

    return false;  // no errors occurred
}

AsciiMessage MagCal::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", magGain00, magGain01, magGain02, magGain10, magGain11, magGain12,
                  magGain20, magGain21, magGain22, magBiasX, magBiasY, magBiasZ);
    return result;
}

bool operator==(const MagCal& lhs, const MagCal& rhs)
{
    return ((lhs.magGain00 == rhs.magGain00) && (lhs.magGain01 == rhs.magGain01) && (lhs.magGain02 == rhs.magGain02) && (lhs.magGain10 == rhs.magGain10) &&
            (lhs.magGain11 == rhs.magGain11) && (lhs.magGain12 == rhs.magGain12) && (lhs.magGain20 == rhs.magGain20) && (lhs.magGain21 == rhs.magGain21) &&
            (lhs.magGain22 == rhs.magGain22) && (lhs.magBiasX == rhs.magBiasX) && (lhs.magBiasY == rhs.magBiasY) && (lhs.magBiasZ == rhs.magBiasZ));
}

bool AccelCal::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 12 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelGain00_t = StringUtils::fromString<float>(start, end);
    if (!accelGain00_t.has_value()) return true;
    accelGain00 = accelGain00_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelGain01_t = StringUtils::fromString<float>(start, end);
    if (!accelGain01_t.has_value()) return true;
    accelGain01 = accelGain01_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelGain02_t = StringUtils::fromString<float>(start, end);
    if (!accelGain02_t.has_value()) return true;
    accelGain02 = accelGain02_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelGain10_t = StringUtils::fromString<float>(start, end);
    if (!accelGain10_t.has_value()) return true;
    accelGain10 = accelGain10_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelGain11_t = StringUtils::fromString<float>(start, end);
    if (!accelGain11_t.has_value()) return true;
    accelGain11 = accelGain11_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelGain12_t = StringUtils::fromString<float>(start, end);
    if (!accelGain12_t.has_value()) return true;
    accelGain12 = accelGain12_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelGain20_t = StringUtils::fromString<float>(start, end);
    if (!accelGain20_t.has_value()) return true;
    accelGain20 = accelGain20_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelGain21_t = StringUtils::fromString<float>(start, end);
    if (!accelGain21_t.has_value()) return true;
    accelGain21 = accelGain21_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelGain22_t = StringUtils::fromString<float>(start, end);
    if (!accelGain22_t.has_value()) return true;
    accelGain22 = accelGain22_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelBiasX_t = StringUtils::fromString<float>(start, end);
    if (!accelBiasX_t.has_value()) return true;
    accelBiasX = accelBiasX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelBiasY_t = StringUtils::fromString<float>(start, end);
    if (!accelBiasY_t.has_value()) return true;
    accelBiasY = accelBiasY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelBiasZ_t = StringUtils::fromString<float>(start, end);
    if (!accelBiasZ_t.has_value()) return true;
    accelBiasZ = accelBiasZ_t.value();

    return false;  // no errors occurred
}

AsciiMessage AccelCal::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", accelGain00, accelGain01, accelGain02, accelGain10, accelGain11,
                  accelGain12, accelGain20, accelGain21, accelGain22, accelBiasX, accelBiasY, accelBiasZ);
    return result;
}

bool operator==(const AccelCal& lhs, const AccelCal& rhs)
{
    return ((lhs.accelGain00 == rhs.accelGain00) && (lhs.accelGain01 == rhs.accelGain01) && (lhs.accelGain02 == rhs.accelGain02) &&
            (lhs.accelGain10 == rhs.accelGain10) && (lhs.accelGain11 == rhs.accelGain11) && (lhs.accelGain12 == rhs.accelGain12) &&
            (lhs.accelGain20 == rhs.accelGain20) && (lhs.accelGain21 == rhs.accelGain21) && (lhs.accelGain22 == rhs.accelGain22) &&
            (lhs.accelBiasX == rhs.accelBiasX) && (lhs.accelBiasY == rhs.accelBiasY) && (lhs.accelBiasZ == rhs.accelBiasZ));
}

bool RefFrameRot::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 9 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto rFR00_t = StringUtils::fromString<float>(start, end);
    if (!rFR00_t.has_value()) return true;
    rFR00 = rFR00_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto rFR01_t = StringUtils::fromString<float>(start, end);
    if (!rFR01_t.has_value()) return true;
    rFR01 = rFR01_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto rFR02_t = StringUtils::fromString<float>(start, end);
    if (!rFR02_t.has_value()) return true;
    rFR02 = rFR02_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto rFR10_t = StringUtils::fromString<float>(start, end);
    if (!rFR10_t.has_value()) return true;
    rFR10 = rFR10_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto rFR11_t = StringUtils::fromString<float>(start, end);
    if (!rFR11_t.has_value()) return true;
    rFR11 = rFR11_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto rFR12_t = StringUtils::fromString<float>(start, end);
    if (!rFR12_t.has_value()) return true;
    rFR12 = rFR12_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto rFR20_t = StringUtils::fromString<float>(start, end);
    if (!rFR20_t.has_value()) return true;
    rFR20 = rFR20_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto rFR21_t = StringUtils::fromString<float>(start, end);
    if (!rFR21_t.has_value()) return true;
    rFR21 = rFR21_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto rFR22_t = StringUtils::fromString<float>(start, end);
    if (!rFR22_t.has_value()) return true;
    rFR22 = rFR22_t.value();

    return false;  // no errors occurred
}

AsciiMessage RefFrameRot::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%f,%f,%f,%f,%f,%f,%f,%f,%f", rFR00, rFR01, rFR02, rFR10, rFR11, rFR12, rFR20, rFR21, rFR22);
    return result;
}

bool operator==(const RefFrameRot& lhs, const RefFrameRot& rhs)
{
    return ((lhs.rFR00 == rhs.rFR00) && (lhs.rFR01 == rhs.rFR01) && (lhs.rFR02 == rhs.rFR02) && (lhs.rFR10 == rhs.rFR10) && (lhs.rFR11 == rhs.rFR11) &&
            (lhs.rFR12 == rhs.rFR12) && (lhs.rFR20 == rhs.rFR20) && (lhs.rFR21 == rhs.rFR21) && (lhs.rFR22 == rhs.rFR22));
}

bool ImuMeas::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 11 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto uncompMagX_t = StringUtils::fromString<float>(start, end);
    if (!uncompMagX_t.has_value()) return true;
    uncompMagX = uncompMagX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto uncompMagY_t = StringUtils::fromString<float>(start, end);
    if (!uncompMagY_t.has_value()) return true;
    uncompMagY = uncompMagY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto uncompMagZ_t = StringUtils::fromString<float>(start, end);
    if (!uncompMagZ_t.has_value()) return true;
    uncompMagZ = uncompMagZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto uncompAccX_t = StringUtils::fromString<float>(start, end);
    if (!uncompAccX_t.has_value()) return true;
    uncompAccX = uncompAccX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto uncompAccY_t = StringUtils::fromString<float>(start, end);
    if (!uncompAccY_t.has_value()) return true;
    uncompAccY = uncompAccY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto uncompAccZ_t = StringUtils::fromString<float>(start, end);
    if (!uncompAccZ_t.has_value()) return true;
    uncompAccZ = uncompAccZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto uncompGyroX_t = StringUtils::fromString<float>(start, end);
    if (!uncompGyroX_t.has_value()) return true;
    uncompGyroX = uncompGyroX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto uncompGyroY_t = StringUtils::fromString<float>(start, end);
    if (!uncompGyroY_t.has_value()) return true;
    uncompGyroY = uncompGyroY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto uncompGyroZ_t = StringUtils::fromString<float>(start, end);
    if (!uncompGyroZ_t.has_value()) return true;
    uncompGyroZ = uncompGyroZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto temperature_t = StringUtils::fromString<float>(start, end);
    if (!temperature_t.has_value()) return true;
    temperature = temperature_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto pressure_t = StringUtils::fromString<float>(start, end);
    if (!pressure_t.has_value()) return true;
    pressure = pressure_t.value();

    return false;  // no errors occurred
}

bool operator==(const ImuMeas& lhs, const ImuMeas& rhs)
{
    return ((lhs.uncompMagX == rhs.uncompMagX) && (lhs.uncompMagY == rhs.uncompMagY) && (lhs.uncompMagZ == rhs.uncompMagZ) &&
            (lhs.uncompAccX == rhs.uncompAccX) && (lhs.uncompAccY == rhs.uncompAccY) && (lhs.uncompAccZ == rhs.uncompAccZ) &&
            (lhs.uncompGyroX == rhs.uncompGyroX) && (lhs.uncompGyroY == rhs.uncompGyroY) && (lhs.uncompGyroZ == rhs.uncompGyroZ) &&
            (lhs.temperature == rhs.temperature) && (lhs.pressure == rhs.pressure));
}

bool DeltaThetaVelocity::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 7 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto deltaTime_t = StringUtils::fromString<float>(start, end);
    if (!deltaTime_t.has_value()) return true;
    deltaTime = deltaTime_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto deltaThetaX_t = StringUtils::fromString<float>(start, end);
    if (!deltaThetaX_t.has_value()) return true;
    deltaThetaX = deltaThetaX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto deltaThetaY_t = StringUtils::fromString<float>(start, end);
    if (!deltaThetaY_t.has_value()) return true;
    deltaThetaY = deltaThetaY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto deltaThetaZ_t = StringUtils::fromString<float>(start, end);
    if (!deltaThetaZ_t.has_value()) return true;
    deltaThetaZ = deltaThetaZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto deltaVelX_t = StringUtils::fromString<float>(start, end);
    if (!deltaVelX_t.has_value()) return true;
    deltaVelX = deltaVelX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto deltaVelY_t = StringUtils::fromString<float>(start, end);
    if (!deltaVelY_t.has_value()) return true;
    deltaVelY = deltaVelY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto deltaVelZ_t = StringUtils::fromString<float>(start, end);
    if (!deltaVelZ_t.has_value()) return true;
    deltaVelZ = deltaVelZ_t.value();

    return false;  // no errors occurred
}

bool operator==(const DeltaThetaVelocity& lhs, const DeltaThetaVelocity& rhs)
{
    return ((lhs.deltaTime == rhs.deltaTime) && (lhs.deltaThetaX == rhs.deltaThetaX) && (lhs.deltaThetaY == rhs.deltaThetaY) &&
            (lhs.deltaThetaZ == rhs.deltaThetaZ) && (lhs.deltaVelX == rhs.deltaVelX) && (lhs.deltaVelY == rhs.deltaVelY) && (lhs.deltaVelZ == rhs.deltaVelZ));
}

bool DeltaThetaVelConfig::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 5 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto integrationFrame_t = StringUtils::fromString<uint8_t>(start, end);
    if (!integrationFrame_t.has_value()) return true;
    integrationFrame = static_cast<DeltaThetaVelConfig::IntegrationFrame>(integrationFrame_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroCompensation_t = StringUtils::fromString<uint8_t>(start, end);
    if (!gyroCompensation_t.has_value()) return true;
    gyroCompensation = static_cast<DeltaThetaVelConfig::GyroCompensation>(gyroCompensation_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelCompensation_t = StringUtils::fromString<uint8_t>(start, end);
    if (!accelCompensation_t.has_value()) return true;
    accelCompensation = static_cast<DeltaThetaVelConfig::AccelCompensation>(accelCompensation_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto earthRateCompensation_t = StringUtils::fromString<uint8_t>(start, end);
    if (!earthRateCompensation_t.has_value()) return true;
    earthRateCompensation = static_cast<DeltaThetaVelConfig::EarthRateCompensation>(earthRateCompensation_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto resv_t = StringUtils::fromString<uint16_t>(start, end);
    if (!resv_t.has_value()) return true;
    resv = resv_t.value();

    return false;  // no errors occurred
}

AsciiMessage DeltaThetaVelConfig::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%u,%u,%u,%u,%u", static_cast<uint8_t>(integrationFrame), static_cast<uint8_t>(gyroCompensation),
                  static_cast<uint8_t>(accelCompensation), static_cast<uint8_t>(earthRateCompensation), resv);
    return result;
}

bool operator==(const DeltaThetaVelConfig& lhs, const DeltaThetaVelConfig& rhs)
{
    return ((lhs.integrationFrame == rhs.integrationFrame) && (lhs.gyroCompensation == rhs.gyroCompensation) &&
            (lhs.accelCompensation == rhs.accelCompensation) && (lhs.earthRateCompensation == rhs.earthRateCompensation) && (lhs.resv == rhs.resv));
}

bool GyroCal::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 12 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroGain00_t = StringUtils::fromString<float>(start, end);
    if (!gyroGain00_t.has_value()) return true;
    gyroGain00 = gyroGain00_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroGain01_t = StringUtils::fromString<float>(start, end);
    if (!gyroGain01_t.has_value()) return true;
    gyroGain01 = gyroGain01_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroGain02_t = StringUtils::fromString<float>(start, end);
    if (!gyroGain02_t.has_value()) return true;
    gyroGain02 = gyroGain02_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroGain10_t = StringUtils::fromString<float>(start, end);
    if (!gyroGain10_t.has_value()) return true;
    gyroGain10 = gyroGain10_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroGain11_t = StringUtils::fromString<float>(start, end);
    if (!gyroGain11_t.has_value()) return true;
    gyroGain11 = gyroGain11_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroGain12_t = StringUtils::fromString<float>(start, end);
    if (!gyroGain12_t.has_value()) return true;
    gyroGain12 = gyroGain12_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroGain20_t = StringUtils::fromString<float>(start, end);
    if (!gyroGain20_t.has_value()) return true;
    gyroGain20 = gyroGain20_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroGain21_t = StringUtils::fromString<float>(start, end);
    if (!gyroGain21_t.has_value()) return true;
    gyroGain21 = gyroGain21_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroGain22_t = StringUtils::fromString<float>(start, end);
    if (!gyroGain22_t.has_value()) return true;
    gyroGain22 = gyroGain22_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroBiasX_t = StringUtils::fromString<float>(start, end);
    if (!gyroBiasX_t.has_value()) return true;
    gyroBiasX = gyroBiasX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroBiasY_t = StringUtils::fromString<float>(start, end);
    if (!gyroBiasY_t.has_value()) return true;
    gyroBiasY = gyroBiasY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroBiasZ_t = StringUtils::fromString<float>(start, end);
    if (!gyroBiasZ_t.has_value()) return true;
    gyroBiasZ = gyroBiasZ_t.value();

    return false;  // no errors occurred
}

AsciiMessage GyroCal::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f", gyroGain00, gyroGain01, gyroGain02, gyroGain10, gyroGain11,
                  gyroGain12, gyroGain20, gyroGain21, gyroGain22, gyroBiasX, gyroBiasY, gyroBiasZ);
    return result;
}

bool operator==(const GyroCal& lhs, const GyroCal& rhs)
{
    return ((lhs.gyroGain00 == rhs.gyroGain00) && (lhs.gyroGain01 == rhs.gyroGain01) && (lhs.gyroGain02 == rhs.gyroGain02) &&
            (lhs.gyroGain10 == rhs.gyroGain10) && (lhs.gyroGain11 == rhs.gyroGain11) && (lhs.gyroGain12 == rhs.gyroGain12) &&
            (lhs.gyroGain20 == rhs.gyroGain20) && (lhs.gyroGain21 == rhs.gyroGain21) && (lhs.gyroGain22 == rhs.gyroGain22) &&
            (lhs.gyroBiasX == rhs.gyroBiasX) && (lhs.gyroBiasY == rhs.gyroBiasY) && (lhs.gyroBiasZ == rhs.gyroBiasZ));
}

bool ImuFilterControl::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 10 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magWindowSize_t = StringUtils::fromString<uint16_t>(start, end);
    if (!magWindowSize_t.has_value()) return true;
    magWindowSize = magWindowSize_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelWindowSize_t = StringUtils::fromString<uint16_t>(start, end);
    if (!accelWindowSize_t.has_value()) return true;
    accelWindowSize = accelWindowSize_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroWindowSize_t = StringUtils::fromString<uint16_t>(start, end);
    if (!gyroWindowSize_t.has_value()) return true;
    gyroWindowSize = gyroWindowSize_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto tempWindowSize_t = StringUtils::fromString<uint16_t>(start, end);
    if (!tempWindowSize_t.has_value()) return true;
    tempWindowSize = tempWindowSize_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto presWindowSize_t = StringUtils::fromString<uint16_t>(start, end);
    if (!presWindowSize_t.has_value()) return true;
    presWindowSize = presWindowSize_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto magFilterMode_t = StringUtils::fromString<uint8_t>(start, end);
    if (!magFilterMode_t.has_value()) return true;
    magFilterMode._value = magFilterMode_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelFilterMode_t = StringUtils::fromString<uint8_t>(start, end);
    if (!accelFilterMode_t.has_value()) return true;
    accelFilterMode._value = accelFilterMode_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroFilterMode_t = StringUtils::fromString<uint8_t>(start, end);
    if (!gyroFilterMode_t.has_value()) return true;
    gyroFilterMode._value = gyroFilterMode_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto tempFilterMode_t = StringUtils::fromString<uint8_t>(start, end);
    if (!tempFilterMode_t.has_value()) return true;
    tempFilterMode._value = tempFilterMode_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto presFilterMode_t = StringUtils::fromString<uint8_t>(start, end);
    if (!presFilterMode_t.has_value()) return true;
    presFilterMode._value = presFilterMode_t.value();

    return false;  // no errors occurred
}

AsciiMessage ImuFilterControl::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%u,%u,%u,%u,%u,%u,%u,%u,%u,%u", magWindowSize, accelWindowSize, gyroWindowSize, tempWindowSize,
                  presWindowSize, magFilterMode._value, accelFilterMode._value, gyroFilterMode._value, tempFilterMode._value, presFilterMode._value);
    return result;
}

bool operator==(const ImuFilterControl& lhs, const ImuFilterControl& rhs)
{
    return ((lhs.magWindowSize == rhs.magWindowSize) && (lhs.accelWindowSize == rhs.accelWindowSize) && (lhs.gyroWindowSize == rhs.gyroWindowSize) &&
            (lhs.tempWindowSize == rhs.tempWindowSize) && (lhs.presWindowSize == rhs.presWindowSize) &&
            (lhs.magFilterMode._value == rhs.magFilterMode._value) && (lhs.accelFilterMode._value == rhs.accelFilterMode._value) &&
            (lhs.gyroFilterMode._value == rhs.gyroFilterMode._value) && (lhs.tempFilterMode._value == rhs.tempFilterMode._value) &&
            (lhs.presFilterMode._value == rhs.presFilterMode._value));
}

}  // namespace IMU

namespace INS
{
bool InsSolLla::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 15 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gpsTow_t = StringUtils::fromString<double>(start, end);
    if (!gpsTow_t.has_value()) return true;
    gpsTow = gpsTow_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto timeGpsWeek_t = StringUtils::fromString<uint16_t>(start, end);
    if (!timeGpsWeek_t.has_value()) return true;
    timeGpsWeek = timeGpsWeek_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto insStatus_t = StringUtils::fromStringHex<uint16_t>(start, end);
    if (!insStatus_t.has_value()) return true;
    insStatus._value = insStatus_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto yaw_t = StringUtils::fromString<float>(start, end);
    if (!yaw_t.has_value()) return true;
    yaw = yaw_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto pitch_t = StringUtils::fromString<float>(start, end);
    if (!pitch_t.has_value()) return true;
    pitch = pitch_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto roll_t = StringUtils::fromString<float>(start, end);
    if (!roll_t.has_value()) return true;
    roll = roll_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posLat_t = StringUtils::fromString<double>(start, end);
    if (!posLat_t.has_value()) return true;
    posLat = posLat_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posLon_t = StringUtils::fromString<double>(start, end);
    if (!posLon_t.has_value()) return true;
    posLon = posLon_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posAlt_t = StringUtils::fromString<double>(start, end);
    if (!posAlt_t.has_value()) return true;
    posAlt = posAlt_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velN_t = StringUtils::fromString<float>(start, end);
    if (!velN_t.has_value()) return true;
    velN = velN_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velE_t = StringUtils::fromString<float>(start, end);
    if (!velE_t.has_value()) return true;
    velE = velE_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velD_t = StringUtils::fromString<float>(start, end);
    if (!velD_t.has_value()) return true;
    velD = velD_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto attUncertainty_t = StringUtils::fromString<float>(start, end);
    if (!attUncertainty_t.has_value()) return true;
    attUncertainty = attUncertainty_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posUncertainty_t = StringUtils::fromString<float>(start, end);
    if (!posUncertainty_t.has_value()) return true;
    posUncertainty = posUncertainty_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velUncertainty_t = StringUtils::fromString<float>(start, end);
    if (!velUncertainty_t.has_value()) return true;
    velUncertainty = velUncertainty_t.value();

    return false;  // no errors occurred
}

bool operator==(const InsSolLla& lhs, const InsSolLla& rhs)
{
    return ((lhs.gpsTow == rhs.gpsTow) && (lhs.timeGpsWeek == rhs.timeGpsWeek) && (lhs.insStatus._value == rhs.insStatus._value) && (lhs.yaw == rhs.yaw) &&
            (lhs.pitch == rhs.pitch) && (lhs.roll == rhs.roll) && (lhs.posLat == rhs.posLat) && (lhs.posLon == rhs.posLon) && (lhs.posAlt == rhs.posAlt) &&
            (lhs.velN == rhs.velN) && (lhs.velE == rhs.velE) && (lhs.velD == rhs.velD) && (lhs.attUncertainty == rhs.attUncertainty) &&
            (lhs.posUncertainty == rhs.posUncertainty) && (lhs.velUncertainty == rhs.velUncertainty));
}

bool InsSolEcef::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 15 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gpsTow_t = StringUtils::fromString<double>(start, end);
    if (!gpsTow_t.has_value()) return true;
    gpsTow = gpsTow_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto timeGpsWeek_t = StringUtils::fromString<uint16_t>(start, end);
    if (!timeGpsWeek_t.has_value()) return true;
    timeGpsWeek = timeGpsWeek_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto insStatus_t = StringUtils::fromStringHex<uint16_t>(start, end);
    if (!insStatus_t.has_value()) return true;
    insStatus._value = insStatus_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto yaw_t = StringUtils::fromString<float>(start, end);
    if (!yaw_t.has_value()) return true;
    yaw = yaw_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto pitch_t = StringUtils::fromString<float>(start, end);
    if (!pitch_t.has_value()) return true;
    pitch = pitch_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto roll_t = StringUtils::fromString<float>(start, end);
    if (!roll_t.has_value()) return true;
    roll = roll_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posEX_t = StringUtils::fromString<double>(start, end);
    if (!posEX_t.has_value()) return true;
    posEX = posEX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posEY_t = StringUtils::fromString<double>(start, end);
    if (!posEY_t.has_value()) return true;
    posEY = posEY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posEZ_t = StringUtils::fromString<double>(start, end);
    if (!posEZ_t.has_value()) return true;
    posEZ = posEZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velEX_t = StringUtils::fromString<float>(start, end);
    if (!velEX_t.has_value()) return true;
    velEX = velEX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velEY_t = StringUtils::fromString<float>(start, end);
    if (!velEY_t.has_value()) return true;
    velEY = velEY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velEZ_t = StringUtils::fromString<float>(start, end);
    if (!velEZ_t.has_value()) return true;
    velEZ = velEZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto attUncertainty_t = StringUtils::fromString<float>(start, end);
    if (!attUncertainty_t.has_value()) return true;
    attUncertainty = attUncertainty_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posUncertainty_t = StringUtils::fromString<float>(start, end);
    if (!posUncertainty_t.has_value()) return true;
    posUncertainty = posUncertainty_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velUncertainty_t = StringUtils::fromString<float>(start, end);
    if (!velUncertainty_t.has_value()) return true;
    velUncertainty = velUncertainty_t.value();

    return false;  // no errors occurred
}

bool operator==(const InsSolEcef& lhs, const InsSolEcef& rhs)
{
    return ((lhs.gpsTow == rhs.gpsTow) && (lhs.timeGpsWeek == rhs.timeGpsWeek) && (lhs.insStatus._value == rhs.insStatus._value) && (lhs.yaw == rhs.yaw) &&
            (lhs.pitch == rhs.pitch) && (lhs.roll == rhs.roll) && (lhs.posEX == rhs.posEX) && (lhs.posEY == rhs.posEY) && (lhs.posEZ == rhs.posEZ) &&
            (lhs.velEX == rhs.velEX) && (lhs.velEY == rhs.velEY) && (lhs.velEZ == rhs.velEZ) && (lhs.attUncertainty == rhs.attUncertainty) &&
            (lhs.posUncertainty == rhs.posUncertainty) && (lhs.velUncertainty == rhs.velUncertainty));
}

bool InsBasicConfig::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 4 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto scenario_t = StringUtils::fromString<uint8_t>(start, end);
    if (!scenario_t.has_value()) return true;
    scenario = static_cast<InsBasicConfig::Scenario>(scenario_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto ahrsAiding_t = StringUtils::fromString<uint8_t>(start, end);
    if (!ahrsAiding_t.has_value()) return true;
    ahrsAiding = static_cast<InsBasicConfig::AhrsAiding>(ahrsAiding_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto estBaseline_t = StringUtils::fromString<uint8_t>(start, end);
    if (!estBaseline_t.has_value()) return true;
    estBaseline = static_cast<InsBasicConfig::EstBaseline>(estBaseline_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto resv_t = StringUtils::fromString<uint8_t>(start, end);
    if (!resv_t.has_value()) return true;
    resv = resv_t.value();

    return false;  // no errors occurred
}

AsciiMessage InsBasicConfig::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%u,%u,%u,%u", static_cast<uint8_t>(scenario), static_cast<uint8_t>(ahrsAiding),
                  static_cast<uint8_t>(estBaseline), resv);
    return result;
}

bool operator==(const InsBasicConfig& lhs, const InsBasicConfig& rhs)
{
    return ((lhs.scenario == rhs.scenario) && (lhs.ahrsAiding == rhs.ahrsAiding) && (lhs.estBaseline == rhs.estBaseline) && (lhs.resv == rhs.resv));
}

bool InsStateLla::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 15 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto yaw_t = StringUtils::fromString<float>(start, end);
    if (!yaw_t.has_value()) return true;
    yaw = yaw_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto pitch_t = StringUtils::fromString<float>(start, end);
    if (!pitch_t.has_value()) return true;
    pitch = pitch_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto roll_t = StringUtils::fromString<float>(start, end);
    if (!roll_t.has_value()) return true;
    roll = roll_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posLat_t = StringUtils::fromString<double>(start, end);
    if (!posLat_t.has_value()) return true;
    posLat = posLat_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posLon_t = StringUtils::fromString<double>(start, end);
    if (!posLon_t.has_value()) return true;
    posLon = posLon_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posAlt_t = StringUtils::fromString<double>(start, end);
    if (!posAlt_t.has_value()) return true;
    posAlt = posAlt_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velN_t = StringUtils::fromString<float>(start, end);
    if (!velN_t.has_value()) return true;
    velN = velN_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velE_t = StringUtils::fromString<float>(start, end);
    if (!velE_t.has_value()) return true;
    velE = velE_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velD_t = StringUtils::fromString<float>(start, end);
    if (!velD_t.has_value()) return true;
    velD = velD_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelX_t = StringUtils::fromString<float>(start, end);
    if (!accelX_t.has_value()) return true;
    accelX = accelX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelY_t = StringUtils::fromString<float>(start, end);
    if (!accelY_t.has_value()) return true;
    accelY = accelY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelZ_t = StringUtils::fromString<float>(start, end);
    if (!accelZ_t.has_value()) return true;
    accelZ = accelZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroX_t = StringUtils::fromString<float>(start, end);
    if (!gyroX_t.has_value()) return true;
    gyroX = gyroX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroY_t = StringUtils::fromString<float>(start, end);
    if (!gyroY_t.has_value()) return true;
    gyroY = gyroY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroZ_t = StringUtils::fromString<float>(start, end);
    if (!gyroZ_t.has_value()) return true;
    gyroZ = gyroZ_t.value();

    return false;  // no errors occurred
}

bool operator==(const InsStateLla& lhs, const InsStateLla& rhs)
{
    return ((lhs.yaw == rhs.yaw) && (lhs.pitch == rhs.pitch) && (lhs.roll == rhs.roll) && (lhs.posLat == rhs.posLat) && (lhs.posLon == rhs.posLon) &&
            (lhs.posAlt == rhs.posAlt) && (lhs.velN == rhs.velN) && (lhs.velE == rhs.velE) && (lhs.velD == rhs.velD) && (lhs.accelX == rhs.accelX) &&
            (lhs.accelY == rhs.accelY) && (lhs.accelZ == rhs.accelZ) && (lhs.gyroX == rhs.gyroX) && (lhs.gyroY == rhs.gyroY) && (lhs.gyroZ == rhs.gyroZ));
}

bool InsStateEcef::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 15 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto yaw_t = StringUtils::fromString<float>(start, end);
    if (!yaw_t.has_value()) return true;
    yaw = yaw_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto pitch_t = StringUtils::fromString<float>(start, end);
    if (!pitch_t.has_value()) return true;
    pitch = pitch_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto roll_t = StringUtils::fromString<float>(start, end);
    if (!roll_t.has_value()) return true;
    roll = roll_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posEX_t = StringUtils::fromString<double>(start, end);
    if (!posEX_t.has_value()) return true;
    posEX = posEX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posEY_t = StringUtils::fromString<double>(start, end);
    if (!posEY_t.has_value()) return true;
    posEY = posEY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto posEZ_t = StringUtils::fromString<double>(start, end);
    if (!posEZ_t.has_value()) return true;
    posEZ = posEZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velEX_t = StringUtils::fromString<float>(start, end);
    if (!velEX_t.has_value()) return true;
    velEX = velEX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velEY_t = StringUtils::fromString<float>(start, end);
    if (!velEY_t.has_value()) return true;
    velEY = velEY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velEZ_t = StringUtils::fromString<float>(start, end);
    if (!velEZ_t.has_value()) return true;
    velEZ = velEZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelX_t = StringUtils::fromString<float>(start, end);
    if (!accelX_t.has_value()) return true;
    accelX = accelX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelY_t = StringUtils::fromString<float>(start, end);
    if (!accelY_t.has_value()) return true;
    accelY = accelY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelZ_t = StringUtils::fromString<float>(start, end);
    if (!accelZ_t.has_value()) return true;
    accelZ = accelZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroX_t = StringUtils::fromString<float>(start, end);
    if (!gyroX_t.has_value()) return true;
    gyroX = gyroX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroY_t = StringUtils::fromString<float>(start, end);
    if (!gyroY_t.has_value()) return true;
    gyroY = gyroY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroZ_t = StringUtils::fromString<float>(start, end);
    if (!gyroZ_t.has_value()) return true;
    gyroZ = gyroZ_t.value();

    return false;  // no errors occurred
}

bool operator==(const InsStateEcef& lhs, const InsStateEcef& rhs)
{
    return ((lhs.yaw == rhs.yaw) && (lhs.pitch == rhs.pitch) && (lhs.roll == rhs.roll) && (lhs.posEX == rhs.posEX) && (lhs.posEY == rhs.posEY) &&
            (lhs.posEZ == rhs.posEZ) && (lhs.velEX == rhs.velEX) && (lhs.velEY == rhs.velEY) && (lhs.velEZ == rhs.velEZ) && (lhs.accelX == rhs.accelX) &&
            (lhs.accelY == rhs.accelY) && (lhs.accelZ == rhs.accelZ) && (lhs.gyroX == rhs.gyroX) && (lhs.gyroY == rhs.gyroY) && (lhs.gyroZ == rhs.gyroZ));
}

bool FilterStartupBias::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 7 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroBiasX_t = StringUtils::fromString<float>(start, end);
    if (!gyroBiasX_t.has_value()) return true;
    gyroBiasX = gyroBiasX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroBiasY_t = StringUtils::fromString<float>(start, end);
    if (!gyroBiasY_t.has_value()) return true;
    gyroBiasY = gyroBiasY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gyroBiasZ_t = StringUtils::fromString<float>(start, end);
    if (!gyroBiasZ_t.has_value()) return true;
    gyroBiasZ = gyroBiasZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelBiasX_t = StringUtils::fromString<float>(start, end);
    if (!accelBiasX_t.has_value()) return true;
    accelBiasX = accelBiasX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelBiasY_t = StringUtils::fromString<float>(start, end);
    if (!accelBiasY_t.has_value()) return true;
    accelBiasY = accelBiasY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto accelBiasZ_t = StringUtils::fromString<float>(start, end);
    if (!accelBiasZ_t.has_value()) return true;
    accelBiasZ = accelBiasZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto presBias_t = StringUtils::fromString<float>(start, end);
    if (!presBias_t.has_value()) return true;
    presBias = presBias_t.value();

    return false;  // no errors occurred
}

AsciiMessage FilterStartupBias::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%f,%f,%f,%f,%f,%f,%f", gyroBiasX, gyroBiasY, gyroBiasZ, accelBiasX, accelBiasY, accelBiasZ, presBias);
    return result;
}

bool operator==(const FilterStartupBias& lhs, const FilterStartupBias& rhs)
{
    return ((lhs.gyroBiasX == rhs.gyroBiasX) && (lhs.gyroBiasY == rhs.gyroBiasY) && (lhs.gyroBiasZ == rhs.gyroBiasZ) && (lhs.accelBiasX == rhs.accelBiasX) &&
            (lhs.accelBiasY == rhs.accelBiasY) && (lhs.accelBiasZ == rhs.accelBiasZ) && (lhs.presBias == rhs.presBias));
}

bool InsRefOffset::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 6 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto refOffsetX_t = StringUtils::fromString<float>(start, end);
    if (!refOffsetX_t.has_value()) return true;
    refOffsetX = refOffsetX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto refOffsetY_t = StringUtils::fromString<float>(start, end);
    if (!refOffsetY_t.has_value()) return true;
    refOffsetY = refOffsetY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto refOffsetZ_t = StringUtils::fromString<float>(start, end);
    if (!refOffsetZ_t.has_value()) return true;
    refOffsetZ = refOffsetZ_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto refUncertX_t = StringUtils::fromString<float>(start, end);
    if (!refUncertX_t.has_value()) return true;
    refUncertX = refUncertX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto refUncertY_t = StringUtils::fromString<float>(start, end);
    if (!refUncertY_t.has_value()) return true;
    refUncertY = refUncertY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto refUncertZ_t = StringUtils::fromString<float>(start, end);
    if (!refUncertZ_t.has_value()) return true;
    refUncertZ = refUncertZ_t.value();

    return false;  // no errors occurred
}

AsciiMessage InsRefOffset::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%f,%f,%f,%f,%f,%f", refOffsetX, refOffsetY, refOffsetZ, refUncertX, refUncertY, refUncertZ);
    return result;
}

bool operator==(const InsRefOffset& lhs, const InsRefOffset& rhs)
{
    return ((lhs.refOffsetX == rhs.refOffsetX) && (lhs.refOffsetY == rhs.refOffsetY) && (lhs.refOffsetZ == rhs.refOffsetZ) &&
            (lhs.refUncertX == rhs.refUncertX) && (lhs.refUncertY == rhs.refUncertY) && (lhs.refUncertZ == rhs.refUncertZ));
}

bool InsGnssSelect::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 6 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto activeReceiverSelect_t = StringUtils::fromString<uint8_t>(start, end);
    if (!activeReceiverSelect_t.has_value()) return true;
    activeReceiverSelect = static_cast<InsGnssSelect::ActiveReceiverSelect>(activeReceiverSelect_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto usedForNavTime_t = StringUtils::fromString<uint8_t>(start, end);
    if (!usedForNavTime_t.has_value()) return true;
    usedForNavTime = usedForNavTime_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto hysteresisTime_t = StringUtils::fromString<uint8_t>(start, end);
    if (!hysteresisTime_t.has_value()) return true;
    hysteresisTime = hysteresisTime_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto useGnssCompass_t = StringUtils::fromString<uint8_t>(start, end);
    if (!useGnssCompass_t.has_value()) return true;
    useGnssCompass = static_cast<InsGnssSelect::UseGnssCompass>(useGnssCompass_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto resv1_t = StringUtils::fromString<uint8_t>(start, end);
    if (!resv1_t.has_value()) return true;
    resv1 = resv1_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto resv2_t = StringUtils::fromString<uint8_t>(start, end);
    if (!resv2_t.has_value()) return true;
    resv2 = resv2_t.value();

    return false;  // no errors occurred
}

AsciiMessage InsGnssSelect::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%u,%u,%u,%u,%u,%u", static_cast<uint8_t>(activeReceiverSelect), usedForNavTime, hysteresisTime,
                  static_cast<uint8_t>(useGnssCompass), resv1, resv2);
    return result;
}

bool operator==(const InsGnssSelect& lhs, const InsGnssSelect& rhs)
{
    return ((lhs.activeReceiverSelect == rhs.activeReceiverSelect) && (lhs.usedForNavTime == rhs.usedForNavTime) &&
            (lhs.hysteresisTime == rhs.hysteresisTime) && (lhs.useGnssCompass == rhs.useGnssCompass) && (lhs.resv1 == rhs.resv1) && (lhs.resv2 == rhs.resv2));
}

}  // namespace INS

namespace System
{
bool UserTag::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 1 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    tag = AsciiMessage(start, end);

    return false;  // no errors occurred
}

AsciiMessage UserTag::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%s", tag.c_str());
    return result;
}

bool operator==(const UserTag& lhs, const UserTag& rhs) { return ((lhs.tag == rhs.tag)); }

bool Model::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 1 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    model = AsciiMessage(start, end);

    return false;  // no errors occurred
}

bool operator==(const Model& lhs, const Model& rhs) { return ((lhs.model == rhs.model)); }

bool HwVer::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 2 + 1)
    {
        // This register may not have an optional field
        if (tokens.size() != 1 + 1) { return true; }
    }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto hwVer_t = StringUtils::fromString<uint32_t>(start, end);
    if (!hwVer_t.has_value()) return true;
    hwVer = hwVer_t.value();

    if (tokens.size() == (2 + 1))
    {
        // This is an optional parameter.
        start = tokens[index++] + 1;
        end = tokens[index];
        const auto hwMinVer_t = StringUtils::fromString<uint32_t>(start, end);
        if (!hwMinVer_t.has_value()) return true;
        hwMinVer = hwMinVer_t.value();
    }

    return false;  // no errors occurred
}

bool operator==(const HwVer& lhs, const HwVer& rhs) { return ((lhs.hwVer == rhs.hwVer) && (lhs.hwMinVer == rhs.hwMinVer)); }

bool Serial::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 1 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto serialNum_t = StringUtils::fromString<uint32_t>(start, end);
    if (!serialNum_t.has_value()) return true;
    serialNum = serialNum_t.value();

    return false;  // no errors occurred
}

bool operator==(const Serial& lhs, const Serial& rhs) { return ((lhs.serialNum == rhs.serialNum)); }

bool FwVer::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 1 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    fwVer = AsciiMessage(start, end);

    return false;  // no errors occurred
}

bool operator==(const FwVer& lhs, const FwVer& rhs) { return ((lhs.fwVer == rhs.fwVer)); }

bool BaudRate::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 2 + 1)
    {
        // This register may not have an optional field
        if (tokens.size() != 1 + 1) { return true; }
    }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto baudRate_t = StringUtils::fromString<uint32_t>(start, end);
    if (!baudRate_t.has_value()) return true;
    baudRate = static_cast<BaudRate::BaudRates>(baudRate_t.value());

    if (tokens.size() == (2 + 1))
    {
        // This is an optional parameter.
        start = tokens[index++] + 1;
        end = tokens[index];
        const auto serialPort_t = StringUtils::fromString<uint8_t>(start, end);
        if (!serialPort_t.has_value()) return true;
        serialPort = static_cast<BaudRate::SerialPort>(serialPort_t.value());
    }

    return false;  // no errors occurred
}

AsciiMessage BaudRate::toString() const
{
    AsciiMessage result = "";
    if (static_cast<uint8_t>(serialPort) == 0) { std::snprintf(result.begin(), result.capacity(), "%d", static_cast<uint32_t>(baudRate)); }
    else
    {
        if (serialPort == SerialPort::Poll) { std::snprintf(result.begin(), result.capacity(), "%u,?", static_cast<uint32_t>(baudRate)); }
        else { std::snprintf(result.begin(), result.capacity(), "%u,%u", static_cast<uint32_t>(baudRate), static_cast<uint8_t>(serialPort)); }
    }
    return result;
}

Command BaudRate::toReadCommand()
{
    AsciiMessage commandString;
    if (serialPort == SerialPort::Poll) { std::snprintf(commandString.begin(), commandString.capacity(), "RRG,%02d,?", _id); }
    else { std::snprintf(commandString.begin(), commandString.capacity(), "RRG,%02d,%1d", _id, static_cast<uint8_t>(serialPort)); }

    Command readCommand(commandString, 6);
    return readCommand;
}

bool operator==(const BaudRate& lhs, const BaudRate& rhs) { return ((lhs.baudRate == rhs.baudRate) && (lhs.serialPort == rhs.serialPort)); }

bool AsyncOutputType::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 2 + 1)
    {
        // This register may not have an optional field
        if (tokens.size() != 1 + 1) { return true; }
    }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto ador_t = StringUtils::fromString<uint32_t>(start, end);
    if (!ador_t.has_value()) return true;
    ador = static_cast<AsyncOutputType::Ador>(ador_t.value());

    if (tokens.size() == (2 + 1))
    {
        // This is an optional parameter.
        start = tokens[index++] + 1;
        end = tokens[index];
        const auto serialPort_t = StringUtils::fromString<uint8_t>(start, end);
        if (!serialPort_t.has_value()) return true;
        serialPort = static_cast<AsyncOutputType::SerialPort>(serialPort_t.value());
    }

    return false;  // no errors occurred
}

AsciiMessage AsyncOutputType::toString() const
{
    AsciiMessage result = "";
    if (static_cast<uint8_t>(serialPort) == 0) { std::snprintf(result.begin(), result.capacity(), "%d", static_cast<uint32_t>(ador)); }
    else
    {
        if (serialPort == SerialPort::Poll) { std::snprintf(result.begin(), result.capacity(), "%u,?", static_cast<uint32_t>(ador)); }
        else { std::snprintf(result.begin(), result.capacity(), "%u,%u", static_cast<uint32_t>(ador), static_cast<uint8_t>(serialPort)); }
    }
    return result;
}

Command AsyncOutputType::toReadCommand()
{
    AsciiMessage commandString;
    if (serialPort == SerialPort::Poll) { std::snprintf(commandString.begin(), commandString.capacity(), "RRG,%02d,?", _id); }
    else { std::snprintf(commandString.begin(), commandString.capacity(), "RRG,%02d,%1d", _id, static_cast<uint8_t>(serialPort)); }

    Command readCommand(commandString, 6);
    return readCommand;
}

bool operator==(const AsyncOutputType& lhs, const AsyncOutputType& rhs) { return ((lhs.ador == rhs.ador) && (lhs.serialPort == rhs.serialPort)); }

bool AsyncOutputFreq::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 2 + 1)
    {
        // This register may not have an optional field
        if (tokens.size() != 1 + 1) { return true; }
    }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto adof_t = StringUtils::fromString<uint32_t>(start, end);
    if (!adof_t.has_value()) return true;
    adof = static_cast<AsyncOutputFreq::Adof>(adof_t.value());

    if (tokens.size() == (2 + 1))
    {
        // This is an optional parameter.
        start = tokens[index++] + 1;
        end = tokens[index];
        const auto serialPort_t = StringUtils::fromString<uint8_t>(start, end);
        if (!serialPort_t.has_value()) return true;
        serialPort = static_cast<AsyncOutputFreq::SerialPort>(serialPort_t.value());
    }

    return false;  // no errors occurred
}

AsciiMessage AsyncOutputFreq::toString() const
{
    AsciiMessage result = "";
    if (static_cast<uint8_t>(serialPort) == 0) { std::snprintf(result.begin(), result.capacity(), "%d", static_cast<uint32_t>(adof)); }
    else
    {
        if (serialPort == SerialPort::Poll) { std::snprintf(result.begin(), result.capacity(), "%u,?", static_cast<uint32_t>(adof)); }
        else { std::snprintf(result.begin(), result.capacity(), "%u,%u", static_cast<uint32_t>(adof), static_cast<uint8_t>(serialPort)); }
    }
    return result;
}

Command AsyncOutputFreq::toReadCommand()
{
    AsciiMessage commandString;
    if (serialPort == SerialPort::Poll) { std::snprintf(commandString.begin(), commandString.capacity(), "RRG,%02d,?", _id); }
    else { std::snprintf(commandString.begin(), commandString.capacity(), "RRG,%02d,%1d", _id, static_cast<uint8_t>(serialPort)); }

    Command readCommand(commandString, 6);
    return readCommand;
}

bool operator==(const AsyncOutputFreq& lhs, const AsyncOutputFreq& rhs) { return ((lhs.adof == rhs.adof) && (lhs.serialPort == rhs.serialPort)); }

bool ProtocolControl::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 7 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto asciiAppendCount_t = StringUtils::fromString<uint8_t>(start, end);
    if (!asciiAppendCount_t.has_value()) return true;
    asciiAppendCount = static_cast<ProtocolControl::AsciiAppendCount>(asciiAppendCount_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto asciiAppendStatus_t = StringUtils::fromString<uint8_t>(start, end);
    if (!asciiAppendStatus_t.has_value()) return true;
    asciiAppendStatus = static_cast<ProtocolControl::AsciiAppendStatus>(asciiAppendStatus_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto spiAppendCount_t = StringUtils::fromString<uint8_t>(start, end);
    if (!spiAppendCount_t.has_value()) return true;
    spiAppendCount = static_cast<ProtocolControl::SpiAppendCount>(spiAppendCount_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto spiAppendStatus_t = StringUtils::fromString<uint8_t>(start, end);
    if (!spiAppendStatus_t.has_value()) return true;
    spiAppendStatus = static_cast<ProtocolControl::SpiAppendStatus>(spiAppendStatus_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto asciiChecksum_t = StringUtils::fromString<uint8_t>(start, end);
    if (!asciiChecksum_t.has_value()) return true;
    asciiChecksum = static_cast<ProtocolControl::AsciiChecksum>(asciiChecksum_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto spiChecksum_t = StringUtils::fromString<uint8_t>(start, end);
    if (!spiChecksum_t.has_value()) return true;
    spiChecksum = static_cast<ProtocolControl::SpiChecksum>(spiChecksum_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto errorMode_t = StringUtils::fromString<uint8_t>(start, end);
    if (!errorMode_t.has_value()) return true;
    errorMode = static_cast<ProtocolControl::ErrorMode>(errorMode_t.value());

    return false;  // no errors occurred
}

AsciiMessage ProtocolControl::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%u,%u,%u,%u,%u,%u,%u", static_cast<uint8_t>(asciiAppendCount), static_cast<uint8_t>(asciiAppendStatus),
                  static_cast<uint8_t>(spiAppendCount), static_cast<uint8_t>(spiAppendStatus), static_cast<uint8_t>(asciiChecksum),
                  static_cast<uint8_t>(spiChecksum), static_cast<uint8_t>(errorMode));
    return result;
}

bool operator==(const ProtocolControl& lhs, const ProtocolControl& rhs)
{
    return ((lhs.asciiAppendCount == rhs.asciiAppendCount) && (lhs.asciiAppendStatus == rhs.asciiAppendStatus) && (lhs.spiAppendCount == rhs.spiAppendCount) &&
            (lhs.spiAppendStatus == rhs.spiAppendStatus) && (lhs.asciiChecksum == rhs.asciiChecksum) && (lhs.spiChecksum == rhs.spiChecksum) &&
            (lhs.errorMode == rhs.errorMode));
}

bool SyncControl::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 9 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto syncInMode_t = StringUtils::fromString<uint8_t>(start, end);
    if (!syncInMode_t.has_value()) return true;
    syncInMode = static_cast<SyncControl::SyncInMode>(syncInMode_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto syncInEdge_t = StringUtils::fromString<uint8_t>(start, end);
    if (!syncInEdge_t.has_value()) return true;
    syncInEdge = static_cast<SyncControl::SyncInEdge>(syncInEdge_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto syncInSkipFactor_t = StringUtils::fromString<uint16_t>(start, end);
    if (!syncInSkipFactor_t.has_value()) return true;
    syncInSkipFactor = syncInSkipFactor_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto resv1_t = StringUtils::fromString<uint32_t>(start, end);
    if (!resv1_t.has_value()) return true;
    resv1 = resv1_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto syncOutMode_t = StringUtils::fromString<uint8_t>(start, end);
    if (!syncOutMode_t.has_value()) return true;
    syncOutMode = static_cast<SyncControl::SyncOutMode>(syncOutMode_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto syncOutPolarity_t = StringUtils::fromString<uint8_t>(start, end);
    if (!syncOutPolarity_t.has_value()) return true;
    syncOutPolarity = static_cast<SyncControl::SyncOutPolarity>(syncOutPolarity_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto syncOutSkipFactor_t = StringUtils::fromString<uint16_t>(start, end);
    if (!syncOutSkipFactor_t.has_value()) return true;
    syncOutSkipFactor = syncOutSkipFactor_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto syncOutPulseWidth_t = StringUtils::fromString<uint32_t>(start, end);
    if (!syncOutPulseWidth_t.has_value()) return true;
    syncOutPulseWidth = syncOutPulseWidth_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto resv2_t = StringUtils::fromString<uint32_t>(start, end);
    if (!resv2_t.has_value()) return true;
    resv2 = resv2_t.value();

    return false;  // no errors occurred
}

AsciiMessage SyncControl::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%u,%u,%u,%u,%u,%u,%u,%u,%u", static_cast<uint8_t>(syncInMode), static_cast<uint8_t>(syncInEdge),
                  syncInSkipFactor, resv1, static_cast<uint8_t>(syncOutMode), static_cast<uint8_t>(syncOutPolarity), syncOutSkipFactor, syncOutPulseWidth,
                  resv2);
    return result;
}

bool operator==(const SyncControl& lhs, const SyncControl& rhs)
{
    return ((lhs.syncInMode == rhs.syncInMode) && (lhs.syncInEdge == rhs.syncInEdge) && (lhs.syncInSkipFactor == rhs.syncInSkipFactor) &&
            (lhs.resv1 == rhs.resv1) && (lhs.syncOutMode == rhs.syncOutMode) && (lhs.syncOutPolarity == rhs.syncOutPolarity) &&
            (lhs.syncOutSkipFactor == rhs.syncOutSkipFactor) && (lhs.syncOutPulseWidth == rhs.syncOutPulseWidth) && (lhs.resv2 == rhs.resv2));
}

bool SyncStatus::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 3 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto syncInCount_t = StringUtils::fromString<uint32_t>(start, end);
    if (!syncInCount_t.has_value()) return true;
    syncInCount = syncInCount_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto syncInTime_t = StringUtils::fromString<uint32_t>(start, end);
    if (!syncInTime_t.has_value()) return true;
    syncInTime = syncInTime_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto syncOutCount_t = StringUtils::fromString<uint32_t>(start, end);
    if (!syncOutCount_t.has_value()) return true;
    syncOutCount = syncOutCount_t.value();

    return false;  // no errors occurred
}

bool operator==(const SyncStatus& lhs, const SyncStatus& rhs)
{
    return ((lhs.syncInCount == rhs.syncInCount) && (lhs.syncInTime == rhs.syncInTime) && (lhs.syncOutCount == rhs.syncOutCount));
}

bool BinaryOutput::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    auto asyncMode_t = StringUtils::fromStringHex<uint16_t>(start, end);
    if (!asyncMode_t.has_value()) return true;
    asyncMode._value = asyncMode_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    auto rateDivisor_t = StringUtils::fromString<uint16_t>(start, end);
    if (!rateDivisor_t.has_value()) return true;
    rateDivisor = rateDivisor_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    auto group_t = StringUtils::fromStringHex<uint16_t>(start, end);
    if (!group_t.has_value()) return true;
    uint8_t group = group_t.value();
    uint8_t groupExt = 0;
    if (groupExt & 0x80)
    {
        start = tokens[index++] + 1;
        end = tokens[index];
        auto groupExt_t = StringUtils::fromStringHex<uint16_t>(start, end);
        if (!groupExt_t.has_value()) return true;
        groupExt = groupExt_t.value();
    }

    if (group & COMMON_BIT)
    {
        start = tokens[index++] + 1;
        end = tokens[index];
        auto common_t = StringUtils::fromStringHex<uint16_t>(start, end);
        if (!common_t.has_value()) return true;
        this->common._value = common_t.value();
        if (this->common._value & 0x8000)
        {
            start = tokens[index++] + 1;
            end = tokens[index];
            this->common._value &= 0x3FFF;
            auto common_t = StringUtils::fromStringHex<uint16_t>(start, end);
            if (!common_t.has_value()) return true;
            this->common._value |= common_t.value() << 7;
        }
    }

    if (group & TIME_BIT)
    {
        start = tokens[index++] + 1;
        end = tokens[index];
        auto time_t = StringUtils::fromStringHex<uint16_t>(start, end);
        if (!time_t.has_value()) return true;
        this->time._value = time_t.value();
        if (this->time._value & 0x8000)
        {
            start = tokens[index++] + 1;
            end = tokens[index];
            this->time._value &= 0x3FFF;
            auto time_t = StringUtils::fromStringHex<uint16_t>(start, end);
            if (!time_t.has_value()) return true;
            this->time._value |= time_t.value() << 7;
        }
    }

    if (group & IMU_BIT)
    {
        start = tokens[index++] + 1;
        end = tokens[index];
        auto imu_t = StringUtils::fromStringHex<uint16_t>(start, end);
        if (!imu_t.has_value()) return true;
        this->imu._value = imu_t.value();
        if (this->imu._value & 0x8000)
        {
            start = tokens[index++] + 1;
            end = tokens[index];
            this->imu._value &= 0x3FFF;
            auto imu_t = StringUtils::fromStringHex<uint16_t>(start, end);
            if (!imu_t.has_value()) return true;
            this->imu._value |= imu_t.value() << 7;
        }
    }

    if (group & GNSS_BIT)
    {
        start = tokens[index++] + 1;
        end = tokens[index];
        auto gnss_t = StringUtils::fromStringHex<uint16_t>(start, end);
        if (!gnss_t.has_value()) return true;
        this->gnss._value = gnss_t.value();
        if (this->gnss._value & 0x8000)
        {
            start = tokens[index++] + 1;
            end = tokens[index];
            this->gnss._value &= 0x3FFF;
            auto gnss_t = StringUtils::fromStringHex<uint16_t>(start, end);
            if (!gnss_t.has_value()) return true;
            this->gnss._value |= gnss_t.value() << 7;
        }
    }

    if (group & ATTITUDE_BIT)
    {
        start = tokens[index++] + 1;
        end = tokens[index];
        auto attitude_t = StringUtils::fromStringHex<uint16_t>(start, end);
        if (!attitude_t.has_value()) return true;
        this->attitude._value = attitude_t.value();
        if (this->attitude._value & 0x8000)
        {
            start = tokens[index++] + 1;
            end = tokens[index];
            this->attitude._value &= 0x3FFF;
            auto attitude_t = StringUtils::fromStringHex<uint16_t>(start, end);
            if (!attitude_t.has_value()) return true;
            this->attitude._value |= attitude_t.value() << 7;
        }
    }

    if (group & INS_BIT)
    {
        start = tokens[index++] + 1;
        end = tokens[index];
        auto ins_t = StringUtils::fromStringHex<uint16_t>(start, end);
        if (!ins_t.has_value()) return true;
        this->ins._value = ins_t.value();
        if (this->ins._value & 0x8000)
        {
            start = tokens[index++] + 1;
            end = tokens[index];
            this->ins._value &= 0x3FFF;
            auto ins_t = StringUtils::fromStringHex<uint16_t>(start, end);
            if (!ins_t.has_value()) return true;
            this->ins._value |= ins_t.value() << 7;
        }
    }

    if (group & GNSS2_BIT)
    {
        start = tokens[index++] + 1;
        end = tokens[index];
        auto gnss2_t = StringUtils::fromStringHex<uint16_t>(start, end);
        if (!gnss2_t.has_value()) return true;
        this->gnss2._value = gnss2_t.value();
        if (this->gnss2._value & 0x8000)
        {
            start = tokens[index++] + 1;
            end = tokens[index];
            this->gnss2._value &= 0x3FFF;
            auto gnss2_t = StringUtils::fromStringHex<uint16_t>(start, end);
            if (!gnss2_t.has_value()) return true;
            this->gnss2._value |= gnss2_t.value() << 7;
        }
    }

    if (groupExt & (GNSS3_BIT >> 8))
    {
        start = tokens[index++] + 1;
        end = tokens[index];
        auto gnss3_t = StringUtils::fromStringHex<uint16_t>(start, end);
        if (!gnss3_t.has_value()) return true;
        this->gnss3._value = gnss3_t.value();
        if (this->gnss3._value & 0x8000)
        {
            start = tokens[index++] + 1;
            end = tokens[index];
            this->gnss3._value &= 0x3FFF;
            auto gnss3_t = StringUtils::fromStringHex<uint16_t>(start, end);
            if (!gnss3_t.has_value()) return true;
            this->gnss3._value |= gnss3_t.value() << 7;
        }
    }

    return false;  // no errors occurred
}

AsciiMessage BinaryOutput::toString() const
{
    AsciiMessage result = "";
    BinaryHeader binaryHeader = toBinaryHeader();
    std::snprintf(result.end(), result.numAvailable(), "%1X,%d%s", asyncMode._value, rateDivisor,
                  binaryHeaderToString<AsciiMessage::capacity()>(binaryHeader).c_str());
    return result;
}

BinaryHeader BinaryOutputMeasurements::toBinaryHeader() const noexcept
{
    BinaryHeader result;
    uint8_t group = 0;
    uint8_t groupExt = 0;
    if (this->common > 0) group += COMMON_BIT;
    if (this->time > 0) group += TIME_BIT;
    if (this->imu > 0) group += IMU_BIT;
    if (this->gnss > 0) group += GNSS_BIT;
    if (this->attitude > 0) group += ATTITUDE_BIT;
    if (this->ins > 0) group += INS_BIT;
    if (this->gnss2 > 0) group += GNSS2_BIT;
    if (this->gnss3 > 0) groupExt += (GNSS3_BIT >> 8);
    if (this->gnss3 > 0) group |= 0x80;
    result.outputGroups.push_back(group);
    if (groupExt > 0) result.outputGroups.push_back(groupExt);

    if (this->common > 0)
    {
        result.outputTypes.push_back(this->common._value);
        if ((this->common & 0xFFFF0000) > 0)
        {
            result.outputTypes.back() |= 0x8000;  // Set extension bit
            result.outputTypes.push_back((this->common._value & 0xFFFF0000) >> 16);
        }
    }

    if (this->time > 0)
    {
        result.outputTypes.push_back(this->time._value);
        if ((this->time & 0xFFFF0000) > 0)
        {
            result.outputTypes.back() |= 0x8000;  // Set extension bit
            result.outputTypes.push_back((this->time._value & 0xFFFF0000) >> 16);
        }
    }

    if (this->imu > 0)
    {
        result.outputTypes.push_back(this->imu._value);
        if ((this->imu & 0xFFFF0000) > 0)
        {
            result.outputTypes.back() |= 0x8000;  // Set extension bit
            result.outputTypes.push_back((this->imu._value & 0xFFFF0000) >> 16);
        }
    }

    if (this->gnss > 0)
    {
        result.outputTypes.push_back(this->gnss._value);
        if ((this->gnss & 0xFFFF0000) > 0)
        {
            result.outputTypes.back() |= 0x8000;  // Set extension bit
            result.outputTypes.push_back((this->gnss._value & 0xFFFF0000) >> 16);
        }
    }

    if (this->attitude > 0)
    {
        result.outputTypes.push_back(this->attitude._value);
        if ((this->attitude & 0xFFFF0000) > 0)
        {
            result.outputTypes.back() |= 0x8000;  // Set extension bit
            result.outputTypes.push_back((this->attitude._value & 0xFFFF0000) >> 16);
        }
    }

    if (this->ins > 0)
    {
        result.outputTypes.push_back(this->ins._value);
        if ((this->ins & 0xFFFF0000) > 0)
        {
            result.outputTypes.back() |= 0x8000;  // Set extension bit
            result.outputTypes.push_back((this->ins._value & 0xFFFF0000) >> 16);
        }
    }

    if (this->gnss2 > 0)
    {
        result.outputTypes.push_back(this->gnss2._value);
        if ((this->gnss2 & 0xFFFF0000) > 0)
        {
            result.outputTypes.back() |= 0x8000;  // Set extension bit
            result.outputTypes.push_back((this->gnss2._value & 0xFFFF0000) >> 16);
        }
    }

    if (this->gnss3 > 0)
    {
        result.outputTypes.push_back(this->gnss3._value);
        if ((this->gnss3 & 0xFFFF0000) > 0)
        {
            result.outputTypes.back() |= 0x8000;  // Set extension bit
            result.outputTypes.push_back((this->gnss3._value & 0xFFFF0000) >> 16);
        }
    }

    return result;
}

bool NmeaOutput1::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 5 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto port_t = StringUtils::fromString<uint8_t>(start, end);
    if (!port_t.has_value()) return true;
    port = static_cast<NmeaOutput1::Port>(port_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto rate_t = StringUtils::fromString<uint8_t>(start, end);
    if (!rate_t.has_value()) return true;
    rate = static_cast<NmeaOutput1::Rate>(rate_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto mode_t = StringUtils::fromString<uint8_t>(start, end);
    if (!mode_t.has_value()) return true;
    mode = static_cast<NmeaOutput1::Mode>(mode_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnssSelect_t = StringUtils::fromString<uint8_t>(start, end);
    if (!gnssSelect_t.has_value()) return true;
    gnssSelect = static_cast<NmeaOutput1::GnssSelect>(gnssSelect_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto msgSelection_t = StringUtils::fromStringHex<uint32_t>(start, end);
    if (!msgSelection_t.has_value()) return true;
    msgSelection._value = msgSelection_t.value();

    return false;  // no errors occurred
}

AsciiMessage NmeaOutput1::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%u,%u,%u,%u,%X", static_cast<uint8_t>(port), static_cast<uint8_t>(rate), static_cast<uint8_t>(mode),
                  static_cast<uint8_t>(gnssSelect), msgSelection._value);
    return result;
}

bool operator==(const NmeaOutput1& lhs, const NmeaOutput1& rhs)
{
    return ((lhs.port == rhs.port) && (lhs.rate == rhs.rate) && (lhs.mode == rhs.mode) && (lhs.gnssSelect == rhs.gnssSelect) &&
            (lhs.msgSelection._value == rhs.msgSelection._value));
}

bool NmeaOutput2::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 5 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto port_t = StringUtils::fromString<uint8_t>(start, end);
    if (!port_t.has_value()) return true;
    port = static_cast<NmeaOutput2::Port>(port_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto rate_t = StringUtils::fromString<uint8_t>(start, end);
    if (!rate_t.has_value()) return true;
    rate = static_cast<NmeaOutput2::Rate>(rate_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto mode_t = StringUtils::fromString<uint8_t>(start, end);
    if (!mode_t.has_value()) return true;
    mode = static_cast<NmeaOutput2::Mode>(mode_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnssSelect_t = StringUtils::fromString<uint8_t>(start, end);
    if (!gnssSelect_t.has_value()) return true;
    gnssSelect = static_cast<NmeaOutput2::GnssSelect>(gnssSelect_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto msgSelection_t = StringUtils::fromStringHex<uint32_t>(start, end);
    if (!msgSelection_t.has_value()) return true;
    msgSelection._value = msgSelection_t.value();

    return false;  // no errors occurred
}

AsciiMessage NmeaOutput2::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%u,%u,%u,%u,%X", static_cast<uint8_t>(port), static_cast<uint8_t>(rate), static_cast<uint8_t>(mode),
                  static_cast<uint8_t>(gnssSelect), msgSelection._value);
    return result;
}

bool operator==(const NmeaOutput2& lhs, const NmeaOutput2& rhs)
{
    return ((lhs.port == rhs.port) && (lhs.rate == rhs.rate) && (lhs.mode == rhs.mode) && (lhs.gnssSelect == rhs.gnssSelect) &&
            (lhs.msgSelection._value == rhs.msgSelection._value));
}

bool LegacyCompatibilitySettings::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 4 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto insLegacy_t = StringUtils::fromString<uint8_t>(start, end);
    if (!insLegacy_t.has_value()) return true;
    insLegacy = static_cast<LegacyCompatibilitySettings::InsLegacy>(insLegacy_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto gnssLegacy_t = StringUtils::fromString<uint8_t>(start, end);
    if (!gnssLegacy_t.has_value()) return true;
    gnssLegacy._value = gnssLegacy_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto imuLegacy_t = StringUtils::fromString<uint8_t>(start, end);
    if (!imuLegacy_t.has_value()) return true;
    imuLegacy = static_cast<LegacyCompatibilitySettings::ImuLegacy>(imuLegacy_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto hwLegacy_t = StringUtils::fromString<uint8_t>(start, end);
    if (!hwLegacy_t.has_value()) return true;
    hwLegacy = static_cast<LegacyCompatibilitySettings::HwLegacy>(hwLegacy_t.value());

    return false;  // no errors occurred
}

AsciiMessage LegacyCompatibilitySettings::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%u,%u,%u,%u", static_cast<uint8_t>(insLegacy), gnssLegacy._value, static_cast<uint8_t>(imuLegacy),
                  static_cast<uint8_t>(hwLegacy));
    return result;
}

bool operator==(const LegacyCompatibilitySettings& lhs, const LegacyCompatibilitySettings& rhs)
{
    return ((lhs.insLegacy == rhs.insLegacy) && (lhs.gnssLegacy._value == rhs.gnssLegacy._value) && (lhs.imuLegacy == rhs.imuLegacy) &&
            (lhs.hwLegacy == rhs.hwLegacy));
}

}  // namespace System

namespace VelocityAiding
{
bool VelAidingMeas::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 3 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velocityX_t = StringUtils::fromString<float>(start, end);
    if (!velocityX_t.has_value()) return true;
    velocityX = velocityX_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velocityY_t = StringUtils::fromString<float>(start, end);
    if (!velocityY_t.has_value()) return true;
    velocityY = velocityY_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velocityZ_t = StringUtils::fromString<float>(start, end);
    if (!velocityZ_t.has_value()) return true;
    velocityZ = velocityZ_t.value();

    return false;  // no errors occurred
}

AsciiMessage VelAidingMeas::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%f,%f,%f", velocityX, velocityY, velocityZ);
    return result;
}

bool operator==(const VelAidingMeas& lhs, const VelAidingMeas& rhs)
{
    return ((lhs.velocityX == rhs.velocityX) && (lhs.velocityY == rhs.velocityY) && (lhs.velocityZ == rhs.velocityZ));
}

bool VelAidingControl::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 3 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velAidEnable_t = StringUtils::fromString<uint8_t>(start, end);
    if (!velAidEnable_t.has_value()) return true;
    velAidEnable = static_cast<VelAidingControl::VelAidEnable>(velAidEnable_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto velUncertTuning_t = StringUtils::fromString<float>(start, end);
    if (!velUncertTuning_t.has_value()) return true;
    velUncertTuning = velUncertTuning_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto resv_t = StringUtils::fromString<float>(start, end);
    if (!resv_t.has_value()) return true;
    resv = resv_t.value();

    return false;  // no errors occurred
}

AsciiMessage VelAidingControl::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%u,%f,%f", static_cast<uint8_t>(velAidEnable), velUncertTuning, resv);
    return result;
}

bool operator==(const VelAidingControl& lhs, const VelAidingControl& rhs)
{
    return ((lhs.velAidEnable == rhs.velAidEnable) && (lhs.velUncertTuning == rhs.velUncertTuning) && (lhs.resv == rhs.resv));
}

}  // namespace VelocityAiding

namespace WorldMagGravityModel
{
bool RefModelConfig::fromString(const AsciiMessage& response)
{
    const auto tokens = findIndexOfFieldSeparators(response);
    if (tokens.size() != 9 + 1) { return true; }

    int index = 0;
    const char* start;
    const char* end;

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto enableMagModel_t = StringUtils::fromString<uint8_t>(start, end);
    if (!enableMagModel_t.has_value()) return true;
    enableMagModel = static_cast<RefModelConfig::EnableMagModel>(enableMagModel_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto enableGravityModel_t = StringUtils::fromString<uint8_t>(start, end);
    if (!enableGravityModel_t.has_value()) return true;
    enableGravityModel = static_cast<RefModelConfig::EnableGravityModel>(enableGravityModel_t.value());

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto resv1_t = StringUtils::fromString<uint8_t>(start, end);
    if (!resv1_t.has_value()) return true;
    resv1 = resv1_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto resv2_t = StringUtils::fromString<uint8_t>(start, end);
    if (!resv2_t.has_value()) return true;
    resv2 = resv2_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto recalcThreshold_t = StringUtils::fromString<uint32_t>(start, end);
    if (!recalcThreshold_t.has_value()) return true;
    recalcThreshold = recalcThreshold_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto year_t = StringUtils::fromString<float>(start, end);
    if (!year_t.has_value()) return true;
    year = year_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto latitude_t = StringUtils::fromString<double>(start, end);
    if (!latitude_t.has_value()) return true;
    latitude = latitude_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto longitude_t = StringUtils::fromString<double>(start, end);
    if (!longitude_t.has_value()) return true;
    longitude = longitude_t.value();

    start = tokens[index++] + 1;
    end = tokens[index];
    const auto altitude_t = StringUtils::fromString<double>(start, end);
    if (!altitude_t.has_value()) return true;
    altitude = altitude_t.value();

    return false;  // no errors occurred
}

AsciiMessage RefModelConfig::toString() const
{
    AsciiMessage result = "";
    std::snprintf(result.begin(), result.capacity(), "%u,%u,%u,%u,%u,%f,%f,%f,%f", static_cast<uint8_t>(enableMagModel),
                  static_cast<uint8_t>(enableGravityModel), resv1, resv2, recalcThreshold, year, latitude, longitude, altitude);
    return result;
}

bool operator==(const RefModelConfig& lhs, const RefModelConfig& rhs)
{
    return ((lhs.enableMagModel == rhs.enableMagModel) && (lhs.enableGravityModel == rhs.enableGravityModel) && (lhs.resv1 == rhs.resv1) &&
            (lhs.resv2 == rhs.resv2) && (lhs.recalcThreshold == rhs.recalcThreshold) && (lhs.year == rhs.year) && (lhs.latitude == rhs.latitude) &&
            (lhs.longitude == rhs.longitude) && (lhs.altitude == rhs.altitude));
}

}  // namespace WorldMagGravityModel

Vector<const char*, Config::PacketFinders::asciiMaxFieldCount> findIndexOfFieldSeparators(const AsciiMessage& input)
{
    Vector<const char*, Config::PacketFinders::asciiMaxFieldCount> indexes;
    size_t offset = 7;
    if (input[0] != '$')
    {
        indexes.push_back(input.begin() - 1);
        offset = 0;
    }

    for (size_t i = input.find(',', offset); i < input.length(); i++)
    {
        if (input[i] == ',' || input[i] == '*') { indexes.push_back(&input[i]); }
    }

    if (input[0] != '$') { indexes.push_back(input.end()); }

    return indexes;
}

}  // namespace Registers
}  // namespace VN
