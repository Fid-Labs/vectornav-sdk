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

#ifndef __REGISTERS_HPP__
#define __REGISTERS_HPP__

#include <array>
#include "Implementation/BinaryHeader.hpp"
#include "Implementation/BinaryMeasurementDefinitions.hpp"
#include "Interface/Commands.hpp"
#include "TemplateLibrary/Vector.hpp"

namespace VN
{
namespace Registers
{
namespace Attitude
{
/**--------------------------------------------------------------------------------------------------
<summary> Register 8 - Yaw Pitch Roll </summary>
<remarks>

Attitude solution as yaw, pitch, and roll in degrees. The yaw, pitch, and roll is given as a 3,2,1
Euler angle rotation sequence describing the orientation of the sensor with respect to the inertial
North East Down (NED) frame.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class YawPitchRoll : public MeasurementRegister
{
public:
    YawPitchRoll() : MeasurementRegister(8) {}
    std::string name() override { return "YawPitchRoll"; };

    float yaw = 0;    ///< Yaw angle.
    float pitch = 0;  ///< Pitch angle.
    float roll = 0;   ///< Roll angle.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const YawPitchRoll& lhs, const YawPitchRoll& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 9 - Quaternion </summary>
<remarks>

Attitude solution as a quaternion.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class Quaternion : public MeasurementRegister
{
public:
    Quaternion() : MeasurementRegister(9) {}
    std::string name() override { return "Quaternion"; };

    float quatX = 0;  ///< First vector component of quaternion.
    float quatY = 0;  ///< Second vector component of quaternion.
    float quatZ = 0;  ///< Third vector component of quaternion.
    float quatS = 0;  ///< Scalar component of quaternion.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const Quaternion& lhs, const Quaternion& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 15 - Quaternion & Compensated IMU </summary>
<remarks>

Quaternion attitude solution, and compensated (Magnetic, Acceleration, Angular Rate) values.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class QuatMagAccelRate : public MeasurementRegister
{
public:
    QuatMagAccelRate() : MeasurementRegister(15) {}
    std::string name() override { return "QuatMagAccelRate"; };

    float quatX = 0;   ///< First vector component of quaternion.
    float quatY = 0;   ///< Second vector component of quaternion.
    float quatZ = 0;   ///< Third vector component of quaternion.
    float quatS = 0;   ///< Scalar component of quaternion.
    float magX = 0;    ///< Compensated magnetometer measurement in the body-frame x-axis.
    float magY = 0;    ///< Compensated magnetometer measurement in the body-frame y-axis.
    float magZ = 0;    ///< Compensated magnetometer measurement in the body-frame z-axis.
    float accelX = 0;  ///< Compensated accelerometer measurement in the body-frame x-axis.
    float accelY = 0;  ///< Compensated accelerometer measurement in the body-frame y-axis.
    float accelZ = 0;  ///< Compensated accelerometer measurement in the body-frame z-axis.
    float gyroX = 0;   ///< Compensated angular rate measurement in the body-frame x-axis.
    float gyroY = 0;   ///< Compensated angular rate measurement in the body-frame y-axis.
    float gyroZ = 0;   ///< Compensated angular rate measurement in the body-frame z-axis.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const QuatMagAccelRate& lhs, const QuatMagAccelRate& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 21 - Magnetic and Gravity Reference Vectors </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class MagGravRefVec : public ConfigurationRegister
{
public:
    MagGravRefVec() : ConfigurationRegister(21) {}
    std::string name() override { return "MagGravRefVec"; };

    float magRefN = 0;
    float magRefE = 0;
    float magRefD = 0;
    float gravRefN = 0;
    float gravRefE = 0;
    float gravRefD = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const MagGravRefVec& lhs, const MagGravRefVec& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 27 - Yaw-Pitch-Roll & Compensated IMU </summary>
<remarks>

Yaw, Pitch, Roll, Accel, and Angular Rates

</remarks>
*-----------------------------------------------------------------------------------------------**/
class YprMagAccelAngularRates : public MeasurementRegister
{
public:
    YprMagAccelAngularRates() : MeasurementRegister(27) {}
    std::string name() override { return "YprMagAccelAngularRates"; };

    float yaw = 0;     ///< Yaw angle.
    float pitch = 0;   ///< Pitch angle.
    float roll = 0;    ///< Roll angle.
    float magX = 0;    ///< Compensated magnetometer measurement in the body-frame x-axis.
    float magY = 0;    ///< Compensated magnetometer measurement in the body-frame y-axis.
    float magZ = 0;    ///< Compensated magnetometer measurement in the body-frame z-axis.
    float accelX = 0;  ///< Compensated accelerometer measurement in the body-frame x-axis.
    float accelY = 0;  ///< Compensated accelerometer measurement in the body-frame y-axis.
    float accelZ = 0;  ///< Compensated accelerometer measurement in the body-frame z-axis.
    float gyroX = 0;   ///< Compensated angular rate measurement in the body-frame x-axis.
    float gyroY = 0;   ///< Compensated angular rate measurement in the body-frame y-axis.
    float gyroZ = 0;   ///< Compensated angular rate measurement in the body-frame z-axis.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const YprMagAccelAngularRates& lhs, const YprMagAccelAngularRates& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 35 - VPE Basic Control </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class VpeBasicControl : public ConfigurationRegister
{
public:
    VpeBasicControl() : ConfigurationRegister(35) {}
    std::string name() override { return "VpeBasicControl"; };

    uint8_t resv = 0;
    enum class HeadingMode : uint8_t
    {
        Absolute = 0,
        Relative = 1,
        Indoor = 2,
    };
    HeadingMode headingMode = static_cast<HeadingMode>(0);
    enum class FilteringMode : uint8_t
    {
        Unfiltered = 0,
        AdaptivelyFiltered = 1,
    };
    FilteringMode filteringMode = static_cast<FilteringMode>(0);
    enum class TuningMode : uint8_t
    {
        Static = 0,
        Adaptive = 1,
    };
    TuningMode tuningMode = static_cast<TuningMode>(0);

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const VpeBasicControl& lhs, const VpeBasicControl& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 36 - VPE Magnetometer Basic Tuning </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class VpeMagBasicTuning : public ConfigurationRegister
{
public:
    VpeMagBasicTuning() : ConfigurationRegister(36) {}
    std::string name() override { return "VpeMagBasicTuning"; };

    float baseTuningX = 0;
    float baseTuningY = 0;
    float baseTuningZ = 0;
    float adaptiveTuningX = 0;
    float adaptiveTuningY = 0;
    float adaptiveTuningZ = 0;
    float adaptiveFilteringX = 0;
    float adaptiveFilteringY = 0;
    float adaptiveFilteringZ = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const VpeMagBasicTuning& lhs, const VpeMagBasicTuning& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 38 - VPE Accelerometer Basic Tuning </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class VpeAccelBasicTuning : public ConfigurationRegister
{
public:
    VpeAccelBasicTuning() : ConfigurationRegister(38) {}
    std::string name() override { return "VpeAccelBasicTuning"; };

    float baseTuningX = 0;
    float baseTuningY = 0;
    float baseTuningZ = 0;
    float adaptiveTuningX = 0;
    float adaptiveTuningY = 0;
    float adaptiveTuningZ = 0;
    float adaptiveFilteringX = 0;
    float adaptiveFilteringY = 0;
    float adaptiveFilteringZ = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const VpeAccelBasicTuning& lhs, const VpeAccelBasicTuning& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 239 - Yaw-Pitch-Roll, Linear Acceleration & Gyro </summary>
<remarks>

Yaw, Pitch, Roll, Linear Body Accel, and Angular Rates.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class YprLinearBodyAccelAngularRates : public MeasurementRegister
{
public:
    YprLinearBodyAccelAngularRates() : MeasurementRegister(239) {}
    std::string name() override { return "YprLinearBodyAccelAngularRates"; };

    float yaw = 0;        ///< Yaw angle.
    float pitch = 0;      ///< Pitch angle.
    float roll = 0;       ///< Roll angle.
    float linAccelX = 0;  ///< Linear acceleration in body-frame x-axis.
    float linAccelY = 0;  ///< Linear acceleration in body-frame y-axis.
    float linAccelZ = 0;  ///< Linear acceleration in body-frame z-axis.
    float gyroX = 0;      ///< Compensated angular rate measurement in the body-frame x-axis.
    float gyroY = 0;      ///< Compensated angular rate measurement in the body-frame y-axis.
    float gyroZ = 0;      ///< Compensated angular rate measurement in the body-frame z-axis.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const YprLinearBodyAccelAngularRates& lhs, const YprLinearBodyAccelAngularRates& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 240 - Yaw-Pitch-Roll, Inertial Linear Acceleration & Gyro </summary>
<remarks>

Yaw, Pitch, Roll, Linear Inertial Accel, and Angular Rates.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class YprLinearInertialAccelAngularRates : public MeasurementRegister
{
public:
    YprLinearInertialAccelAngularRates() : MeasurementRegister(240) {}
    std::string name() override { return "YprLinearInertialAccelAngularRates"; };

    float yaw = 0;        ///< Yaw angle.
    float pitch = 0;      ///< Pitch angle.
    float roll = 0;       ///< Roll angle.
    float linAccelN = 0;  ///< Linear acceleration in North direction.
    float linAccelE = 0;  ///< Linear acceleration in East direction.
    float linAccelD = 0;  ///< Linear acceleration in Down direction.
    float gyroX = 0;      ///< Compensated angular rate measurement in the body-frame x-axis.
    float gyroY = 0;      ///< Compensated angular rate measurement in the body-frame y-axis.
    float gyroZ = 0;      ///< Compensated angular rate measurement in the body-frame z-axis.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const YprLinearInertialAccelAngularRates& lhs, const YprLinearInertialAccelAngularRates& rhs);
}  // namespace Attitude

namespace GNSS
{
/**--------------------------------------------------------------------------------------------------
<summary> Register 55 - GNSS Basic Configuration </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class GnssBasicConfig : public ConfigurationRegister
{
public:
    GnssBasicConfig() : ConfigurationRegister(55) {}
    std::string name() override { return "GnssBasicConfig"; };

    enum class ReceiverEnable : uint8_t
    {
        Internal = 0,
        VnWrite = 1,
        VnAdor = 2,
        GnssA = 100,
        VnWriteAndGnssA = 101,
        VnAdorAndGnssA = 102,
        VnWriteAndGnssAB = 201,
        VnAdorAndGnssAB = 202,
    };
    ReceiverEnable receiverEnable = static_cast<ReceiverEnable>(0);
    enum class PpsSource : uint8_t
    {
        PpsPinRising = 0,
        PpsPinFalling = 1,
        SyncInRising = 2,
        SyncInFalling = 3,
        None = 6,
    };
    PpsSource ppsSource = static_cast<PpsSource>(0);
    enum class Rate : uint8_t
    {
        Rate1Hz = 1,
        Rate5Hz = 5,
    };
    Rate rate = static_cast<Rate>(0);
    uint8_t resv4 = 0;
    enum class AntPower : uint8_t
    {
        Off = 0,
        Internal = 1,
        External = 2,
    };
    AntPower antPower = static_cast<AntPower>(0);

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const GnssBasicConfig& lhs, const GnssBasicConfig& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 57 - GNSS Internal A Antenna Offset </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class GnssAOffset : public ConfigurationRegister
{
public:
    GnssAOffset() : ConfigurationRegister(57) {}
    std::string name() override { return "GnssAOffset"; };

    float positionX = 0;
    float positionY = 0;
    float positionZ = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const GnssAOffset& lhs, const GnssAOffset& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 58 - GNSS Solution - LLA </summary>
<remarks>

Primary GNSS receiver measurement with lat/lon/alt position and velocity in NED frame.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class GnssSolLla : public MeasurementRegister
{
public:
    GnssSolLla() : MeasurementRegister(58) {}
    std::string name() override { return "GnssSolLla"; };

    double gpsTow = 0;      ///< GPS time of week.
    uint16_t gps1Week = 0;  ///< The current GPS week.
    enum class Gnss1Fix : uint8_t
    {
        NoFix = 0,     ///< Not tracking.
        TimeFix = 1,   ///< Time-only fix.
        Fix2D = 2,     ///< Time plus 2D position fix.
        Fix3D = 3,     ///< Time plus full 3D position fix.
        SBAS = 4,      ///< 3D Fix using SBAS corrections.
        RtkFloat = 7,  ///< 3D Fix using an approximate RTK solution.
        RtkFix = 8,    ///< 3D Fix using a fixed RTK solution.
    };
    Gnss1Fix gnss1Fix = static_cast<Gnss1Fix>(0);
    uint8_t gnss1NumSats = 0;        ///< Number of satellites tracked by GNSS receiver.
    double lat = 0;                  ///< GNSS geodetic latitude.
    double lon = 0;                  ///< GNSS longitude.
    double alt = 0;                  ///< GNSS altitude above WGS84 ellipsoid.
    float velN = 0;                  ///< GNSS velocity in North direction.
    float velE = 0;                  ///< GNSS velocity in East direction.
    float velD = 0;                  ///< GNSS velocity in Down direction.
    float posUncertaintyN = 0;       ///< GNSS position uncertainty, North direction.
    float posUncertaintyE = 0;       ///< GNSS position uncertainty, East direction.
    float posUncertaintyD = 0;       ///< GNSS position uncertainty, Down direction.
    float gnss1VelUncertainty = 0;   ///< GNSS velocity uncertainty (scalar).
    float gnss1TimeUncertainty = 0;  ///< GNSS time uncertainty.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const GnssSolLla& lhs, const GnssSolLla& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 59 - GNSS Solution - ECEF </summary>
<remarks>

Primary GNSS receiver measurement in ECEF frame.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class GnssSolEcef : public MeasurementRegister
{
public:
    GnssSolEcef() : MeasurementRegister(59) {}
    std::string name() override { return "GnssSolEcef"; };

    double gpsTow = 0;      ///< GPS time of week.
    uint16_t gps1Week = 0;  ///< The current GPS week.
    enum class Gnss1Fix : uint8_t
    {
        NoFix = 0,     ///< Not tracking.
        TimeFix = 1,   ///< Time-only fix.
        Fix2D = 2,     ///< Time plus 2D position fix.
        Fix3D = 3,     ///< Time plus full 3D position fix.
        SBAS = 4,      ///< 3D Fix using SBAS corrections.
        RtkFloat = 7,  ///< 3D Fix using an approximate RTK solution.
        RtkFix = 8,    ///< 3D Fix using a fixed RTK solution.
    };
    Gnss1Fix gnss1Fix = static_cast<Gnss1Fix>(0);
    uint8_t gnss1NumSats = 0;        ///< Number of satellites tracked by GNSS receiver.
    double posX = 0;                 ///< GNSS position in ECEF-frame x-axis.
    double posY = 0;                 ///< GNSS position in ECEF-frame y-axis.
    double posZ = 0;                 ///< GNSS position in ECEF-frame z-axis.
    float velX = 0;                  ///< GNSS velocity in ECEF-frame x-axis.
    float velY = 0;                  ///< GNSS velocity in ECEF-frame y-axis.
    float velZ = 0;                  ///< GNSS velocity in ECEF-frame z-axis.
    float posUncertaintyX = 0;       ///< GNSS position uncertainty ECEF X.
    float posUncertaintyY = 0;       ///< GNSS position uncertainty ECEF Y.
    float posUncertaintyZ = 0;       ///< GNSS position uncertainty ECEF Z.
    float gnss1VelUncertainty = 0;   ///< GNSS velocity uncertainty (scalar).
    float gnss1TimeUncertainty = 0;  ///< GNSS time uncertainty.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const GnssSolEcef& lhs, const GnssSolEcef& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 99 - GNSS System Configuration </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class GnssSystemConfig : public ConfigurationRegister
{
public:
    GnssSystemConfig() : ConfigurationRegister(99) {}
    std::string name() override { return "GnssSystemConfig"; };

    union Systems
    {
        uint16_t _value;
        struct
        {
            uint16_t gps : 1;
            uint16_t sbas : 1;
            uint16_t glonass : 1;
            uint16_t beidou : 1;
            uint16_t galileo : 1;
            uint16_t imes : 1;
            uint16_t qzssL1Ca : 1;
            uint16_t qzssL1Saif : 1;
        };
        Systems(const uint16_t& other) : _value(other) {}
        operator uint16_t() const { return _value; }
        Systems& operator=(const uint16_t& other)
        {
            _value = other;
            return *this;
        }
    } systems = 0;
    uint8_t minCno = 0;
    uint8_t minElev = 0;
    uint8_t maxSats = 0;
    union SbasMode
    {
        uint8_t _value;
        struct
        {
            uint8_t ranging : 1;
            uint8_t diffCorr : 1;
            uint8_t integrity : 1;
            uint8_t testMode : 1;
        };
        SbasMode(const uint8_t& other) : _value(other) {}
        operator uint8_t() const { return _value; }
        SbasMode& operator=(const uint8_t& other)
        {
            _value = other;
            return *this;
        }
    } sbasMode = 0;
    union SbasSelect1
    {
        uint16_t _value;
        struct
        {
            uint16_t sbas120 : 1;
            uint16_t sbas121 : 1;
            uint16_t sbas122 : 1;
            uint16_t sbas123 : 1;
            uint16_t sbas124 : 1;
            uint16_t sbas125 : 1;
            uint16_t sbas126 : 1;
            uint16_t sbas127 : 1;
            uint16_t sbas128 : 1;
            uint16_t sbas129 : 1;
            uint16_t sbas130 : 1;
            uint16_t sbas131 : 1;
            uint16_t sbas132 : 1;
            uint16_t sbas133 : 1;
            uint16_t sbas134 : 1;
            uint16_t sbas135 : 1;
        };
        SbasSelect1(const uint16_t& other) : _value(other) {}
        operator uint16_t() const { return _value; }
        SbasSelect1& operator=(const uint16_t& other)
        {
            _value = other;
            return *this;
        }
    } sbasSelect1 = 0;
    union SbasSelect2
    {
        uint16_t _value;
        struct
        {
            uint16_t sbas136 : 1;
            uint16_t sbas137 : 1;
            uint16_t sbas138 : 1;
            uint16_t sbas139 : 1;
            uint16_t sbas140 : 1;
            uint16_t sbas141 : 1;
            uint16_t sbas142 : 1;
            uint16_t sbas143 : 1;
            uint16_t sbas144 : 1;
            uint16_t sbas145 : 1;
            uint16_t sbas146 : 1;
            uint16_t sbas147 : 1;
            uint16_t sbas148 : 1;
            uint16_t sbas149 : 1;
            uint16_t sbas150 : 1;
            uint16_t sbas151 : 1;
        };
        SbasSelect2(const uint16_t& other) : _value(other) {}
        operator uint16_t() const { return _value; }
        SbasSelect2& operator=(const uint16_t& other)
        {
            _value = other;
            return *this;
        }
    } sbasSelect2 = 0;
    union SbasSelect3
    {
        uint16_t _value;
        struct
        {
            uint16_t sbas152 : 1;
            uint16_t sbas153 : 1;
            uint16_t sbas154 : 1;
            uint16_t sbas155 : 1;
            uint16_t sbas156 : 1;
            uint16_t sbas157 : 1;
            uint16_t sbas158 : 1;
        };
        SbasSelect3(const uint16_t& other) : _value(other) {}
        operator uint16_t() const { return _value; }
        SbasSelect3& operator=(const uint16_t& other)
        {
            _value = other;
            return *this;
        }
    } sbasSelect3 = 0;
    enum class ReceiverSelect : uint8_t
    {
        GnssAB = 0,
        GnssA = 1,
        GnssB = 2,
    };
    ReceiverSelect receiverSelect = static_cast<ReceiverSelect>(0);

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
    Command toReadCommand() override;
};

bool operator==(const GnssSystemConfig& lhs, const GnssSystemConfig& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 100 - GNSS Sync Configuration </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class GnssSyncConfig : public ConfigurationRegister
{
public:
    GnssSyncConfig() : ConfigurationRegister(100) {}
    std::string name() override { return "GnssSyncConfig"; };

    enum class GnssSyncEnable : uint8_t
    {
        Off = 0,
        AlwaysOn = 1,
        OnWhenLocked = 2,
    };
    GnssSyncEnable gnssSyncEnable = static_cast<GnssSyncEnable>(0);
    enum class Polarity : uint8_t
    {
        FallingEdge = 0,
        RisingEdge = 1,
    };
    Polarity polarity = static_cast<Polarity>(0);
    enum class SpecType : uint8_t
    {
        PeriodPulseWidth = 0,
        FreqDutyCycle = 1,
    };
    SpecType specType = static_cast<SpecType>(0);
    uint8_t resv = 0;
    uint32_t period = 0;
    uint32_t pulseWidth = 0;
    int32_t offset = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const GnssSyncConfig& lhs, const GnssSyncConfig& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 103 - GNSS 2 Solution - LLA </summary>
<remarks>

Estimated GNSS 2 Solution with lat/lon/alt position. This register is deprecated and will be
removed in future firmware versions.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class Gnss2SolLla : public MeasurementRegister
{
public:
    Gnss2SolLla() : MeasurementRegister(103) {}
    std::string name() override { return "Gnss2SolLla"; };

    double gpsTow = 0;      ///< GPS time of week.
    uint16_t gps2Week = 0;  ///< The current GPS week.
    enum class Gnss2Fix : uint8_t
    {
        NoFix = 0,     ///< Not tracking.
        TimeFix = 1,   ///< Time-only fix.
        Fix2D = 2,     ///< Time plus 2D position fix.
        Fix3D = 3,     ///< Time plus full 3D position fix.
        SBAS = 4,      ///< 3D Fix using SBAS corrections.
        RtkFloat = 7,  ///< 3D Fix using an approximate RTK solution.
        RtkFix = 8,    ///< 3D Fix using a fixed RTK solution.
    };
    Gnss2Fix gnss2Fix = static_cast<Gnss2Fix>(0);
    uint8_t gnss2NumSats = 0;        ///< Number of satellites tracked by GNSS receiver.
    double lat = 0;                  ///< GNSS geodetic latitude.
    double lon = 0;                  ///< GNSS longitude.
    double alt = 0;                  ///< GNSS altitude above WGS84 ellipsoid.
    float velN = 0;                  ///< GNSS velocity in North direction.
    float velE = 0;                  ///< GNSS velocity in East direction.
    float velD = 0;                  ///< GNSS velocity in Down direction.
    float posUncertaintyN = 0;       ///< GNSS position uncertainty, North direction.
    float posUncertaintyE = 0;       ///< GNSS position uncertainty, East direction.
    float posUncertaintyD = 0;       ///< GNSS position uncertainty, Down direction.
    float gnss2VelUncertainty = 0;   ///< GNSS velocity uncertainty (scalar).
    float gnss2TimeUncertainty = 0;  ///< GNSS time uncertainty.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const Gnss2SolLla& lhs, const Gnss2SolLla& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 104 - GNSS 2 Solution - ECEF </summary>
<remarks>

Estimated GNSS 2 Solution with ECEF position. This register is deprecated and will be removed in
future firmware versions.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class Gnss2SolEcef : public MeasurementRegister
{
public:
    Gnss2SolEcef() : MeasurementRegister(104) {}
    std::string name() override { return "Gnss2SolEcef"; };

    double gpsTow = 0;      ///< GPS time of week.
    uint16_t gps2Week = 0;  ///< The current GPS week.
    enum class Gnss2Fix : uint8_t
    {
        NoFix = 0,     ///< Not tracking.
        TimeFix = 1,   ///< Time-only fix.
        Fix2D = 2,     ///< Time plus 2D position fix.
        Fix3D = 3,     ///< Time plus full 3D position fix.
        SBAS = 4,      ///< 3D Fix using SBAS corrections.
        RtkFloat = 7,  ///< 3D Fix using an approximate RTK solution.
        RtkFix = 8,    ///< 3D Fix using a fixed RTK solution.
    };
    Gnss2Fix gnss2Fix = static_cast<Gnss2Fix>(0);
    uint8_t gnss2NumSats = 0;        ///< Number of satellites tracked by GNSS receiver.
    double posX = 0;                 ///< GNSS position in ECEF-frame x-axis.
    double posY = 0;                 ///< GNSS position in ECEF-frame y-axis.
    double posZ = 0;                 ///< GNSS position in ECEF-frame z-axis.
    float velX = 0;                  ///< GNSS velocity in ECEF-frame x-axis.
    float velY = 0;                  ///< GNSS velocity in ECEF-frame y-axis.
    float velZ = 0;                  ///< GNSS velocity in ECEF-frame z-axis.
    float posUncertaintyX = 0;       ///< GNSS position uncertainty ECEF X.
    float posUncertaintyY = 0;       ///< GNSS position uncertainty ECEF Y.
    float posUncertaintyZ = 0;       ///< GNSS position uncertainty ECEF Z.
    float gnss2VelUncertainty = 0;   ///< GNSS velocity uncertainty (scalar).
    float gnss2TimeUncertainty = 0;  ///< GNSS time uncertainty.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const Gnss2SolEcef& lhs, const Gnss2SolEcef& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 157 - External GNSS Antenna Offset </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class ExtGnssOffset : public ConfigurationRegister
{
public:
    ExtGnssOffset() : ConfigurationRegister(157) {}
    std::string name() override { return "ExtGnssOffset"; };

    float positionX = 0;
    float positionY = 0;
    float positionZ = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const ExtGnssOffset& lhs, const ExtGnssOffset& rhs);
}  // namespace GNSS

namespace GNSSCompass
{
/**--------------------------------------------------------------------------------------------------
<summary> Register 86 - GNSS Compass Signal Health Status </summary>
<remarks>

Provides several indicators that serve as an overall health status of the GNSS compass subsystem.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class GnssCompassSignalHealthStatus : public MeasurementRegister
{
public:
    GnssCompassSignalHealthStatus() : MeasurementRegister(86) {}
    std::string name() override { return "GnssCompassSignalHealthStatus"; };

    float numSatsPvtA = 0;    ///< Number of satellites available for PVT solution for receiver A.
    float numSatsRtkA = 0;    ///< Number of satellites available for RTK solution for receiver A.
    float highestCn0A = 0;    ///< Highest CN0 reported on receiver A.
    float numSatsPvtB = 0;    ///< Number of satellites available for the PVT solution of receiver B.
    float numSatsRtkB = 0;    ///< Number of satellites available for the RTK solution of receiver B.
    float highestCn0B = 0;    ///< Highest CN0 reported on receiver B.
    float numComSatsPVT = 0;  ///< The number of common satellites that are used in the PVT solutions of both receiver A and receiver B.
    float numComSatsRTK = 0;  ///< The number of common satellites that are used in the RTK solutions of both receiver A and receiver B.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const GnssCompassSignalHealthStatus& lhs, const GnssCompassSignalHealthStatus& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 93 - GNSS Compass Antenna Baseline </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class GnssCompassBaseline : public ConfigurationRegister
{
public:
    GnssCompassBaseline() : ConfigurationRegister(93) {}
    std::string name() override { return "GnssCompassBaseline"; };

    float positionX = 0;
    float positionY = 0;
    float positionZ = 0;
    float uncertaintyX = 0;
    float uncertaintyY = 0;
    float uncertaintyZ = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const GnssCompassBaseline& lhs, const GnssCompassBaseline& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 97 - GNSS Compass Estimated Baseline </summary>
<remarks>

Provides the estimated GNSS compass baseline measurement. The estimated position offset and
measurement uncertainty is for antenna B relative to antenna A in the body reference frame.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class GnssCompassEstBaseline : public MeasurementRegister
{
public:
    GnssCompassEstBaseline() : MeasurementRegister(97) {}
    std::string name() override { return "GnssCompassEstBaseline"; };

    uint8_t estBaselineComplete = 0;  ///< Set to 1 when baseline estimation has complete.
    uint8_t resv = 0;                 ///< Reserved.
    uint16_t numMeas = 0;             ///< Number of measurements used by the estimated solution.
    float positionX = 0;              ///< Position of GNSS antenna B with respect to antenna A in the body-frame x-axis.
    float positionY = 0;              ///< Position of GNSS antenna B with respect to antenna A in the body-frame y-axis.
    float positionZ = 0;              ///< Position of GNSS antenna B with respect to antenna A in the body-frame z-axis.
    float uncertaintyX = 0;           ///< Uncertainty in the x-axis position measurement.
    float uncertaintyY = 0;           ///< Uncertainty in the y-axis position measurement.
    float uncertaintyZ = 0;           ///< Uncertainty in the z-axis position measurement.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const GnssCompassEstBaseline& lhs, const GnssCompassEstBaseline& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 98 - GNSS Compass Startup Status </summary>
<remarks>

Provides status information on the GNSS compass startup process.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class GnssCompassStartupStatus : public MeasurementRegister
{
public:
    GnssCompassStartupStatus() : MeasurementRegister(98) {}
    std::string name() override { return "GnssCompassStartupStatus"; };

    uint8_t percentComplete = 0;  ///< The estimated startup process completion percent.
    float currentHeading = 0;     ///< The current GNSS compass heading estimate.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const GnssCompassStartupStatus& lhs, const GnssCompassStartupStatus& rhs);
}  // namespace GNSSCompass

namespace HardSoftIronEstimator
{
/**--------------------------------------------------------------------------------------------------
<summary> Register 44 - Real-Time HSI Control </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class RealTimeHsiControl : public ConfigurationRegister
{
public:
    RealTimeHsiControl() : ConfigurationRegister(44) {}
    std::string name() override { return "RealTimeHsiControl"; };

    enum class Mode : uint8_t
    {
        Off = 0,
        Run = 1,
        Reset = 2,
    };
    Mode mode = static_cast<Mode>(0);
    enum class ApplyCompensation : uint8_t
    {
        Disable = 1,
        Enable = 3,
    };
    ApplyCompensation applyCompensation = static_cast<ApplyCompensation>(0);
    uint8_t convergeRate = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const RealTimeHsiControl& lhs, const RealTimeHsiControl& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 47 - Real-Time HSI Results </summary>
<remarks>

Magnetometer calibration values calculated by the real-time HSI calibration filter.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class EstMagCal : public MeasurementRegister
{
public:
    EstMagCal() : MeasurementRegister(47) {}
    std::string name() override { return "EstMagCal"; };

    float magGain00 = 0;  ///< Magnetometer calibration gain matrix, row 0, colum 0.
    float magGain01 = 0;  ///< Magnetometer calibration gain matrix, row 0, colum 1.
    float magGain02 = 0;  ///< Magnetometer calibration gain matrix, row 0, colum 2.
    float magGain10 = 0;  ///< Magnetometer calibration gain matrix, row 1, colum 0.
    float magGain11 = 0;  ///< Magnetometer calibration gain matrix, row 1, colum 1.
    float magGain12 = 0;  ///< Magnetometer calibration gain matrix, row 1, colum 2.
    float magGain20 = 0;  ///< Magnetometer calibration gain matrix, row 2, colum 0.
    float magGain21 = 0;  ///< Magnetometer calibration gain matrix, row 2, colum 1.
    float magGain22 = 0;  ///< Magnetometer calibration gain matrix, row 2, colum 2.
    float magBiasX = 0;   ///< Magnetometer bias calibration, sensor-frame x-axis.
    float magBiasY = 0;   ///< Magnetometer bias calibration, sensor-frame y-axis.
    float magBiasZ = 0;   ///< Magnetometer bias calibration, sensor-frame z-axis.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const EstMagCal& lhs, const EstMagCal& rhs);
}  // namespace HardSoftIronEstimator

namespace Heave
{
/**--------------------------------------------------------------------------------------------------
<summary> Register 115 - Heave and Heave Rate </summary>
<remarks>

Real-time heave and heave-rate estimates, plus a delayed-heave estimate.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class HeaveOutputs : public MeasurementRegister
{
public:
    HeaveOutputs() : MeasurementRegister(115) {}
    std::string name() override { return "HeaveOutputs"; };

    float heave = 0;         ///< Real-time heave estimate.
    float heaveRate = 0;     ///< Real-time heave rate estimate.
    float delayedHeave = 0;  ///< Delayed heave. Higher accuracy than real-time heave, but 50 seconds old.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const HeaveOutputs& lhs, const HeaveOutputs& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 116 - Heave Basic Configuration </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class HeaveBasicConfig : public ConfigurationRegister
{
public:
    HeaveBasicConfig() : ConfigurationRegister(116) {}
    std::string name() override { return "HeaveBasicConfig"; };

    float initialWavePeriod = 0;
    float initialWaveAmplitude = 0;
    float maxWavePeriod = 0;
    float minWaveAmplitude = 0;
    float delayedHeaveCutoffFreq = 0;
    float heaveCutoffFreq = 0;
    float heaveRateCutoffFreq = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const HeaveBasicConfig& lhs, const HeaveBasicConfig& rhs);
}  // namespace Heave

namespace IMU
{
/**--------------------------------------------------------------------------------------------------
<summary> Register 17 - Compensated Magnetometer </summary>
<remarks>

Compensated magnetometer measurements.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class Mag : public MeasurementRegister
{
public:
    Mag() : MeasurementRegister(17) {}
    std::string name() override { return "Mag"; };

    float magX = 0;  ///< Compensated magnetometer measurement in the body-frame x-axis.
    float magY = 0;  ///< Compensated magnetometer measurement in the body-frame y-axis.
    float magZ = 0;  ///< Compensated magnetometer measurement in the body-frame z-axis.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const Mag& lhs, const Mag& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 18 - Compensated Accelerometer </summary>
<remarks>

Compensated acceleration measurements

</remarks>
*-----------------------------------------------------------------------------------------------**/
class Accel : public MeasurementRegister
{
public:
    Accel() : MeasurementRegister(18) {}
    std::string name() override { return "Accel"; };

    float accelX = 0;  ///< Compensated accelerometer measurement in the body-frame x-axis.
    float accelY = 0;  ///< Compensated accelerometer measurement in the body-frame y-axis.
    float accelZ = 0;  ///< Compensated accelerometer measurement in the body-frame z-axis.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const Accel& lhs, const Accel& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 19 - Compensated Gyro </summary>
<remarks>

Compensated angular rate measurements.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class Gyro : public MeasurementRegister
{
public:
    Gyro() : MeasurementRegister(19) {}
    std::string name() override { return "Gyro"; };

    float gyroX = 0;  ///< Compensated angular rate measurement in the body-frame x-axis.
    float gyroY = 0;  ///< Compensated angular rate measurement in the body-frame y-axis.
    float gyroZ = 0;  ///< Compensated angular rate measurement in the body-frame z-axis.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const Gyro& lhs, const Gyro& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 20 - Compensated IMU </summary>
<remarks>

Compensated magnetic, acceleration, and angular rate measurements.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class MagAccelGyro : public MeasurementRegister
{
public:
    MagAccelGyro() : MeasurementRegister(20) {}
    std::string name() override { return "MagAccelGyro"; };

    float magX = 0;    ///< Compensated magnetometer measurement in the body-frame x-axis.
    float magY = 0;    ///< Compensated magnetometer measurement in the body-frame y-axis.
    float magZ = 0;    ///< Compensated magnetometer measurement in the body-frame z-axis.
    float accelX = 0;  ///< Compensated accelerometer measurement in the body-frame x-axis.
    float accelY = 0;  ///< Compensated accelerometer measurement in the body-frame y-axis.
    float accelZ = 0;  ///< Compensated accelerometer measurement in the body-frame z-axis.
    float gyroX = 0;   ///< Compensated angular rate measurement in the body-frame x-axis.
    float gyroY = 0;   ///< Compensated angular rate measurement in the body-frame y-axis.
    float gyroZ = 0;   ///< Compensated angular rate measurement in the body-frame z-axis.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const MagAccelGyro& lhs, const MagAccelGyro& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 23 - Magnetometer Calibration </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class MagCal : public ConfigurationRegister
{
public:
    MagCal() : ConfigurationRegister(23) {}
    std::string name() override { return "MagCal"; };

    float magGain00 = 0;
    float magGain01 = 0;
    float magGain02 = 0;
    float magGain10 = 0;
    float magGain11 = 0;
    float magGain12 = 0;
    float magGain20 = 0;
    float magGain21 = 0;
    float magGain22 = 0;
    float magBiasX = 0;
    float magBiasY = 0;
    float magBiasZ = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const MagCal& lhs, const MagCal& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 25 - Accel Calibration </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class AccelCal : public ConfigurationRegister
{
public:
    AccelCal() : ConfigurationRegister(25) {}
    std::string name() override { return "AccelCal"; };

    float accelGain00 = 0;
    float accelGain01 = 0;
    float accelGain02 = 0;
    float accelGain10 = 0;
    float accelGain11 = 0;
    float accelGain12 = 0;
    float accelGain20 = 0;
    float accelGain21 = 0;
    float accelGain22 = 0;
    float accelBiasX = 0;
    float accelBiasY = 0;
    float accelBiasZ = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const AccelCal& lhs, const AccelCal& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 26 - Reference Frame Rotation </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class RefFrameRot : public ConfigurationRegister
{
public:
    RefFrameRot() : ConfigurationRegister(26) {}
    std::string name() override { return "RefFrameRot"; };

    float rFR00 = 0;
    float rFR01 = 0;
    float rFR02 = 0;
    float rFR10 = 0;
    float rFR11 = 0;
    float rFR12 = 0;
    float rFR20 = 0;
    float rFR21 = 0;
    float rFR22 = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const RefFrameRot& lhs, const RefFrameRot& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 54 - IMU Measurements </summary>
<remarks>

Provides the calibrated IMU measurements including barometric pressure.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class ImuMeas : public MeasurementRegister
{
public:
    ImuMeas() : MeasurementRegister(54) {}
    std::string name() override { return "ImuMeas"; };

    float uncompMagX = 0;   ///< Magnetometer body-frame x-axis.
    float uncompMagY = 0;   ///< Magnetometer body-frame y-axis.
    float uncompMagZ = 0;   ///< Magnetometer body-frame z-axis.
    float uncompAccX = 0;   ///< Accelerometer body-frame x-axis.
    float uncompAccY = 0;   ///< Accelerometer body-frame y-axis.
    float uncompAccZ = 0;   ///< Accelerometer body-frame z-axis.
    float uncompGyroX = 0;  ///< Angular rate body-frame x-axis.
    float uncompGyroY = 0;  ///< Angular rate body-frame y-axis.
    float uncompGyroZ = 0;  ///< Angular rate body-frame z-axis.
    float temperature = 0;  ///< Sensor temperature.
    float pressure = 0;     ///< Barometric pressure.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const ImuMeas& lhs, const ImuMeas& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 80 - Delta Theta and Delta Velocity </summary>
<remarks>

This register contains the output values of the onboard coning and sculling algorithm.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class DeltaThetaVelocity : public MeasurementRegister
{
public:
    DeltaThetaVelocity() : MeasurementRegister(80) {}
    std::string name() override { return "DeltaThetaVelocity"; };

    float deltaTime = 0;    ///< Duration of integration interval.
    float deltaThetaX = 0;  ///< Integrated rotation vector x-axis.
    float deltaThetaY = 0;  ///< Integrated rotation vector y-axis.
    float deltaThetaZ = 0;  ///< Integrated rotation vector z-axis.
    float deltaVelX = 0;    ///< Integrated velocity x-axis.
    float deltaVelY = 0;    ///< Integrated velocity y-axis.
    float deltaVelZ = 0;    ///< Integrated velocity z-axis.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const DeltaThetaVelocity& lhs, const DeltaThetaVelocity& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 82 - Delta Theta and Delta Velocity Configuration </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class DeltaThetaVelConfig : public ConfigurationRegister
{
public:
    DeltaThetaVelConfig() : ConfigurationRegister(82) {}
    std::string name() override { return "DeltaThetaVelConfig"; };

    enum class IntegrationFrame : uint8_t
    {
        Body = 0,
        NED = 1,
    };
    IntegrationFrame integrationFrame = static_cast<IntegrationFrame>(0);
    enum class GyroCompensation : uint8_t
    {
        None = 0,
        Bias = 1,
    };
    GyroCompensation gyroCompensation = static_cast<GyroCompensation>(0);
    enum class AccelCompensation : uint8_t
    {
        None = 0,
        Gravity = 1,
        Bias = 2,
        BiasAndGravity = 3,
    };
    AccelCompensation accelCompensation = static_cast<AccelCompensation>(0);
    enum class EarthRateCompensation : uint8_t
    {
        None = 0,
        GyroRate = 1,
        CoriolisAccel = 2,
        RateAndCoriolis = 3,
    };
    EarthRateCompensation earthRateCompensation = static_cast<EarthRateCompensation>(0);
    uint16_t resv = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const DeltaThetaVelConfig& lhs, const DeltaThetaVelConfig& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 84 - Gyro Calibration </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class GyroCal : public ConfigurationRegister
{
public:
    GyroCal() : ConfigurationRegister(84) {}
    std::string name() override { return "GyroCal"; };

    float gyroGain00 = 0;
    float gyroGain01 = 0;
    float gyroGain02 = 0;
    float gyroGain10 = 0;
    float gyroGain11 = 0;
    float gyroGain12 = 0;
    float gyroGain20 = 0;
    float gyroGain21 = 0;
    float gyroGain22 = 0;
    float gyroBiasX = 0;
    float gyroBiasY = 0;
    float gyroBiasZ = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const GyroCal& lhs, const GyroCal& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 85 - IMU Filtering Configuration </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class ImuFilterControl : public ConfigurationRegister
{
public:
    ImuFilterControl() : ConfigurationRegister(85) {}
    std::string name() override { return "ImuFilterControl"; };

    uint16_t magWindowSize = 0;
    uint16_t accelWindowSize = 0;
    uint16_t gyroWindowSize = 0;
    uint16_t tempWindowSize = 0;
    uint16_t presWindowSize = 0;
    union MagFilterMode
    {
        uint8_t _value;
        struct
        {
            uint8_t uncomp : 1;
            uint8_t comp : 1;
        };
        MagFilterMode(const uint8_t& other) : _value(other) {}
        operator uint8_t() const { return _value; }
        MagFilterMode& operator=(const uint8_t& other)
        {
            _value = other;
            return *this;
        }
    } magFilterMode = 0;
    union AccelFilterMode
    {
        uint8_t _value;
        struct
        {
            uint8_t uncomp : 1;
            uint8_t comp : 1;
        };
        AccelFilterMode(const uint8_t& other) : _value(other) {}
        operator uint8_t() const { return _value; }
        AccelFilterMode& operator=(const uint8_t& other)
        {
            _value = other;
            return *this;
        }
    } accelFilterMode = 0;
    union GyroFilterMode
    {
        uint8_t _value;
        struct
        {
            uint8_t uncomp : 1;
            uint8_t comp : 1;
        };
        GyroFilterMode(const uint8_t& other) : _value(other) {}
        operator uint8_t() const { return _value; }
        GyroFilterMode& operator=(const uint8_t& other)
        {
            _value = other;
            return *this;
        }
    } gyroFilterMode = 0;
    union TempFilterMode
    {
        uint8_t _value;
        struct
        {
            uint8_t uncomp : 1;
            uint8_t comp : 1;
        };
        TempFilterMode(const uint8_t& other) : _value(other) {}
        operator uint8_t() const { return _value; }
        TempFilterMode& operator=(const uint8_t& other)
        {
            _value = other;
            return *this;
        }
    } tempFilterMode = 0;
    union PresFilterMode
    {
        uint8_t _value;
        struct
        {
            uint8_t uncomp : 1;
            uint8_t comp : 1;
        };
        PresFilterMode(const uint8_t& other) : _value(other) {}
        operator uint8_t() const { return _value; }
        PresFilterMode& operator=(const uint8_t& other)
        {
            _value = other;
            return *this;
        }
    } presFilterMode = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const ImuFilterControl& lhs, const ImuFilterControl& rhs);
}  // namespace IMU

namespace INS
{
/**--------------------------------------------------------------------------------------------------
<summary> Register 63 - INS Solution - LLA </summary>
<remarks>

Estimated INS solution with lat/lon/alt position.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class InsSolLla : public MeasurementRegister
{
public:
    InsSolLla() : MeasurementRegister(63) {}
    std::string name() override { return "InsSolLla"; };

    double gpsTow = 0;         ///< GPS time of week.
    uint16_t timeGpsWeek = 0;  ///< The current GPS week.
    InsStatus insStatus = 0;   ///< Ins status bitfield.
    float yaw = 0;             ///< Yaw angle.
    float pitch = 0;           ///< Pitch angle.
    float roll = 0;            ///< Roll angle.
    double posLat = 0;         ///< Geodetic latitude.
    double posLon = 0;         ///< Longitude.
    double posAlt = 0;         ///< Altitude above WGS84 ellipsoid.
    float velN = 0;            ///< Velocity in North direction.
    float velE = 0;            ///< Velocity in East direction.
    float velD = 0;            ///< Velocity in Down direction.
    float attUncertainty = 0;  ///< Filter estimated attitude uncertainty.
    float posUncertainty = 0;  ///< Filter estimated position uncertainty.
    float velUncertainty = 0;  ///< Filter estimated velocity uncertainty.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const InsSolLla& lhs, const InsSolLla& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 64 - INS Solution - ECEF </summary>
<remarks>

Estimated INS Solution with ECEF position

</remarks>
*-----------------------------------------------------------------------------------------------**/
class InsSolEcef : public MeasurementRegister
{
public:
    InsSolEcef() : MeasurementRegister(64) {}
    std::string name() override { return "InsSolEcef"; };

    double gpsTow = 0;         ///< GPS time of week.
    uint16_t timeGpsWeek = 0;  ///< The current GPS week.
    InsStatus insStatus = 0;   ///< Ins status bitfield.
    float yaw = 0;             ///< Yaw angle.
    float pitch = 0;           ///< Pitch angle.
    float roll = 0;            ///< Roll angle.
    double posEX = 0;          ///< Position in ECEF-frame x-axis.
    double posEY = 0;          ///< Position in ECEF-frame y-axis.
    double posEZ = 0;          ///< Position in ECEF-frame z-axis.
    float velEX = 0;           ///< Velocity in ECEF-frame x-axis.
    float velEY = 0;           ///< Velocity in ECEF-frame y-axis.
    float velEZ = 0;           ///< Velocity in ECEF-frame z-axis.
    float attUncertainty = 0;  ///< Filter estimated attitude uncertainty.
    float posUncertainty = 0;  ///< Filter estimated position uncertainty.
    float velUncertainty = 0;  ///< Filter estimated velocity uncertainty.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const InsSolEcef& lhs, const InsSolEcef& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 67 - INS Basic Configuration </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class InsBasicConfig : public ConfigurationRegister
{
public:
    InsBasicConfig() : ConfigurationRegister(67) {}
    std::string name() override { return "InsBasicConfig"; };

    enum class Scenario : uint8_t
    {
        Ahrs = 0,
        GnssInsWithPressure = 1,
        GnssInsNoPressure = 2,
        DualGnssNoPressure = 3,
        DualGnssWithPressure = 5,
    };
    Scenario scenario = static_cast<Scenario>(0);
    enum class AhrsAiding : uint8_t
    {
        Disable = 0,
        Enable = 1,
    };
    AhrsAiding ahrsAiding = static_cast<AhrsAiding>(0);
    enum class EstBaseline : uint8_t
    {
        Disable = 0,
        Enable = 1,
    };
    EstBaseline estBaseline = static_cast<EstBaseline>(0);
    uint8_t resv = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const InsBasicConfig& lhs, const InsBasicConfig& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 72 - INS State - LLA </summary>
<remarks>

Estimated INS state with lat/lon/alt position.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class InsStateLla : public MeasurementRegister
{
public:
    InsStateLla() : MeasurementRegister(72) {}
    std::string name() override { return "InsStateLla"; };

    float yaw = 0;      ///< Yaw angle.
    float pitch = 0;    ///< Pitch angle.
    float roll = 0;     ///< Roll angle.
    double posLat = 0;  ///< Geodetic latitude.
    double posLon = 0;  ///< Longitude.
    double posAlt = 0;  ///< Altitude above WGS84 ellipsoid.
    float velN = 0;     ///< Velocity in North direction.
    float velE = 0;     ///< Velocity in East direction.
    float velD = 0;     ///< Velocity in Down direction.
    float accelX = 0;   ///< Compensated accelerometer measurement in the body-frame x-axis.
    float accelY = 0;   ///< Compensated accelerometer measurement in the body-frame y-axis.
    float accelZ = 0;   ///< Compensated accelerometer measurement in the body-frame z-axis.
    float gyroX = 0;    ///< Compensated angular rate measurement in the body-frame x-axis.
    float gyroY = 0;    ///< Compensated angular rate measurement in the body-frame y-axis.
    float gyroZ = 0;    ///< Compensated angular rate measurement in the body-frame z-axis.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const InsStateLla& lhs, const InsStateLla& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 73 - INS State - ECEF </summary>
<remarks>

Estimated INS state with ECEF position.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class InsStateEcef : public MeasurementRegister
{
public:
    InsStateEcef() : MeasurementRegister(73) {}
    std::string name() override { return "InsStateEcef"; };

    float yaw = 0;     ///< Yaw angle.
    float pitch = 0;   ///< Pitch angle.
    float roll = 0;    ///< Roll angle.
    double posEX = 0;  ///< Position in ECEF-frame x-axis.
    double posEY = 0;  ///< Position in ECEF-frame y-axis.
    double posEZ = 0;  ///< Position in ECEF-frame z-axis.
    float velEX = 0;   ///< Velocity in ECEF-frame x-axis.
    float velEY = 0;   ///< Velocity in ECEF-frame y-axis.
    float velEZ = 0;   ///< Velocity in ECEF-frame z-axis.
    float accelX = 0;  ///< Compensated accelerometer measurement in the body-frame x-axis.
    float accelY = 0;  ///< Compensated accelerometer measurement in the body-frame y-axis.
    float accelZ = 0;  ///< Compensated accelerometer measurement in the body-frame z-axis.
    float gyroX = 0;   ///< Compensated angular rate measurement in the body-frame x-axis.
    float gyroY = 0;   ///< Compensated angular rate measurement in the body-frame y-axis.
    float gyroZ = 0;   ///< Compensated angular rate measurement in the body-frame z-axis.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const InsStateEcef& lhs, const InsStateEcef& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 74 - Filter Startup Bias </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class FilterStartupBias : public ConfigurationRegister
{
public:
    FilterStartupBias() : ConfigurationRegister(74) {}
    std::string name() override { return "FilterStartupBias"; };

    float gyroBiasX = 0;
    float gyroBiasY = 0;
    float gyroBiasZ = 0;
    float accelBiasX = 0;
    float accelBiasY = 0;
    float accelBiasZ = 0;
    float presBias = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const FilterStartupBias& lhs, const FilterStartupBias& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 105 - INS Reference Point Offset </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class InsRefOffset : public ConfigurationRegister
{
public:
    InsRefOffset() : ConfigurationRegister(105) {}
    std::string name() override { return "InsRefOffset"; };

    float refOffsetX = 0;
    float refOffsetY = 0;
    float refOffsetZ = 0;
    float refUncertX = 0;
    float refUncertY = 0;
    float refUncertZ = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const InsRefOffset& lhs, const InsRefOffset& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 144 - INS GNSS Select </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class InsGnssSelect : public ConfigurationRegister
{
public:
    InsGnssSelect() : ConfigurationRegister(144) {}
    std::string name() override { return "InsGnssSelect"; };

    enum class ActiveReceiverSelect : uint8_t
    {
        None = 0,
        PrimaryReceiver = 1,
        SecondaryReceiver = 2,
        TertiaryReceiver = 3,
        FallbackOnFailure = 5,
    };
    ActiveReceiverSelect activeReceiverSelect = static_cast<ActiveReceiverSelect>(0);
    uint8_t usedForNavTime = 0;
    uint8_t hysteresisTime = 0;
    enum class UseGnssCompass : uint8_t
    {
        Off = 0,
        On = 1,
    };
    UseGnssCompass useGnssCompass = static_cast<UseGnssCompass>(0);
    uint8_t resv1 = 0;
    uint8_t resv2 = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const InsGnssSelect& lhs, const InsGnssSelect& rhs);
}  // namespace INS

namespace System
{
/**--------------------------------------------------------------------------------------------------
<summary> Register 0 - User Tag </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class UserTag : public ConfigurationRegister
{
public:
    UserTag() : ConfigurationRegister(0) {}
    std::string name() override { return "UserTag"; };

    AsciiMessage tag = "";

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const UserTag& lhs, const UserTag& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 1 - Model </summary>
<remarks>

Product model string.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class Model : public MeasurementRegister
{
public:
    Model() : MeasurementRegister(1) {}
    std::string name() override { return "Model"; };

    AsciiMessage model = "";  ///< Product model number, maximum length 24 characters.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const Model& lhs, const Model& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 2 - Hardware Version </summary>
<remarks>

Hardware version number.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class HwVer : public MeasurementRegister
{
public:
    HwVer() : MeasurementRegister(2) {}
    std::string name() override { return "HwVer"; };

    uint32_t hwVer = 0;     ///< Hardware version number.
    uint32_t hwMinVer = 0;  ///< Hardware minor version number.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const HwVer& lhs, const HwVer& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 3 - Serial Number </summary>
<remarks>

Device serial number.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class Serial : public MeasurementRegister
{
public:
    Serial() : MeasurementRegister(3) {}
    std::string name() override { return "Serial"; };

    uint32_t serialNum = 0;  ///< The unit's serial number.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const Serial& lhs, const Serial& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 4 - Firmware Version </summary>
<remarks>

Firmware version number.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class FwVer : public MeasurementRegister
{
public:
    FwVer() : MeasurementRegister(4) {}
    std::string name() override { return "FwVer"; };

    AsciiMessage fwVer = "";  ///< Firmware version.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const FwVer& lhs, const FwVer& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 5 - Baud Rate </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class BaudRate : public ConfigurationRegister
{
public:
    BaudRate() : ConfigurationRegister(5) {}
    std::string name() override { return "BaudRate"; };

    enum class BaudRates : uint32_t
    {
        Baud9600 = 9600,
        Baud19200 = 19200,
        Baud38400 = 38400,
        Baud57600 = 57600,
        Baud115200 = 115200,
        Baud128000 = 128000,
        Baud230400 = 230400,
        Baud460800 = 460800,
        Baud921600 = 921600,
    };
    BaudRates baudRate = static_cast<BaudRates>(0);
    enum class SerialPort : uint8_t
    {
        ActiveSerial = 0,  ///< Set value on active serial port.
        Serial1 = 1,
        Serial2 = 2,
        Poll = static_cast<uint8_t>('?')
    };
    SerialPort serialPort = static_cast<SerialPort>(0);

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
    Command toReadCommand() override;
};

bool operator==(const BaudRate& lhs, const BaudRate& rhs);

inline std::ostream& operator<<(std::ostream& outStream, const BaudRate::BaudRates& baudrate) noexcept { return outStream << static_cast<uint32_t>(baudrate); }
/**--------------------------------------------------------------------------------------------------
<summary> Register 6 - Async Data Output Type </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class AsyncOutputType : public ConfigurationRegister
{
public:
    AsyncOutputType() : ConfigurationRegister(6) {}
    std::string name() override { return "AsyncOutputType"; };

    enum class Ador : uint32_t
    {
        OFF = 0,
        YPR = 1,
        QTN = 2,
        QMR = 8,
        MAG = 10,
        ACC = 11,
        GYR = 12,
        MAR = 13,
        YMR = 14,
        YBA = 16,
        YIA = 17,
        IMU = 19,
        GPS = 20,
        GPE = 21,
        INS = 22,
        INE = 23,
        ISL = 28,
        ISE = 29,
        DTV = 30,
        G2S = 32,
        G2E = 33,
        HVE = 34,
    };
    Ador ador = static_cast<Ador>(0);
    enum class SerialPort : uint8_t
    {
        ActiveSerial = 0,  ///< Set value on active serial port.
        Serial1 = 1,
        Serial2 = 2,
        Poll = static_cast<uint8_t>('?')
    };
    SerialPort serialPort = static_cast<SerialPort>(0);

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
    Command toReadCommand() override;
};

bool operator==(const AsyncOutputType& lhs, const AsyncOutputType& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 7 - Async Data Output Freq </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class AsyncOutputFreq : public ConfigurationRegister
{
public:
    AsyncOutputFreq() : ConfigurationRegister(7) {}
    std::string name() override { return "AsyncOutputFreq"; };

    enum class Adof : uint32_t
    {
        Rate0Hz = 0,
        Rate1Hz = 1,
        Rate2Hz = 2,
        Rate4Hz = 4,
        Rate5Hz = 5,
        Rate10Hz = 10,
        Rate20Hz = 20,
        Rate25Hz = 25,
        Rate40Hz = 40,
        Rate50Hz = 50,
        Rate100Hz = 100,
        Rate200Hz = 200,
    };
    Adof adof = static_cast<Adof>(0);
    enum class SerialPort : uint8_t
    {
        ActiveSerial = 0,  ///< Set value on active serial port.
        Serial1 = 1,
        Serial2 = 2,
        Poll = static_cast<uint8_t>('?')
    };
    SerialPort serialPort = static_cast<SerialPort>(0);

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
    Command toReadCommand() override;
};

bool operator==(const AsyncOutputFreq& lhs, const AsyncOutputFreq& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 30 - Communication Protocol Control </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class ProtocolControl : public ConfigurationRegister
{
public:
    ProtocolControl() : ConfigurationRegister(30) {}
    std::string name() override { return "ProtocolControl"; };

    enum class AsciiAppendCount : uint8_t
    {
        None = 0,
        SyncInCount = 1,
        SyncInTime = 2,
        SyncOutCount = 3,
        GpsPps = 4,
        GpsTow = 5,
    };
    AsciiAppendCount asciiAppendCount = static_cast<AsciiAppendCount>(0);
    enum class AsciiAppendStatus : uint8_t
    {
        None = 0,
        Ahrs = 1,
        Ins = 2,
        Imu = 3,
        Gnss1 = 4,
        Gnss2 = 5,
        Gnss3 = 6,
    };
    AsciiAppendStatus asciiAppendStatus = static_cast<AsciiAppendStatus>(0);
    enum class SpiAppendCount : uint8_t
    {
        None = 0,
        SyncInCount = 1,
        SyncInTime = 2,
        SyncOutCount = 3,
        GpsPps = 4,
        GpsTow = 5,
    };
    SpiAppendCount spiAppendCount = static_cast<SpiAppendCount>(0);
    enum class SpiAppendStatus : uint8_t
    {
        None = 0,
        Ahrs = 1,
        Ins = 2,
        Imu = 3,
        Gnss1 = 4,
        Gnss2 = 5,
        Gnss3 = 6,
    };
    SpiAppendStatus spiAppendStatus = static_cast<SpiAppendStatus>(0);
    enum class AsciiChecksum : uint8_t
    {
        Checksum8bit = 1,
        Crc16bit = 3,
    };
    AsciiChecksum asciiChecksum = static_cast<AsciiChecksum>(0);
    enum class SpiChecksum : uint8_t
    {
        Off = 0,
        Checksum8bit = 1,
        Crc16bit = 3,
    };
    SpiChecksum spiChecksum = static_cast<SpiChecksum>(0);
    enum class ErrorMode : uint8_t
    {
        Ignore = 0,
        SendError = 1,
        AdorOff = 2,
    };
    ErrorMode errorMode = static_cast<ErrorMode>(0);

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const ProtocolControl& lhs, const ProtocolControl& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 32 - Synchronization Control </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class SyncControl : public ConfigurationRegister
{
public:
    SyncControl() : ConfigurationRegister(32) {}
    std::string name() override { return "SyncControl"; };

    enum class SyncInMode : uint8_t
    {
        Disable = 0,
        Count = 3,
        ImuSample = 4,
        AsyncAll = 5,
        Async0 = 6,
    };
    SyncInMode syncInMode = static_cast<SyncInMode>(0);
    enum class SyncInEdge : uint8_t
    {
        RisingEdge = 0,
        FallingEdge = 1,
    };
    SyncInEdge syncInEdge = static_cast<SyncInEdge>(0);
    uint16_t syncInSkipFactor = 0;
    uint32_t resv1 = 0;
    enum class SyncOutMode : uint8_t
    {
        None = 0,
        ImuStart = 1,
        ImuReady = 2,
        NavReady = 3,
        GpsPps = 6,
    };
    SyncOutMode syncOutMode = static_cast<SyncOutMode>(0);
    enum class SyncOutPolarity : uint8_t
    {
        NegativePulse = 0,
        PositivePulse = 1,
    };
    SyncOutPolarity syncOutPolarity = static_cast<SyncOutPolarity>(0);
    uint16_t syncOutSkipFactor = 0;
    uint32_t syncOutPulseWidth = 0;
    uint32_t resv2 = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const SyncControl& lhs, const SyncControl& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 33 - Synchronization Status </summary>
<remarks>

Contains counters based on the SyncIn and SyncOut events.

</remarks>
*-----------------------------------------------------------------------------------------------**/
class SyncStatus : public MeasurementRegister
{
public:
    SyncStatus() : MeasurementRegister(33) {}
    std::string name() override { return "SyncStatus"; };

    uint32_t syncInCount = 0;   ///< Amount of SyncIn Events that have occurred.
    uint32_t syncInTime = 0;    ///< The amount of time that has elapsed since the last SyncIn Event.
    uint32_t syncOutCount = 0;  ///< Keeps track of the number of times that a SyncOut pulse has occurred.

    bool fromString(const AsciiMessage& sensorResponse) override;
    using MeasurementRegister::toString;
};

bool operator==(const SyncStatus& lhs, const SyncStatus& rhs);
struct BinaryOutputMeasurements
{
    union Common
    {
        uint32_t _value;
        struct
        {
            uint32_t timeStartup : 1;
            uint32_t timeGps : 1;
            uint32_t timeSyncIn : 1;
            uint32_t ypr : 1;
            uint32_t quaternion : 1;
            uint32_t angularRate : 1;
            uint32_t posLla : 1;
            uint32_t velNed : 1;
            uint32_t accel : 1;
            uint32_t imu : 1;
            uint32_t magPres : 1;
            uint32_t deltas : 1;
            uint32_t insStatus : 1;
            uint32_t syncInCnt : 1;
            uint32_t timeGpsPps : 1;
        };
        Common() {}
        Common(const uint32_t& other) : _value(other) {}
        operator uint32_t() const
        {
            {
                return _value;
            }
        }
        Common& operator=(const uint32_t& other)
        {
            _value = other;
            return *this;
        }
    } common = 0;
    union Time
    {
        uint32_t _value;
        struct
        {
            uint32_t timeStartup : 1;
            uint32_t timeGps : 1;
            uint32_t timeGpsTow : 1;
            uint32_t timeGpsWeek : 1;
            uint32_t timeSyncIn : 1;
            uint32_t timeGpsPps : 1;
            uint32_t timeUtc : 1;
            uint32_t syncInCnt : 1;
            uint32_t syncOutCnt : 1;
            uint32_t timeStatus : 1;
        };
        Time() {}
        Time(const uint32_t& other) : _value(other) {}
        operator uint32_t() const
        {
            {
                return _value;
            }
        }
        Time& operator=(const uint32_t& other)
        {
            _value = other;
            return *this;
        }
    } time = 0;
    union Imu
    {
        uint32_t _value;
        struct
        {
            uint32_t imuStatus : 1;
            uint32_t uncompMag : 1;
            uint32_t uncompAccel : 1;
            uint32_t uncompGyro : 1;
            uint32_t temperature : 1;
            uint32_t pressure : 1;
            uint32_t deltaTheta : 1;
            uint32_t deltaVel : 1;
            uint32_t mag : 1;
            uint32_t accel : 1;
            uint32_t angularRate : 1;
            uint32_t sensSat : 1;
        };
        Imu() {}
        Imu(const uint32_t& other) : _value(other) {}
        operator uint32_t() const
        {
            {
                return _value;
            }
        }
        Imu& operator=(const uint32_t& other)
        {
            _value = other;
            return *this;
        }
    } imu = 0;
    union Gnss
    {
        uint32_t _value;
        struct
        {
            uint32_t gnss1TimeUtc : 1;
            uint32_t gps1Tow : 1;
            uint32_t gps1Week : 1;
            uint32_t gnss1NumSats : 1;
            uint32_t gnss1Fix : 1;
            uint32_t gnss1PosLla : 1;
            uint32_t gnss1PosEcef : 1;
            uint32_t gnss1VelNed : 1;
            uint32_t gnss1VelEcef : 1;
            uint32_t gnss1PosUncertainty : 1;
            uint32_t gnss1VelUncertainty : 1;
            uint32_t gnss1TimeUncertainty : 1;
            uint32_t gnss1TimeInfo : 1;
            uint32_t gnss1Dop : 1;
            uint32_t gnss1SatInfo : 1;
            uint32_t : 1;  ///< 15: Extension Bit
            uint32_t gnss1RawMeas : 1;
            uint32_t gnss1Status : 1;
            uint32_t gnss1AltMSL : 1;
        };
        Gnss() {}
        Gnss(const uint32_t& other) : _value(other) {}
        operator uint32_t() const
        {
            {
                return _value;
            }
        }
        Gnss& operator=(const uint32_t& other)
        {
            _value = other;
            return *this;
        }
    } gnss = 0;
    union Attitude
    {
        uint32_t _value;
        struct
        {
            uint32_t : 1;  ///< 0: padding
            uint32_t ypr : 1;
            uint32_t quaternion : 1;
            uint32_t dcm : 1;
            uint32_t magNed : 1;
            uint32_t accelNed : 1;
            uint32_t linBodyAcc : 1;
            uint32_t linAccelNed : 1;
            uint32_t yprU : 1;
            uint32_t : 1;  ///< 9: padding
            uint32_t : 1;  ///< 10: padding
            uint32_t : 1;  ///< 11: padding
            uint32_t heave : 1;
            uint32_t attU : 1;
        };
        Attitude() {}
        Attitude(const uint32_t& other) : _value(other) {}
        operator uint32_t() const
        {
            {
                return _value;
            }
        }
        Attitude& operator=(const uint32_t& other)
        {
            _value = other;
            return *this;
        }
    } attitude = 0;
    union Ins
    {
        uint32_t _value;
        struct
        {
            uint32_t insStatus : 1;
            uint32_t posLla : 1;
            uint32_t posEcef : 1;
            uint32_t velBody : 1;
            uint32_t velNed : 1;
            uint32_t velEcef : 1;
            uint32_t magEcef : 1;
            uint32_t accelEcef : 1;
            uint32_t linAccelEcef : 1;
            uint32_t posU : 1;
            uint32_t velU : 1;
        };
        Ins() {}
        Ins(const uint32_t& other) : _value(other) {}
        operator uint32_t() const
        {
            {
                return _value;
            }
        }
        Ins& operator=(const uint32_t& other)
        {
            _value = other;
            return *this;
        }
    } ins = 0;
    union Gnss2
    {
        uint32_t _value;
        struct
        {
            uint32_t gnss2TimeUtc : 1;
            uint32_t gps2Tow : 1;
            uint32_t gps2Week : 1;
            uint32_t gnss2NumSats : 1;
            uint32_t gnss2Fix : 1;
            uint32_t gnss2PosLla : 1;
            uint32_t gnss2PosEcef : 1;
            uint32_t gnss2VelNed : 1;
            uint32_t gnss2VelEcef : 1;
            uint32_t gnss2PosUncertainty : 1;
            uint32_t gnss2VelUncertainty : 1;
            uint32_t gnss2TimeUncertainty : 1;
            uint32_t gnss2TimeInfo : 1;
            uint32_t gnss2Dop : 1;
            uint32_t gnss2SatInfo : 1;
            uint32_t : 1;  ///< 15: Extension Bit
            uint32_t gnss2RawMeas : 1;
            uint32_t gnss2Status : 1;
            uint32_t gnss2AltMSL : 1;
        };
        Gnss2() {}
        Gnss2(const uint32_t& other) : _value(other) {}
        operator uint32_t() const
        {
            {
                return _value;
            }
        }
        Gnss2& operator=(const uint32_t& other)
        {
            _value = other;
            return *this;
        }
    } gnss2 = 0;
    union Gnss3
    {
        uint32_t _value;
        struct
        {
            uint32_t gnss3TimeUtc : 1;
            uint32_t gps3Tow : 1;
            uint32_t gps3Week : 1;
            uint32_t gnss3NumSats : 1;
            uint32_t gnss3Fix : 1;
            uint32_t gnss3PosLla : 1;
            uint32_t gnss3PosEcef : 1;
            uint32_t gnss3VelNed : 1;
            uint32_t gnss3VelEcef : 1;
            uint32_t gnss3PosUncertainty : 1;
            uint32_t gnss3VelUncertainty : 1;
            uint32_t gnss3TimeUncertainty : 1;
            uint32_t gnss3TimeInfo : 1;
            uint32_t gnss3Dop : 1;
            uint32_t gnss3SatInfo : 1;
            uint32_t : 1;  ///< 15: Extension Bit
            uint32_t gnss3RawMeas : 1;
            uint32_t gnss3Status : 1;
            uint32_t gnss3AltMSL : 1;
        };
        Gnss3() {}
        Gnss3(const uint32_t& other) : _value(other) {}
        operator uint32_t() const
        {
            {
                return _value;
            }
        }
        Gnss3& operator=(const uint32_t& other)
        {
            _value = other;
            return *this;
        }
    } gnss3 = 0;

    BinaryHeader toBinaryHeader() const noexcept;
};
class BinaryOutput : public ConfigurationRegister, public BinaryOutputMeasurements
{
public:
    using ConfigurationRegister::ConfigurationRegister;

    union AsyncMode
    {
        uint16_t _value;
        struct
        {
            uint16_t serial1 : 1;
            uint16_t serial2 : 1;
            uint16_t spi : 1;
        };
        AsyncMode() {}
        AsyncMode(const uint16_t& other) : _value(other) {}
        operator uint16_t() const { return _value; }
        AsyncMode& operator=(const uint16_t& other)
        {
            _value = other;
            return *this;
        }
    } asyncMode = 0;
    uint16_t rateDivisor = 0;

    bool fromString(const AsciiMessage& sensorResponsePayload) override;
    AsciiMessage toString() const override;

protected:
    using Register::_id;
};
inline bool operator==(const BinaryOutput& lhs, const BinaryOutput& rhs)
{
    return ((lhs.asyncMode._value == rhs.asyncMode._value) && (lhs.rateDivisor == rhs.rateDivisor) && (lhs.common._value == rhs.common._value) &&
            (lhs.time._value == rhs.time._value) && (lhs.imu._value == rhs.imu._value) && (lhs.gnss._value == rhs.gnss._value) &&
            (lhs.attitude._value == rhs.attitude._value) && (lhs.ins._value == rhs.ins._value) && (lhs.gnss2._value == rhs.gnss2._value) &&
            (lhs.gnss3._value == rhs.gnss3._value));
}

class BinaryOutput1 : public BinaryOutput
{
public:
    BinaryOutput1() : BinaryOutput{75} {}
    std::string name() override { return "BinaryOutput1"; };
};

class BinaryOutput2 : public BinaryOutput
{
public:
    BinaryOutput2() : BinaryOutput{76} {}
    std::string name() override { return "BinaryOutput2"; };
};

class BinaryOutput3 : public BinaryOutput
{
public:
    BinaryOutput3() : BinaryOutput{77} {}
    std::string name() override { return "BinaryOutput3"; };
};

/**--------------------------------------------------------------------------------------------------
<summary> Register 101 - NMEA Output 1 </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class NmeaOutput1 : public ConfigurationRegister
{
public:
    NmeaOutput1() : ConfigurationRegister(101) {}
    std::string name() override { return "NmeaOutput1"; };

    enum class Port : uint8_t
    {
        None = 0,
        Serial1 = 1,
        Serial2 = 2,
        Serial1And2 = 3,
    };
    Port port = static_cast<Port>(0);
    enum class Rate : uint8_t
    {
        Rate0Hz = 0,
        Rate1Hz = 1,
        Rate5Hz = 5,
        Rate10Hz = 10,
        Rate20Hz = 20,
    };
    Rate rate = static_cast<Rate>(0);
    enum class Mode : uint8_t
    {
        V41GPID = 0,
        V23GPID = 1,
        V41SYSID = 2,
    };
    Mode mode = static_cast<Mode>(0);
    enum class GnssSelect : uint8_t
    {
        GnssA = 0,
        ActiveReceiver = 1,
    };
    GnssSelect gnssSelect = static_cast<GnssSelect>(0);
    union MsgSelection
    {
        uint32_t _value;
        struct
        {
            uint32_t rmcGnss : 1;
            uint32_t rmcIns : 1;
            uint32_t ggaGnss : 1;
            uint32_t ggaIns : 1;
            uint32_t gllGnss : 1;
            uint32_t gllIns : 1;
            uint32_t gsaGnss : 1;
            uint32_t gsvGnss : 1;
            uint32_t hdgIns : 1;
            uint32_t hdtIns : 1;
            uint32_t thsIns : 1;
            uint32_t vtgGnss : 1;
            uint32_t vtgIns : 1;
            uint32_t zdaGnss : 1;
            uint32_t zdaIns : 1;
            uint32_t pashrIns : 1;
            uint32_t tss1Ins : 1;
            uint32_t indyn : 1;
        };
        MsgSelection(const uint32_t& other) : _value(other) {}
        operator uint32_t() const { return _value; }
        MsgSelection& operator=(const uint32_t& other)
        {
            _value = other;
            return *this;
        }
    } msgSelection = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const NmeaOutput1& lhs, const NmeaOutput1& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 102 - NMEA Output 2 </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class NmeaOutput2 : public ConfigurationRegister
{
public:
    NmeaOutput2() : ConfigurationRegister(102) {}
    std::string name() override { return "NmeaOutput2"; };

    enum class Port : uint8_t
    {
        None = 0,
        Serial1 = 1,
        Serial2 = 2,
        Serial1And2 = 3,
    };
    Port port = static_cast<Port>(0);
    enum class Rate : uint8_t
    {
        Rate0Hz = 0,
        Rate1Hz = 1,
        Rate5Hz = 5,
        Rate10Hz = 10,
        Rate20Hz = 20,
    };
    Rate rate = static_cast<Rate>(0);
    enum class Mode : uint8_t
    {
        V41GPID = 0,
        V23GPID = 1,
        V41SYSID = 2,
    };
    Mode mode = static_cast<Mode>(0);
    enum class GnssSelect : uint8_t
    {
        GnssA = 0,
        ActiveReceiver = 1,
    };
    GnssSelect gnssSelect = static_cast<GnssSelect>(0);
    union MsgSelection
    {
        uint32_t _value;
        struct
        {
            uint32_t rmcGnss : 1;
            uint32_t rmcIns : 1;
            uint32_t ggaGnss : 1;
            uint32_t ggaIns : 1;
            uint32_t gllGnss : 1;
            uint32_t gllIns : 1;
            uint32_t gsaGnss : 1;
            uint32_t gsvGnss : 1;
            uint32_t hdgIns : 1;
            uint32_t hdtIns : 1;
            uint32_t thsIns : 1;
            uint32_t vtgGnss : 1;
            uint32_t vtgIns : 1;
            uint32_t zdaGnss : 1;
            uint32_t zdaIns : 1;
            uint32_t pashrIns : 1;
            uint32_t tss1Ins : 1;
            uint32_t indyn : 1;
        };
        MsgSelection(const uint32_t& other) : _value(other) {}
        operator uint32_t() const { return _value; }
        MsgSelection& operator=(const uint32_t& other)
        {
            _value = other;
            return *this;
        }
    } msgSelection = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const NmeaOutput2& lhs, const NmeaOutput2& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 206 - Legacy Compatibility Settings </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class LegacyCompatibilitySettings : public ConfigurationRegister
{
public:
    LegacyCompatibilitySettings() : ConfigurationRegister(206) {}
    std::string name() override { return "LegacyCompatibilitySettings"; };

    enum class InsLegacy : uint8_t
    {
        Current = 0,
        Legacy = 1,
    };
    InsLegacy insLegacy = static_cast<InsLegacy>(0);
    union GnssLegacy
    {
        uint8_t _value;
        struct
        {
            uint8_t legacyGnssFix : 1;
            uint8_t requireReg55Reset : 1;
            uint8_t alwaysPpsPulse : 1;
        };
        GnssLegacy(const uint8_t& other) : _value(other) {}
        operator uint8_t() const { return _value; }
        GnssLegacy& operator=(const uint8_t& other)
        {
            _value = other;
            return *this;
        }
    } gnssLegacy = 0;
    enum class ImuLegacy : uint8_t
    {
        Current = 0,
        Legacy = 1,
    };
    ImuLegacy imuLegacy = static_cast<ImuLegacy>(0);
    enum class HwLegacy : uint8_t
    {
        Current = 0,
        Legacy = 1,
    };
    HwLegacy hwLegacy = static_cast<HwLegacy>(0);

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const LegacyCompatibilitySettings& lhs, const LegacyCompatibilitySettings& rhs);
}  // namespace System

namespace VelocityAiding
{
/**--------------------------------------------------------------------------------------------------
<summary> Register 50 - Velocity Aiding Measurement </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class VelAidingMeas : public ConfigurationRegister
{
public:
    VelAidingMeas() : ConfigurationRegister(50) {}
    std::string name() override { return "VelAidingMeas"; };

    float velocityX = 0;
    float velocityY = 0;
    float velocityZ = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const VelAidingMeas& lhs, const VelAidingMeas& rhs);
/**--------------------------------------------------------------------------------------------------
<summary> Register 51 - Velocity Aiding Control </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class VelAidingControl : public ConfigurationRegister
{
public:
    VelAidingControl() : ConfigurationRegister(51) {}
    std::string name() override { return "VelAidingControl"; };

    enum class VelAidEnable : uint8_t
    {
        Disable = 0,
        Enable = 1,
    };
    VelAidEnable velAidEnable = static_cast<VelAidEnable>(0);
    float velUncertTuning = 0;
    float resv = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const VelAidingControl& lhs, const VelAidingControl& rhs);
}  // namespace VelocityAiding

namespace WorldMagGravityModel
{
/**--------------------------------------------------------------------------------------------------
<summary> Register 83 - Reference Model Configuration </summary>
<remarks>
</remarks>
*-----------------------------------------------------------------------------------------------**/
class RefModelConfig : public ConfigurationRegister
{
public:
    RefModelConfig() : ConfigurationRegister(83) {}
    std::string name() override { return "RefModelConfig"; };

    enum class EnableMagModel : uint8_t
    {
        Disabled = 0,
        Enabled = 1,
    };
    EnableMagModel enableMagModel = static_cast<EnableMagModel>(0);
    enum class EnableGravityModel : uint8_t
    {
        Disabled = 0,
        Enabled = 1,
    };
    EnableGravityModel enableGravityModel = static_cast<EnableGravityModel>(0);
    uint8_t resv1 = 0;
    uint8_t resv2 = 0;
    uint32_t recalcThreshold = 0;
    float year = 0;
    double latitude = 0;
    double longitude = 0;
    double altitude = 0;

    bool fromString(const AsciiMessage& sensorResponse) override;
    AsciiMessage toString() const override;
};

bool operator==(const RefModelConfig& lhs, const RefModelConfig& rhs);
}  // namespace WorldMagGravityModel

namespace ById
{
using Reg0 = Registers::System::UserTag;
using Reg1 = Registers::System::Model;
using Reg2 = Registers::System::HwVer;
using Reg3 = Registers::System::Serial;
using Reg4 = Registers::System::FwVer;
using Reg5 = Registers::System::BaudRate;
using Reg6 = Registers::System::AsyncOutputType;
using Reg7 = Registers::System::AsyncOutputFreq;
using Reg8 = Registers::Attitude::YawPitchRoll;
using Reg9 = Registers::Attitude::Quaternion;
using Reg15 = Registers::Attitude::QuatMagAccelRate;
using Reg17 = Registers::IMU::Mag;
using Reg18 = Registers::IMU::Accel;
using Reg19 = Registers::IMU::Gyro;
using Reg20 = Registers::IMU::MagAccelGyro;
using Reg21 = Registers::Attitude::MagGravRefVec;
using Reg23 = Registers::IMU::MagCal;
using Reg25 = Registers::IMU::AccelCal;
using Reg26 = Registers::IMU::RefFrameRot;
using Reg27 = Registers::Attitude::YprMagAccelAngularRates;
using Reg30 = Registers::System::ProtocolControl;
using Reg32 = Registers::System::SyncControl;
using Reg33 = Registers::System::SyncStatus;
using Reg35 = Registers::Attitude::VpeBasicControl;
using Reg36 = Registers::Attitude::VpeMagBasicTuning;
using Reg38 = Registers::Attitude::VpeAccelBasicTuning;
using Reg44 = Registers::HardSoftIronEstimator::RealTimeHsiControl;
using Reg47 = Registers::HardSoftIronEstimator::EstMagCal;
using Reg50 = Registers::VelocityAiding::VelAidingMeas;
using Reg51 = Registers::VelocityAiding::VelAidingControl;
using Reg54 = Registers::IMU::ImuMeas;
using Reg55 = Registers::GNSS::GnssBasicConfig;
using Reg57 = Registers::GNSS::GnssAOffset;
using Reg58 = Registers::GNSS::GnssSolLla;
using Reg59 = Registers::GNSS::GnssSolEcef;
using Reg63 = Registers::INS::InsSolLla;
using Reg64 = Registers::INS::InsSolEcef;
using Reg67 = Registers::INS::InsBasicConfig;
using Reg72 = Registers::INS::InsStateLla;
using Reg73 = Registers::INS::InsStateEcef;
using Reg74 = Registers::INS::FilterStartupBias;
using Reg75 = Registers::System::BinaryOutput1;
using Reg76 = Registers::System::BinaryOutput2;
using Reg77 = Registers::System::BinaryOutput3;
using Reg80 = Registers::IMU::DeltaThetaVelocity;
using Reg82 = Registers::IMU::DeltaThetaVelConfig;
using Reg83 = Registers::WorldMagGravityModel::RefModelConfig;
using Reg84 = Registers::IMU::GyroCal;
using Reg85 = Registers::IMU::ImuFilterControl;
using Reg86 = Registers::GNSSCompass::GnssCompassSignalHealthStatus;
using Reg93 = Registers::GNSSCompass::GnssCompassBaseline;
using Reg97 = Registers::GNSSCompass::GnssCompassEstBaseline;
using Reg98 = Registers::GNSSCompass::GnssCompassStartupStatus;
using Reg99 = Registers::GNSS::GnssSystemConfig;
using Reg100 = Registers::GNSS::GnssSyncConfig;
using Reg101 = Registers::System::NmeaOutput1;
using Reg102 = Registers::System::NmeaOutput2;
using Reg103 = Registers::GNSS::Gnss2SolLla;
using Reg104 = Registers::GNSS::Gnss2SolEcef;
using Reg105 = Registers::INS::InsRefOffset;
using Reg115 = Registers::Heave::HeaveOutputs;
using Reg116 = Registers::Heave::HeaveBasicConfig;
using Reg144 = Registers::INS::InsGnssSelect;
using Reg157 = Registers::GNSS::ExtGnssOffset;
using Reg206 = Registers::System::LegacyCompatibilitySettings;
using Reg239 = Registers::Attitude::YprLinearBodyAccelAngularRates;
using Reg240 = Registers::Attitude::YprLinearInertialAccelAngularRates;
}  // namespace ById

inline bool isConfigEqual(const VN::ConfigurationRegister& lhs, const VN::ConfigurationRegister& rhs)
{
    using namespace VN::Registers::ById;

    if (lhs.id() != rhs.id()) { return false; }

    switch (lhs.id())
    {
        case 0:
            return static_cast<const Reg0&>(lhs) == static_cast<const Reg0&>(rhs);
        case 5:
            return static_cast<const Reg5&>(lhs) == static_cast<const Reg5&>(rhs);
        case 6:
            return static_cast<const Reg6&>(lhs) == static_cast<const Reg6&>(rhs);
        case 7:
            return static_cast<const Reg7&>(lhs) == static_cast<const Reg7&>(rhs);
        case 21:
            return static_cast<const Reg21&>(lhs) == static_cast<const Reg21&>(rhs);
        case 23:
            return static_cast<const Reg23&>(lhs) == static_cast<const Reg23&>(rhs);
        case 25:
            return static_cast<const Reg25&>(lhs) == static_cast<const Reg25&>(rhs);
        case 26:
            return static_cast<const Reg26&>(lhs) == static_cast<const Reg26&>(rhs);
        case 30:
            return static_cast<const Reg30&>(lhs) == static_cast<const Reg30&>(rhs);
        case 32:
            return static_cast<const Reg32&>(lhs) == static_cast<const Reg32&>(rhs);
        case 35:
            return static_cast<const Reg35&>(lhs) == static_cast<const Reg35&>(rhs);
        case 36:
            return static_cast<const Reg36&>(lhs) == static_cast<const Reg36&>(rhs);
        case 38:
            return static_cast<const Reg38&>(lhs) == static_cast<const Reg38&>(rhs);
        case 44:
            return static_cast<const Reg44&>(lhs) == static_cast<const Reg44&>(rhs);
        case 50:
            return static_cast<const Reg50&>(lhs) == static_cast<const Reg50&>(rhs);
        case 51:
            return static_cast<const Reg51&>(lhs) == static_cast<const Reg51&>(rhs);
        case 55:
            return static_cast<const Reg55&>(lhs) == static_cast<const Reg55&>(rhs);
        case 57:
            return static_cast<const Reg57&>(lhs) == static_cast<const Reg57&>(rhs);
        case 67:
            return static_cast<const Reg67&>(lhs) == static_cast<const Reg67&>(rhs);
        case 74:
            return static_cast<const Reg74&>(lhs) == static_cast<const Reg74&>(rhs);
        case 75:
            return static_cast<const VN::Registers::System::BinaryOutput&>(lhs) == static_cast<const VN::Registers::System::BinaryOutput&>(rhs);
        case 76:
            return static_cast<const VN::Registers::System::BinaryOutput&>(lhs) == static_cast<const VN::Registers::System::BinaryOutput&>(rhs);
        case 77:
            return static_cast<const VN::Registers::System::BinaryOutput&>(lhs) == static_cast<const VN::Registers::System::BinaryOutput&>(rhs);
        case 82:
            return static_cast<const Reg82&>(lhs) == static_cast<const Reg82&>(rhs);
        case 83:
            return static_cast<const Reg83&>(lhs) == static_cast<const Reg83&>(rhs);
        case 84:
            return static_cast<const Reg84&>(lhs) == static_cast<const Reg84&>(rhs);
        case 85:
            return static_cast<const Reg85&>(lhs) == static_cast<const Reg85&>(rhs);
        case 93:
            return static_cast<const Reg93&>(lhs) == static_cast<const Reg93&>(rhs);
        case 99:
            return static_cast<const Reg99&>(lhs) == static_cast<const Reg99&>(rhs);
        case 100:
            return static_cast<const Reg100&>(lhs) == static_cast<const Reg100&>(rhs);
        case 101:
            return static_cast<const Reg101&>(lhs) == static_cast<const Reg101&>(rhs);
        case 102:
            return static_cast<const Reg102&>(lhs) == static_cast<const Reg102&>(rhs);
        case 105:
            return static_cast<const Reg105&>(lhs) == static_cast<const Reg105&>(rhs);
        case 116:
            return static_cast<const Reg116&>(lhs) == static_cast<const Reg116&>(rhs);
        case 144:
            return static_cast<const Reg144&>(lhs) == static_cast<const Reg144&>(rhs);
        case 157:
            return static_cast<const Reg157&>(lhs) == static_cast<const Reg157&>(rhs);
        case 206:
            return static_cast<const Reg206&>(lhs) == static_cast<const Reg206&>(rhs);
        default:
            abort();
            break;
    }

    return false;
}  // bool operator==
}  // namespace Registers
}  // namespace VN

#endif  //__REGISTERS_HPP__
