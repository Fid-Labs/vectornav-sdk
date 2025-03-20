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

#include <string>
#include <cstdint>
#include <functional>
#include <vector>
#include <filesystem>

#include "Interface/Sensor.hpp"
#include "Interface/Registers.hpp"

namespace VN
{
class SensorConfigurator
{
public:
    SensorConfigurator(Sensor& sensor, std::string port) : sensor(sensor), port(port) {};

    bool ConfigureSensor(std::vector<std::unique_ptr<ConfigurationRegister>>& config);

    bool SaveConfiguration(std::filesystem::path path);

    bool LoadConfiguration(std::filesystem::path path);

private:
    Sensor& sensor;
    std::string port;

    static constexpr bool SUCCESS = false;
    static constexpr bool FAIL = true;

    std::unique_ptr<VN::ConfigurationRegister> GetRegisterByIndex(uint32_t idx);

    std::vector<std::unique_ptr<VN::ConfigurationRegister>> registerScan();
};
}  // namespace VN
