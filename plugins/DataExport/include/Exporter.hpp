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

#ifndef EXPORTER_HEADER_H
#define EXPORTER_HEADER_H

#include <array>
#include <atomic>
#include <cstdio>
#include <HAL/Thread.hpp>

#include "HAL/Thread.hpp"
#include "Implementation/FaPacketProtocol.hpp"
#include "Implementation/Packet.hpp"
#include "Implementation/QueueDefinitions.hpp"

namespace VN
{

class Exporter
{
public:
    Exporter(const size_t& packetCapacity) : _queue{packetCapacity} {}

    virtual ~Exporter() {};

    virtual void addPacketToProcess([[maybe_unused]] const std::shared_ptr<Packet> ptr) = 0;

    virtual void exportToFile() = 0;

    bool start()
    {
        if (_thread == nullptr) { _logging = true; }
        else { return true; }
        _thread = std::make_unique<Thread>(&Exporter::_export, this);
        return false;
    }

    void stop()
    {
        _logging = false;
        _thread->join();
    }

    bool isLogging() const { return _logging; }

    PacketQueue_Interface* getQueuePtr() { return &_queue; }

protected:
    std::atomic<bool> _logging = false;
    std::unique_ptr<Thread> _thread = nullptr;
    PacketQueue<1000> _queue;

private:
    void _export()
    {
        while (_logging || (!_queue.isEmpty()))
        {
            thisThread::sleepFor(1ms);
            exportToFile();
        }
    }
};

template <class T>
std::string extractToString(FaPacketExtractor& extractor, const size_t numToExtract) noexcept
{
    std::string out;
    for (size_t i = 0; i < numToExtract; i++) { out += std::to_string(extractor.extract_unchecked<T>()) + ","; }
    return out;
}

enum Type : uint8_t
{
    U8,
    U16,
    U32,
    U64,
    UTC,
    FLO,
    DUB,
    SAT,
    RAW,
    UNK,
    NON,
};

using VN::Type;

static std::array<std::array<std::pair<Type, uint8_t>, 19>, 7> dataTypes{{
    {std::pair{U64, 1},
     {U64, 1},
     {U64, 1},
     {FLO, 3},
     {FLO, 4},
     {FLO, 3},
     {DUB, 3},
     {FLO, 3},
     {FLO, 3},
     {FLO, 6},
     {FLO, 5},
     {FLO, 7},
     {U16, 1},
     {U32, 1},
     {U64, 1},
     {NON, 0},
     {NON, 0},
     {NON, 0},
     {NON, 0}},  // COMMON
    {std::pair{U64, 1},
     {U64, 1},
     {U64, 1},
     {U16, 1},
     {U64, 1},
     {U64, 1},
     {UTC, 1},
     {U32, 1},
     {U32, 1},
     {U8, 1},
     {NON, 0},
     {NON, 0},
     {NON, 0},
     {NON, 0},
     {NON, 0},
     {NON, 0},
     {NON, 0},
     {NON, 0},
     {NON, 0}},  // TIME
    {std::pair{U16, 1},
     {FLO, 3},
     {FLO, 3},
     {FLO, 3},
     {FLO, 1},
     {FLO, 1},
     {FLO, 4},
     {FLO, 3},
     {FLO, 3},
     {FLO, 3},
     {FLO, 3},
     {U16, 1},
     {NON, 0},
     {NON, 0},
     {NON, 0},
     {NON, 0},
     {NON, 0},
     {NON, 0},
     {NON, 0}},  // IMU
    {std::pair{UTC, 1},
     {U64, 1},
     {U16, 1},
     {U8, 1},
     {U8, 1},
     {DUB, 3},
     {DUB, 3},
     {FLO, 3},
     {FLO, 3},
     {FLO, 3},
     {FLO, 1},
     {FLO, 1},
     {U16, 1},
     {FLO, 7},
     {SAT, 1},
     {NON, 0},
     {RAW, 1},
     {U16, 1},
     {DUB, 1}},  // GNSS
    {std::pair{U16, 1},
     {FLO, 3},
     {FLO, 4},
     {FLO, 9},
     {FLO, 3},
     {FLO, 3},
     {FLO, 3},
     {FLO, 3},
     {FLO, 3},
     {UNK, 12},
     {NON, 0},
     {NON, 0},
     {FLO, 3},
     {FLO, 1},
     {NON, 0},
     {NON, 0},
     {NON, 0},
     {NON, 0},
     {NON, 0}},  // ATTITUDE
    {std::pair{U16, 1},
     {DUB, 3},
     {DUB, 3},
     {FLO, 3},
     {FLO, 3},
     {FLO, 3},
     {FLO, 3},
     {FLO, 3},
     {FLO, 3},
     {FLO, 1},
     {FLO, 1},
     {NON, 0},
     {NON, 0},
     {NON, 0},
     {NON, 0},
     {NON, 0},
     {NON, 0},
     {NON, 0},
     {NON, 0}},  // INS
    {std::pair{UTC, 1},
     {U64, 1},
     {U16, 1},
     {U8, 1},
     {U8, 1},
     {DUB, 3},
     {DUB, 3},
     {FLO, 3},
     {FLO, 3},
     {FLO, 3},
     {FLO, 1},
     {FLO, 1},
     {U16, 1},
     {FLO, 7},
     {SAT, 1},
     {NON, 0},
     {RAW, 1},
     {U16, 1},
     {DUB, 1}},  // GNSS2
}};

}  // namespace VN

#endif  // EXPORTER_HEADER_H
