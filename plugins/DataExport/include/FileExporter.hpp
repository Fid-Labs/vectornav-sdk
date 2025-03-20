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

#ifndef FILE_EXPORTER_HPP
#define FILE_EXPORTER_HPP

#include <memory>
#include <filesystem>

#include "Exporter.hpp"
#include "HAL/File.hpp"
#include "TemplateLibrary/Vector.hpp"
#include "Interface/Registers.hpp"
#include "Implementation/PacketSynchronizer.hpp"
#include "Implementation/AsciiPacketDispatcher.hpp"
#include "Implementation/FaPacketDispatcher.hpp"
#include "Implementation/FbPacketDispatcher.hpp"

namespace VN
{

constexpr uint8_t EXPORTER_CAPACITY = 5;

class FileExporter
{
public:
    FileExporter() = default;

    bool processFile(const Filesystem::FilePath& fileName)
    {
        const auto fileSizeInBytes = std::filesystem::file_size(fileName.c_str());
        auto buffer = std::make_unique<uint8_t[]>(fileSizeInBytes);

        InputFile inputFile(fileName);

        inputFile.read(reinterpret_cast<char*>(buffer.get()), fileSizeInBytes);

        ByteBuffer byteBuffer(buffer.get(), fileSizeInBytes, fileSizeInBytes);

        PacketSynchronizer packetSynchronizer(byteBuffer);
        packetSynchronizer.addDispatcher(&_asciiPacketDispatcher);
        packetSynchronizer.addDispatcher(&_faPacketDispatcher);
        packetSynchronizer.addDispatcher(&_fbPacketDispatcher);

        Registers::System::BinaryOutputMeasurements bor;

        for (auto& e : _exporters)
        {
            _faPacketDispatcher.addSubscriber(e->getQueuePtr(), bor.toBinaryHeader().toMeasurementHeader(), FaPacketDispatcher::SubscriberFilterType::AnyMatch);
            _asciiPacketDispatcher.addSubscriber(e->getQueuePtr(), "VN", AsciiPacketDispatcher::SubscriberFilterType::StartsWith);

            e->start();
        }

        while (!packetSynchronizer.dispatchNextPacket()) {};

        for (auto& e : _exporters) { e->stop(); }
        return false;
    }

    bool addExporter(std::unique_ptr<Exporter>&& exporterToAdd) { return _exporters.push_back(std::move(exporterToAdd)); };

private:
    MeasurementQueue _measurementQueue{Config::PacketDispatchers::compositeDataQueueCapacity};

    CommandProcessor _commandProcessor{[]([[maybe_unused]] AsyncError&& error) {}};

    AsciiPacketDispatcher _asciiPacketDispatcher{&_measurementQueue, Config::PacketDispatchers::cdEnabledMeasTypes, &_commandProcessor};

    FaPacketDispatcher _faPacketDispatcher{&_measurementQueue, Config::PacketDispatchers::cdEnabledMeasTypes};
    FbPacketDispatcher _fbPacketDispatcher{&_faPacketDispatcher, Config::PacketFinders::fbBufferCapacity};

    Vector<std::unique_ptr<Exporter>, EXPORTER_CAPACITY> _exporters;
};

}  // namespace VN

#endif  // FILE_EXPORTER_HPP
