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

#ifndef FILES_EXPORT_ASCII_HEADER_HPP
#define FILES_EXPORT_ASCII_HEADER_HPP

#include <vector>
#include <memory>

#include "Exporter.hpp"
#include "HAL/File.hpp"
#include "Implementation/Packet.hpp"
#include "Interface/Command.hpp"

namespace VN
{

struct AsciiInfo
{
    const AsciiHeader header;
    OutputFile file;
};

class ExporterAscii : public Exporter
{
public:
    ExporterAscii(const Filesystem::FilePath& outputDir) : Exporter(2048), _filePath(outputDir)
    {
        if (!_filePath.empty() && _filePath.to_string().back() != std::filesystem::path::preferred_separator)
        {
            _filePath = _filePath + std::filesystem::path::preferred_separator;
        }
    }

    ~ExporterAscii() {}

    void addPacketToProcess([[maybe_unused]] const std::shared_ptr<Packet> ptr) override
    {
        // _queue.push_back(ptr);
    }

    void exportToFile() override
    {
        while (!_queue.isEmpty())
        {
            const auto p = _queue.get();
            if (!p || (p->details.syncByte != PacketDetails::SyncByte::Ascii)) { return; }

            OutputFile& ascii = getFileHandle(p->details.asciiMetadata.header);

            ascii.write(reinterpret_cast<const char*>(p->buffer), p->details.asciiMetadata.length);
        }
    }

private:
    Filesystem::FilePath _filePath;

    std::vector<AsciiInfo> _asciiInfo;

    OutputFile& getFileHandle(const AsciiHeader& header)
    {
        for (auto& tmp : _asciiInfo)
        {
            if (tmp.header == header) { return tmp.file; }
        }

        // if we don't find the header we need to init a new text file
        AsciiMessage fileName = _filePath + header + ".txt";

        _asciiInfo.emplace_back(AsciiInfo{header, OutputFile(fileName)});

        return _asciiInfo.back().file;
    }
};

}  // namespace VN

#endif  // FILES_EXPORT_ASCII_HEADER_HPP
