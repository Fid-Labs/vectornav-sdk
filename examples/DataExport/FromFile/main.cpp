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

#include <memory>
#include <filesystem>

#include "FileExporter.hpp"
#include "ExporterCsv.hpp"
#include "ExporterAscii.hpp"

using namespace VN;
namespace fs = std::filesystem;

FileExporter fileExporter;
std::string outputDirectory = fs::path(__FILE__).parent_path().string();
std::string filePath = (fs::path(__FILE__).parent_path() / "DataExportFromFile.bin").string();

int main([[maybe_unused]] int argc, [[maybe_unused]] char* argv[])
{
    fileExporter.addExporter(std::make_unique<ExporterCsv>(outputDirectory));
    fileExporter.addExporter(std::make_unique<ExporterAscii>(outputDirectory));

    if (fileExporter.processFile(filePath))
    {
        std::cout << "Error: File processing failed." << std::endl;
        return 1;
    }
    std::cout << "FileExport example complete." << std::endl;
}
