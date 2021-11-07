/*******************************************************************************

 CoLoRd
 Copyright (C) 2021, M. Kokot, S. Deorowicz, and A. Gudys
 https://github.com/refresh-bio/CoLoRd

 This program is free software: you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation; either version 3 of the License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 A PARTICULAR PURPOSE. See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with this
 program. If not, see https://www.gnu.org/licenses/.

******************************************************************************/
#pragma once
#include <exception>
#include <string>
#include <memory>
#include <vector>
#include <ctime>
#include <ostream>

namespace colord
{
    class DecompressionRecord
    {
        friend class CToAPIDecompressedStreamConsumer;
        bool finished = false;
        std::string read_header, read, qual_header, qual;
    public:
        operator bool() const
        {
            return !finished;
        }
        const std::string& ReadHeader() const
        {
            return read_header;
        }
        const std::string& Read() const
        {
            return read;
        }
        const std::string& QualHeader() const
        {
            return qual_header;
        }
        const std::string& Qual() const
        {
            return qual;
        }
    };
    enum class ReadsSource { ONT, PBRaw, PBHiFi };    
    enum class QualityCompressionMode {
        Original,
        QuinaryAverage,
        QuadAverage,
        BinaryAverage,
        QuinaryThreshold,
        QuadThreshold,
        BinaryThreshold,
        Average,
        None
    };

    enum class HeaderCompressionMode { Original, Main, None };

    struct Info
    {
        bool isFastq;
        uint32_t versionMajor;
        uint32_t versionMinor;
        uint32_t versionPatch;
        uint64_t totalBytes;
        uint64_t totalBases;
        uint32_t totalReads;
        uint64_t time;
        std::string fullCommandLine;

        int32_t compressionLevel{};
        ReadsSource readsSource{};
        QualityCompressionMode qualityCompressionMode{};
        HeaderCompressionMode headerCompressionMode;

        std::vector<uint32_t> qualityReverseThresholds;

        void ToOstream(std::ostream& oss) const;
    };
    class DecompressionStream
    {
        class DecompressionStreamImpl;
        const std::unique_ptr<DecompressionStreamImpl> pImpl;
    public:
        explicit DecompressionStream(const std::string& inputFilePath);
        explicit DecompressionStream(const std::string& inputFilePath, const std::string& refGenomePath);
        Info GetInfo() const;
        DecompressionRecord NextRecord();
        ~DecompressionStream();
    };
}
