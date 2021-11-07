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
#include "colord_api.h"
#include "decompression_common.h"


namespace colord
{
    std::string ReadsSourceToString(ReadsSource readsSource)
    {
        switch (readsSource)
        {
        case colord::ReadsSource::ONT:
            return "Oxford Nanopore";
        case colord::ReadsSource::PBRaw:
            return "PacBio raw";
        case colord::ReadsSource::PBHiFi:
            return "PacBio HiFi";
        default:
            std::ostringstream oss;
            oss << "Unknown data source, this should never happen, please contact authors: " << __FILE__ << ": " << __LINE__;
            throw std::runtime_error(oss.str());
        }
    }

    std::string QualityCompressionModeToString(QualityCompressionMode qualityCompressionMode)
    {
        switch (qualityCompressionMode)
        {
        case colord::QualityCompressionMode::Original:
            return "Original";            
        case colord::QualityCompressionMode::QuinaryAverage:
            return "Quinary average";
        case colord::QualityCompressionMode::QuadAverage:
            return "Quad average";
        case colord::QualityCompressionMode::BinaryAverage:
            return "Binary average";
        case colord::QualityCompressionMode::QuinaryThreshold:
            return "Quinary threshold";
        case colord::QualityCompressionMode::QuadThreshold:
            return "Quad threshold";
        case colord::QualityCompressionMode::BinaryThreshold:
            return "Binary threshold";
        case colord::QualityCompressionMode::Average:
            return "Average";
        case colord::QualityCompressionMode::None:
            return "None";
        default:
            std::ostringstream oss;
            oss << "Unknown quality compression mode, this should never happen, please contact authors: " << __FILE__ << ": " << __LINE__;
            throw std::runtime_error(oss.str());
        }
    }

    std::string HeaderCompressionModeToString(HeaderCompressionMode headerCompressionMode)
    {
        switch (headerCompressionMode)
        {
        case colord::HeaderCompressionMode::Original:
            return "Original";
        case colord::HeaderCompressionMode::Main:
            return "Main";
        case colord::HeaderCompressionMode::None:
            return "None";
        default:
            std::ostringstream oss;
            oss << "Unknown header compression mode, this should never happen, please contact authors: " << __FILE__ << ": " << __LINE__;
            throw std::runtime_error(oss.str());
        }
    }

    void Info::ToOstream(std::ostream& oss) const
    {
        oss << "is fastq: " << std::boolalpha << isFastq << "\n";
        oss << "colord archive version: " << versionMajor << "." << versionMinor << "." << versionPatch << "\n";
        oss << "total reads: " << totalReads << "\n";
        time_t time = this->time;
        oss << "colord archive creaton datetime: " << asctime(localtime(&time));
        oss << "command line used to create colord archive: " << fullCommandLine << "\n";
        oss << "compression level: " << compressionLevel << "\n";
        oss << "reads source: " << ReadsSourceToString(readsSource) << "\n";
        oss << "quality compression mode: " << QualityCompressionModeToString(qualityCompressionMode) << "\n";
        oss << "header compression mode: " << HeaderCompressionModeToString(headerCompressionMode) << "\n";

        oss << "quality reverse thresholds: ";
        for (auto v : qualityReverseThresholds)
            oss << v << " ";
        oss << "\n";
    }

    class CNullLogger : public ILogger
    {
    public:
        virtual void Log(const std::string& msg) override
        {
            
        }
    };

    class CExceptonErrorHandler : public IErrorHandler
    {
    public:
        virtual void LogError(const std::string& msg) override
        {
            throw std::runtime_error(msg);            
        }
    };

    class CToAPIDecompressedStreamConsumer : public IDecompressedStreamConsumer
    {        
        CParallelQueue<decomp_read_pack_t>* read_decompr_queue{};
        CParallelQueue<header_pack_t>* header_decompr_queue{};
        CParallelQueue<decomp_qual_pack_t>* qual_decompr_queue{};

        decomp_read_pack_t read_pack;
        read_t read;
        uint32_t read_pack_pos{};

        header_pack_t header_pack;
        header_elem_t header;
        uint32_t header_pack_pos{};

        decomp_qual_pack_t qual_pack;
        qual_t qual;
        uint32_t qual_pack_pos{};
        
        template<typename T_ELEM, typename PACK_T>
        bool nextImpl(T_ELEM& elem, CParallelQueue<PACK_T>& queue, PACK_T& pack, uint32_t& pos)
        {
            if (pos == pack.size())
            {
                pos = 0;
                if (!queue.Pop(pack))
                    return false;
            }
            elem = std::move(pack[pos++]);
            return true;
        }
        bool nextRead(read_t& read)
        {
            return nextImpl(read, *read_decompr_queue, read_pack, read_pack_pos);
        }
        bool nextQual(qual_t& qual)
        {
            return nextImpl(qual, *qual_decompr_queue, qual_pack, qual_pack_pos);
        }
        bool nextHeader(header_elem_t& header)
        {
            return nextImpl(header, *header_decompr_queue, header_pack, header_pack_pos);
        }
    public:        
        virtual void ConsumeFasta(CParallelQueue<decomp_read_pack_t>& read_decompr_queue, CParallelQueue<header_pack_t>& header_decompr_queue) override
        {
            this->read_decompr_queue = &read_decompr_queue;
            this->header_decompr_queue = &header_decompr_queue;
        }
        virtual void ConsumeFastq(CParallelQueue<decomp_read_pack_t>& read_decompr_queue, CParallelQueue<header_pack_t>& header_decompr_queue, CParallelQueue<decomp_qual_pack_t>& qual_decompr_queue) override
        {
            this->read_decompr_queue = &read_decompr_queue;
            this->header_decompr_queue = &header_decompr_queue;
            this->qual_decompr_queue = &qual_decompr_queue;
        }

        DecompressionRecord NextRecordFasta()
        {
            DecompressionRecord res;

            bool hasNextRead = nextRead(read);
            bool hasNextHeader = nextHeader(header);
            
            if (hasNextRead && hasNextHeader)
            {
                res.read_header = std::move(header.first);
                
                res.read.reserve(read.size());                
                for (uint32_t i = 0; i < read.size() - 1; ++i)
                    res.read.push_back("ACGTN"[read[i]]);
                return res;
            }
            if (hasNextRead + hasNextHeader != 0)
            {
                std::ostringstream oss;
                oss << "cirtical, contact authors, file: " << __FILE__ << ", line: " << __LINE__;
                throw std::runtime_error(oss.str());
            }
            res.finished = true;
            return res;
        }

        DecompressionRecord NextRecordFastq()
        {
            DecompressionRecord res;

            bool hasNextRead = nextRead(read);
            bool hasNextQual = nextQual(qual);
            bool hasNextHeader = nextHeader(header);
            
            while (hasNextRead && hasNextQual && hasNextHeader)
            {
                res.read_header = std::move(header.first);

                res.read.reserve(read.size());

                res.read.reserve(read.size());
                for (uint32_t i = 0; i < read.size() - 1; ++i)
                    res.read.push_back("ACGTN"[read[i]]);
                
                if (header.second == qual_header_type::eq_read_header)
                    res.qual_header = res.read_header;
                
                res.qual.reserve(qual.size());
                for (uint32_t i = 0; i < qual.size(); ++i)
                    res.qual.push_back(static_cast<char>(qual[i]));

                return res;
            }
            if (hasNextRead + hasNextQual + hasNextHeader != 0)
            {
                std::ostringstream oss;
                oss << "cirtical, contact authors, file: " << __FILE__ << ", line: " << __LINE__;
                throw std::runtime_error(oss.str());                
            }
            res.finished = true;
            return res;
        }
        
        DecompressionRecord NextRecord()
        {
            if (qual_decompr_queue)
                return NextRecordFastq();
            return NextRecordFasta();
        }
    };

    class DecompressionStream::DecompressionStreamImpl
    {
        CNullLogger null_logger;
        CExceptonErrorHandler exception_error_handler;
        CDecmpressionModule decompression_module;
        CToAPIDecompressedStreamConsumer consumer;
        Info info;

        void readInfo()
        {
            info.isFastq = decompression_module.GetIsFastq();

            info.versionMajor = decompression_module.GetInfo().version_major;
            info.versionMinor = decompression_module.GetInfo().version_minor;
            info.versionPatch = decompression_module.GetInfo().version_patch;
            info.totalBytes = decompression_module.GetInfo().total_bytes;
            info.totalBases = decompression_module.GetInfo().total_bases;
            info.totalReads = decompression_module.GetInfo().total_reads;
            info.time = decompression_module.GetInfo().time;
            info.fullCommandLine = decompression_module.GetInfo().full_command_line;

            info.compressionLevel = decompression_module.GetMetaData().compressionLevel;
            switch (decompression_module.GetMetaData().dataSource)
            {
            case DataSource::ONT:
                info.readsSource = ReadsSource::ONT;
                break;
            case DataSource::PBHiFi:
                info.readsSource = ReadsSource::PBHiFi;
                break;
            case DataSource::PBRaw:
                info.readsSource = ReadsSource::PBRaw;
                break;
            default:
                std::ostringstream oss;
                oss << "Unknown data source, this should never happen, please contact authors: " << __FILE__ << ": " << __LINE__;
                throw std::runtime_error(oss.str());
            }

            switch (decompression_module.GetMetaData().headerComprMode)
            {
            case HeaderComprMode::Original:
                info.headerCompressionMode = HeaderCompressionMode::Original;
                break;
            case HeaderComprMode::Main:
                info.headerCompressionMode = HeaderCompressionMode::Main;
                break;
            case HeaderComprMode::None:
                info.headerCompressionMode = HeaderCompressionMode::None;
                break;
            default:
                std::ostringstream oss;
                oss << "Unknown header compression mode, this should never happen, please contact authors: " << __FILE__ << ": " << __LINE__;
                throw std::runtime_error(oss.str());
            }

            switch (decompression_module.GetMetaData().qualityComprMode)
            {
            case QualityComprMode::Original:
                info.qualityCompressionMode = QualityCompressionMode::Original;
                break;
            case QualityComprMode::QuinaryAverage:
                info.qualityCompressionMode = QualityCompressionMode::QuinaryAverage;
                break;
            case QualityComprMode::QuadAverage:
                info.qualityCompressionMode = QualityCompressionMode::QuadAverage;
                break;
            case QualityComprMode::BinaryAverage:
                info.qualityCompressionMode = QualityCompressionMode::BinaryAverage;
                break;
            case QualityComprMode::QuinaryThreshold:
                info.qualityCompressionMode = QualityCompressionMode::QuinaryThreshold;
                break;
            case QualityComprMode::QuadThreshold:
                info.qualityCompressionMode = QualityCompressionMode::QuadThreshold;
                break;
            case QualityComprMode::BinaryThreshold:
                info.qualityCompressionMode = QualityCompressionMode::BinaryThreshold;
                break;
            case QualityComprMode::Average:
                info.qualityCompressionMode = QualityCompressionMode::Average;
                break;
            case QualityComprMode::None:
                info.qualityCompressionMode = QualityCompressionMode::None;
                break;
            default:
                std::ostringstream oss;
                oss << "Unknown quality compression mode, this should never happen, please contact authors: " << __FILE__ << ": " << __LINE__;
                throw std::runtime_error(oss.str());
            }

            info.qualityReverseThresholds = decompression_module.GetQualityRevThresholds();
        }
    public:
        DecompressionStreamImpl(const std::string& inputFilePath, const std::string& refGenomePath) :
            decompression_module(inputFilePath, refGenomePath, null_logger, false, exception_error_handler, consumer, true)
        {
            decompression_module.Run();
            
            readInfo();
        }

        Info GetInfo() const
        {
            return info;
        }

        DecompressionRecord NextRecord()
        {
            return consumer.NextRecord();
        }

        ~DecompressionStreamImpl()
        {
            decompression_module.Cancel();
        }
    };

    DecompressionStream::DecompressionStream(const std::string& inputFilePath) 
        : DecompressionStream(inputFilePath, "")
    {

    }
    DecompressionStream::DecompressionStream(const std::string& inputFilePath, const std::string& refGenomePath)
        :
        pImpl(std::make_unique<DecompressionStreamImpl>(inputFilePath, refGenomePath))
    {

    }
    DecompressionStream::~DecompressionStream() = default;

    Info DecompressionStream::GetInfo() const
    {
        return pImpl->GetInfo();
    }
    DecompressionRecord DecompressionStream::NextRecord()
    {
        return pImpl->NextRecord();
    }
}