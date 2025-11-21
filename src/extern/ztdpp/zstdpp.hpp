#pragma once

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>
#include <cstdint>
#include <fstream>
#include <zstd.h>

namespace zstdpp {

using byte_t = std::uint8_t;
using compress_level_t = std::uint8_t;
using threads_number_t = std::uint8_t;
using buffer_t = std::vector<byte_t>;
using string_t = std::string;
using size_buffer_t = std::size_t;

namespace utils{
    
    inline buffer_t to_bytes(std::string const& str){
        buffer_t bytes{};
        bytes.reserve(str.size());
    
        for (auto const& c : str) {
            bytes.push_back(c);
        }
        return bytes;
    }
    
    inline string_t to_string(buffer_t const& bytes){
        string_t str{};
        str.reserve(bytes.size());
        for (auto const& byte : bytes) {
            str.push_back(byte);
        }
        return str;
    }
    
} // namespace utils

namespace stream{
    
    /* Resources Management Structure */
    struct Resources{
        
        Resources(){
            buffIn.resize(buffInSize);
            buffOut.resize(buffOutSize);
        }
        
        size_t getToRead(){ return buffInSize; }
        size_t getToWrite(){ return buffOutSize; }
        
        size_t readFrom(std::istream& in){
            in.read((char*)buffIn.data(), buffInSize);
            return in.gcount();
        }
        
        void writeTo(std::ostream& out, size_t pos = 0){
            out.write((char*)buffOut.data(), pos);
        }
        
        void* getRawInData(){ return buffIn.data(); }
        void* getRawOutData(){ return buffOut.data(); }
        
        private:
            buffer_t buffIn, buffOut;
            size_t const buffInSize = ZSTD_CStreamInSize();
            size_t const buffOutSize = ZSTD_CStreamOutSize();
    };
    
    struct Context{
        
        /// Default is decompression
        Context(): compress_ctx(NULL), decompress_ctx(ZSTD_createDCtx()) {
            if (decompress_ctx == NULL) {
                throw std::runtime_error("ZSTD_createCCtx() failed!");
            }
        }
        
        /// Compression with specified level
        Context(compress_level_t compress_level, threads_number_t nThreads)
        : compress_ctx(ZSTD_createCCtx()), decompress_ctx(NULL) {
            if (compress_ctx == NULL) {
                throw std::runtime_error("ZSTD_createCCtx() failed!");
            }
            
            /* Set the compression level, and enable the checksum. */
            ZSTD_CCtx_setParameter(compress_ctx, ZSTD_c_compressionLevel, 3);
            ZSTD_CCtx_setParameter(compress_ctx, ZSTD_c_checksumFlag, 1);
            
            /* Config if required workers */
            size_t const r = ZSTD_CCtx_setParameter(compress_ctx, ZSTD_c_nbWorkers, nThreads);
            if (ZSTD_isError(r)) {
                std::cerr << "Note: the linked libzstd library doesn't support multithreading. \n"
                          << "\tReverting to single-thread mode. \n" << std::endl;
            }
        }
        
        ~Context(){
            // Free the context
            if (compress_ctx != NULL){
                ZSTD_freeCCtx(compress_ctx);
            }
            
            if (decompress_ctx != NULL) {
                ZSTD_freeDCtx(decompress_ctx);
            }
        }
        
        size_t operator()(
            ZSTD_inBuffer& in,
            ZSTD_outBuffer& out,
            ZSTD_EndDirective const& mode
        ){
            return ZSTD_compressStream2(
                compress_ctx, 
                &out, 
                &in, 
                mode
            );
        }
        
        size_t operator()(
            ZSTD_inBuffer& in,
            ZSTD_outBuffer& out
        ){ return ZSTD_decompressStream(decompress_ctx, &out , &in); }
        
        private:
            ZSTD_CCtx* const compress_ctx;
            ZSTD_DCtx* const decompress_ctx;
    };
    
    inline void compress(
        std::istream& in, 
        std::ostream& out, 
        threads_number_t nThreads = 1,
        compress_level_t compress_level = 3
    ){
        
        Resources res{};
        Context ctx(compress_level, nThreads);
        
        /* Loop for read chunks & write to output */
        size_t const toRead = res.getToRead();
        while (true) {
            size_t const read = res.readFrom(in);
            auto isLastChunk = read < toRead;
            ZSTD_EndDirective const mode = isLastChunk ? ZSTD_e_end : ZSTD_e_continue;
            
            /* Set the input buffer to what we just read.
            * We compress until the input buffer is empty, each time flushing the
            * output.
            */
            ZSTD_inBuffer input = { res.getRawInData(), read, 0 };
            int finished{0};
            do{
                ZSTD_outBuffer output = { res.getRawOutData(), res.getToWrite(), 0 };
                size_t const remaining = ctx(input,output,mode); // perform compression
                
                res.writeTo(out, output.pos);
                
                /* Verify that the output was written correctly */
                finished = isLastChunk ? (remaining == 0) : (input.pos == input.size);
            }while(!finished);
            
            if (input.pos != input.size){
                throw std::runtime_error("Error: zstd only returns 0 when the input is completely consumed!");
            }
            
            if (isLastChunk) {
                break;
            }
            
        }
        
    }
    
    inline void decompress(
        std::istream& in, 
        std::ostream& out, 
        threads_number_t nThreads = 1
    ){
        Resources res{};
        Context ctx{};
        
        size_t const toRead = res.getToRead();
        size_t read;
        size_t lastRet = 0;
        int isEmpty = 0;
        
        while (!isEmpty) {
            read = res.readFrom(in);
            isEmpty = read == 0;
            
            ZSTD_inBuffer input = { res.getRawInData(), read, 0 };
            
            /* Given a valid frame, zstd won't consume the last byte of the frame
            * until it has flushed all of the decompressed data of the frame.
            * Therefore, instead of checking if the return code is 0, we can
            * decompress just check if input.pos < input.size.
            */
            while (input.pos < input.size) {
                ZSTD_outBuffer output = { res.getRawOutData(), res.getToWrite(), 0 };
                size_t const ret = ctx(input, output); // perform decompression
                
                res.writeTo(out, output.pos);
                lastRet = ret;
            }
        }
        
        if (lastRet != 0) {
            throw std::runtime_error("Error: zstd only returns 0 when the input is completely consumed!");
        }
        
    }
    
    
    
} // namespace stream

/* In-place compression */
namespace inplace{
    
    inline size_buffer_t compress(
        buffer_t const& data,
        buffer_t& buffer,
        compress_level_t compress_level = 3
    ) {
      size_t est_compress_size = ZSTD_compressBound(data.size());
    
      buffer.resize(est_compress_size);
    
      auto compress_size = ZSTD_compress((void*)buffer.data(), est_compress_size,
                                         data.data(), data.size(), compress_level);
    
      buffer.resize(compress_size);
      buffer.shrink_to_fit();
    
      return buffer.size();
    }
    
    inline size_buffer_t decompress(buffer_t &data, buffer_t& out_buffer) {
      auto const est_decomp_size =
          ZSTD_getFrameContentSize(data.data(), data.size());
      out_buffer.resize(est_decomp_size);
    
      size_t const decomp_size = ZSTD_decompress(
          (void*)out_buffer.data(), est_decomp_size, data.data(), data.size());
    
      out_buffer.resize(decomp_size);
      out_buffer.shrink_to_fit();
      return decomp_size;
    }
}

/* Streaming Functions */

inline void stream_compress(
    string_t const& in, 
    string_t const& out, 
    threads_number_t nThreads = 1,
    compress_level_t compress_level = 3
){
    std::ifstream in_file(in, std::ios::binary);
    std::ofstream out_file(out, std::ios::binary);
    stream::compress(in_file, out_file, nThreads, compress_level);
}

inline void stream_decompress(
    string_t const& in, 
    string_t const& out
){
    std::ifstream in_file(in, std::ios::binary);
    std::ofstream out_file(out, std::ios::binary);
    stream::decompress(in_file, out_file);
}

/* Principal functions using inplace functions */

inline buffer_t compress( buffer_t const& data, compress_level_t compress_level = 3 ){
    buffer_t comp_buffer{};
    size_t est_compress_size = inplace::compress(data, comp_buffer);
    return comp_buffer;
}

inline buffer_t decompress(buffer_t &data) {
  buffer_t decomp_buffer{};
  auto const est_decomp_size = inplace::decompress(data, decomp_buffer);
  return decomp_buffer;
}

/* Re-Using compress/decompress functions with buffer_t for receive strings like data */

inline buffer_t compress(const string_t data, compress_level_t compress_level) {
    auto comp_buffer = utils::to_bytes(data);
    return compress(comp_buffer, compress_level);
}

inline buffer_t decompress(string_t data) {
    auto buffer = utils::to_bytes(data);
    return decompress(buffer);
}


} // namespace zstdpp
