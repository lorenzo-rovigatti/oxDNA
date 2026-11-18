/*
 * CompressedFileHeader.h
 *
 * The 8-byte header of zstd-compressed oxDNA files (trajectories, configurations and observable outputs):
 *
 *   bytes 0-3: magic 'O', 'X', 'D', 0x01
 *   byte  4  : header version
 *   byte  5  : zstd compression level
 *   bytes 6-7: reserved (0)
 *
 * followed by one zstd frame per configuration or output. The format is defined here only: both the code that writes
 * these files and the code that reads them use these definitions. Python tools that read them (oat's RyeReader,
 * utils/decompress_zstd.py and utils/fix_zstd_trajectory.py) have their own copies of the magic and version.
 */

#ifndef COMPRESSEDFILEHEADER_H_
#define COMPRESSEDFILEHEADER_H_

#include <array>
#include <cstddef>
#include <cstdint>

namespace CompressedFileHeader {

constexpr std::size_t SIZE = 8;
constexpr std::size_t MAGIC_SIZE = 4;
constexpr std::array<std::uint8_t, MAGIC_SIZE> MAGIC = { { 'O', 'X', 'D', 0x01 } };
/// the version written by this code, and the only one it can read
constexpr std::uint8_t VERSION = 0x01;

using Header = std::array<std::uint8_t, SIZE>;

/// the header of a new file compressed with the given zstd level
inline Header make(int compression_level) {
	return Header { { MAGIC[0], MAGIC[1], MAGIC[2], MAGIC[3], VERSION, static_cast<std::uint8_t>(compression_level), 0x00, 0x00 } };
}

/// true if the first n bytes of data (n may be smaller than SIZE) start with the magic
inline bool has_magic(const void *data, std::size_t n) {
	if(n < MAGIC_SIZE) {
		return false;
	}
	const std::uint8_t *bytes = static_cast<const std::uint8_t *>(data);
	for(std::size_t i = 0; i < MAGIC_SIZE; i++) {
		if(bytes[i] != MAGIC[i]) {
			return false;
		}
	}
	return true;
}

/// the version byte of a complete header
inline std::uint8_t version(const void *header) {
	return static_cast<const std::uint8_t *>(header)[4];
}

/// the compression level stored in a complete header
inline int compression_level(const void *header) {
	return static_cast<const std::uint8_t *>(header)[5];
}

} // namespace CompressedFileHeader

#endif /* COMPRESSEDFILEHEADER_H_ */
