/*
 * ObservableOutput.h
 *
 *  Created on: Feb 11, 2013
 *      Author: rovigatti
 */

#ifndef OBSERVABLEOUTPUT_H_
#define OBSERVABLEOUTPUT_H_

#include <vector>
#include <fstream>

#include "BaseObservable.h"

#include <zstd.h>
class ZstdOStream : public std::ostream {
public:
    ZstdOStream(std::ostream& dest, int level = 1)
        : std::ostream(&buffer_), buffer_(dest, level) {}

    ~ZstdOStream() {
        flush();
        buffer_.finish();
    }

private:
    class ZstdStreamBuf : public std::streambuf {
    public:
        ZstdStreamBuf(std::ostream& dest, int level)
            : dest_(dest) {
            cstream_ = ZSTD_createCStream();
            ZSTD_initCStream(cstream_, level);

            outBuffer_.resize(131072); // 128 KB
        }

        ~ZstdStreamBuf() {
            ZSTD_freeCStream(cstream_);
        }

        int overflow(int ch) override {
            if (ch != EOF) {
                inBuffer_.push_back((char)ch);
                if (inBuffer_.size() >= 65536) { // 64 KB chunk
                    compressChunk(false);
                }
            }
            return ch;
        }

        int sync() override {
            compressChunk(false);
            dest_.flush();
            return 0;
        }

        void finish() {
            compressChunk(true);
            dest_.flush();
        }

    private:
        void compressChunk(bool end) {
            ZSTD_inBuffer in{ inBuffer_.data(), inBuffer_.size(), 0 };

            while (in.pos < in.size || end) {
                ZSTD_outBuffer out{ outBuffer_.data(), outBuffer_.size(), 0 };

                size_t ret = end
                    ? ZSTD_endStream(cstream_, &out)
                    : ZSTD_compressStream(cstream_, &out, &in);

                dest_.write((char*)out.dst, out.pos);

                if (!end && in.pos == in.size)
                    break;

                if (ZSTD_isError(ret))
                    throw std::runtime_error(ZSTD_getErrorName(ret));

                if (end && ret == 0)
                    break;
            }
            inBuffer_.clear();
        }

        std::ostream& dest_;
        ZSTD_CStream* cstream_;
        std::vector<char> inBuffer_;
        std::vector<char> outBuffer_;
    };

    ZstdStreamBuf buffer_;
};


/**
 * @brief Manages a single output stream.
 *
 * Instances of this class manage all the output streams associated to the simulation or to the analysis.
 * It stores a list of observables (derived from BaseObservable), whose results are then computed and printed
 * periodically. The constructor requires a string containing all the key=values associated to the instance.
 * The supported syntax is (optional values are between [])
 * @verbatim
 name = <string> (name of the output stream. stdout or stderr are accepted values)
 print_every = <integer> (frequency of output, in number steps for oxDNA, in number of configurations for DNAnalysis)
 [start_from = <integer> (start outputing from the given step, defaults to 0)]
 [stop_at = <integer> (stop outputing at this step, defaults to -1 (which means never))]
 [only_last = <bool> (if true, the output will not be appended to the stream, but it will overwrite the previous output each time, defaults to false)]
 [binary = <bool> (if true, the output will be printed in binary, defaults to false)]
 [linear = <bool> (if true the OutputObservable will save in linear scale, otherwise will use the logline scale by FS. Defaults to true)]
 [update_name_with_time = <bool> (if true the output filename will be changed by using the 'name' key as a prefix and the current step as a suffix. Defaults to false)]
 col_<n> = {\ntype = name of the first observable\n[other observable options as lines of 'key = value']\n} (this syntax specifies the column of the output file. Note that <n> is the column index and should start from 1)
 @endverbatim
 * You can have as many observables (columns) as you want. You can put them in any order, but there has to be a col_1 key.
 * In addition, in order to have, e.g., a col_4, you have to also have col_1, col_2 and col_3.
 *
 * SimBackend stores a list of {@link ObservableOutput ObservableOutputs}. For oxDNA, there are two default
 * streams: one writes to stdout, the other one writes to the energy file. For both oxDNA and DNAnalysis,
 * you can add custom streams by adding them to the input file as options or as command line arguments.
 * For oxDNA, streams can be specified with 'data_output_\<n\>' keys, whereas for DNAnalysis, streams can be specified
 * with 'analysis_data_output_\<n\>' keys. In both cases, \<n\> is an integer starting from 1.
 */
class ObservableOutput {
protected:
	std::vector<ObservablePtr> _obss;
	std::ofstream _output_stream;
	std::ostream *_output;
	llint _print_every;
	bool _is_binary;
	int _append;
	int _only_last;
	llint _start_from;
	llint _stop_at;
	std::string _prefix;
	std::string _output_name;
	std::string _base_name;
	bool _linear;
	int _log_ppc;
	int _log_n0;
	number _log_fact;
	llint _log_next;
	llint _log_tot_cycle;
	int _log_pos_in_cycle;
	int _log_n_cycle;
	bool _update_name_with_time;

	// Zstd compressed stream, so that we can manage its lifetime
	std::unique_ptr<std::ostream> _compressed_stream;
	bool _use_zstd = false;
	int  _zstd_level = 3;

	llint _bytes_written;

	void _open_output();
	void _set_next_log_step();

public:
	/**
	 * @brief Constructor. The two arguments are mandatory (i.e. there are no other constructors in this class)
	 *
	 * @param stream_string a string containing all the key=values lines related to the object and to its associated observables
	 * @param sim_inp simulation input file
	 */
	ObservableOutput(std::string &stream_string);
	virtual ~ObservableOutput();

	/**
	 * @brief Initialize the object
	 */
	void init();

	/**
	 * @brief Add the given observable to this output's list of observables
	 *
	 * @param obs_string
	 */
	void add_observable(ObservablePtr new_obs);

	/**
	 * @brief Add the observable defined by obs_string to the list
	 *
	 * The observable will be appended to the observable list. The output order will reflect the order with which
	 * observable have been added to the list.
	 * @param obs_string A string containing the definition of the observable. The accepted syntax is
	 * type = observable name
	 * [additional, observable-dependent options as lines of 'key = value']
	 */
	void add_observable(std::string obs_string);

	/**
	 * @brief Clear the output of all the observables
	 */
	void clear();

	/**
	 * @brief Prints out the output generated by the stored observables
	 * @param step simulation step
	 */
	void print_output(llint step);

	/**
	 * @brief Changes the output stream. It only supports files (and not other kind of streams such as stdout or stderr)
	 *
	 * @param new_filename
	 */
	void change_output_file(std::string new_filename);

	/**
	 * @brief Checks whether the object is ready to print, and returns the result
	 * @param step simulation step
	 * @return true if the object is ready to print, false otherwise
	 */
	bool is_ready(llint step);

	/**
	 * @brief Returns the number of bytes written to the output file
	 *
	 */
	llint get_bytes_written() {
		return _bytes_written;
	}

	/**
	 * @brief Returns the name of the file being written to as a string
	 *
	 */
	std::string get_output_name() {
		return _output_name;
	}
};

using ObservableOutputPtr = std::shared_ptr<ObservableOutput>;

#endif /* OBSERVABLEOUTPUT_H_ */
