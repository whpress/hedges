#include <cstddef>
#include <iostream>
#include <string>

#include "schifra/schifra_galois_field.hpp"
#include "schifra/schifra_galois_field_polynomial.hpp"
#include "schifra/schifra_sequential_root_generator_polynomial_creator.hpp"
#include "schifra/schifra_reed_solomon_encoder.hpp"
#include "schifra/schifra_reed_solomon_decoder.hpp"
#include "schifra/schifra_reed_solomon_block.hpp"
#include "schifra/schifra_error_processes.hpp"

union b8 {
	unsigned long long w;
	unsigned char b[8];
	b8() {}
	b8(long long ww) { w = ww; }
	b8(const VecUchar &vec) {
		if (vec.size() != 8) throw("bad conversion of VecUchar to b8");
		for (int i = 0; i < 8; i++) b[i] = vec[i];
	}
	inline bool operator==(const b8 &b) { return w == b.w; }
	inline bool operator!=(const b8 &b) { return w != b.w; }
	inline bool operator>(const b8 &b) { return w > b.w; }
	inline bool operator<(const b8 &b) { return w < b.w; }
	inline bool operator>=(const b8 &b) { return w >= b.w; }
	inline bool operator<=(const b8 &b) { return w <= b.w; }
};


using namespace schifra;
template<size_t code_length, size_t fec_length>
struct SchifraCode { // see Schifra's example08

	static const size_t data_length = code_length - fec_length;
	static const size_t generator_polynomial_root_count = fec_length;
	static const size_t field_descriptor = 8;
	static const size_t generator_polynomial_index = 120;

	typedef reed_solomon::shortened_encoder<code_length, fec_length, data_length> encoder_t;
	typedef reed_solomon::shortened_decoder<code_length, fec_length, data_length> decoder_t;
	typedef const galois::field field_t;
	typedef galois::field_polynomial field_polynomial_t;
	typedef reed_solomon::block<code_length, fec_length> block_t;

	encoder_t* encoder;
	decoder_t* decoder;
	field_t* field;
	field_polynomial_t* generator_polynomial;
	block_t block;

	SchifraCode() {
		/* Instantiate Finite Field, Generator Polynomials, Encoder, Decoder */
		const Uint pp6size = 9; // see schifra_galois_field.hpp
		const Uint pp6[] = { 1, 1, 1, 0, 0, 0, 0, 1, 1 };
		field = new field_t(field_descriptor, pp6size, pp6);
		generator_polynomial = new field_polynomial_t(*field);
		if (!make_sequential_root_generator_polynomial(*field,
			generator_polynomial_index,
			generator_polynomial_root_count,
			*generator_polynomial)
			) {
			throw("Error - Failed to create sequential root generator!");
		}
		encoder = new encoder_t(*field, *generator_polynomial);
		decoder = new decoder_t(*field, generator_polynomial_index);
	}
	VecUchar encode(VecUchar &mess) {
		string message(code_length, 0x00);
		for (int i = 0; i < data_length; i++) message[i] = mess[i];
		if (!encoder->encode(message, block)) {
			cout << "Error - Critical encoding failure! "
				<< "Msg: " << block.error_as_string() << endl;
			exit(-1);
		}
		VecUchar codeword(code_length);
		for (int i = 0; i < code_length; i++) codeword[i] = static_cast<Uchar>(block.data[i]);
		return codeword;
	}
	b8 encode(b8 mess) {
		string message(code_length, 0x00);
		for (int i = 0; i < data_length; i++) message[i] = mess.b[i];
		if (!encoder->encode(message, block)) {
			cout << "Error - Critical encoding failure! "
				<< "Msg: " << block.error_as_string() << endl;
			exit(-1);
		}
		b8 codeword;
		for (int i = 0; i < code_length; i++) codeword.b[i] = static_cast<Uchar>(block.data[i]);
		return codeword;
	}
	VecUchar decode(VecUchar &codeword, Int &errs_detected, Int &errs_corrected,
		Int &err_number, bool &recoverable) {
		block.reset(); // not needed?
		for (int i = 0; i < code_length; i++) block.data[i] =
			static_cast<galois::field_symbol>(codeword[i]);
		decoder->decode(block);
		err_number = block.error;
		VecUchar message(code_length);
		for (int i = 0; i < code_length; i++) message[i] = static_cast<Uchar>(block.data[i]);
		errs_detected = Int(block.errors_detected);
		errs_corrected = Int(block.errors_corrected);
		recoverable = !block.unrecoverable;
		return message;
	}
	VecUchar decode(VecUchar &codeword, VecInt &erasures, Int &errs_detected, Int &errs_corrected,
		Int &err_number, bool &recoverable) {
		Int esize = erasures.size();
		std::vector<std::size_t> erasures_long(esize); // in Windows must be <Ullong>  why?
		for (int i = 0; i < esize; i++) erasures_long[i] = std::size_t(erasures[i]);
		block.reset(); // not needed?
		for (int i = 0; i < code_length; i++) block.data[i] =
			static_cast<galois::field_symbol>(codeword[i]);
		decoder->decode(block,erasures_long);
		err_number = block.error;
		VecUchar message(code_length);
		for (int i = 0; i < code_length; i++) message[i] = static_cast<Uchar>(block.data[i]);
		errs_detected = Int(block.errors_detected);
		errs_corrected = Int(block.errors_corrected);
		recoverable = !block.unrecoverable;
		return message;
	}
	b8 decode(b8 codeword, Int &errs_detected, Int &errs_corrected,
		Int &err_number, bool &recoverable) {
		block.reset(); // not needed?
		for (int i = 0; i < code_length; i++) block.data[i] =
			static_cast<galois::field_symbol>(codeword.b[i]);
		decoder->decode(block);
		err_number = block.error;
		b8 message;
		for (int i = 0; i < code_length; i++) message.b[i] = static_cast<Uchar>(block.data[i]);
		errs_detected = Int(block.errors_detected);
		errs_corrected = Int(block.errors_corrected);
		recoverable = !block.unrecoverable;
		return message;
	}
	b8 decode(b8 codeword, vector<Ullong> &erasures, Int &errs_detected, Int &errs_corrected,
		Int &err_number, bool &recoverable) {
		block.reset(); // not needed?
		for (int i = 0; i < code_length; i++) block.data[i] =
			static_cast<galois::field_symbol>(codeword.b[i]);
		decoder->decode(block,erasures);
		err_number = block.error;
		b8 message;
		for (int i = 0; i < code_length; i++) message.b[i] = static_cast<Uchar>(block.data[i]);
		errs_detected = Int(block.errors_detected);
		errs_corrected = Int(block.errors_corrected);
		recoverable = !block.unrecoverable;
		return message;
	}

	/*
	0 "No Error";
	1 "Invalid Encoder";
	2 "Incompatible Generator Polynomial";
	3 "Invalid Decoder";
	4 "Decoder Failure - Non-zero Syndrome";
	5 "Decoder Failure - Too Many Errors/Erasures";
	6 "Decoder Failure - Invalid Symbol Correction";
	7 "Decoder Failure - Invalid Codeword Correction";
	*/
	Int code_len() { return code_length; }
	Int fec_len() { return fec_length; }
	Int data_len() { return code_length - fec_length; }
};
