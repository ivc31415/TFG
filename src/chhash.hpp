#pragma once

#include <stdint.h>
#include <string>

#include "q_ary.hpp"

struct _chhash_tree;

class CHHash {
	private:
		_chhash_tree* tree = nullptr;

		//Generate the linear code
		bool generate_linear_code(const uint64_t& quantity, const uint64_t& code_length, const uint64_t& q, matrix<uint64_t>* A, const uint64_t& seed);
		std::string to_c_code(_chhash_tree* branch, const std::string& preamble);
		void delete_tree(_chhash_tree* t = nullptr);

	public:
		q_ary* C = nullptr;
		double time_code = 0, time_chhash = 0;

		//Create a Linear Code with a target number of codewords (quantity). A matrix loaded from the directory data_folder.
		CHHash(const uint64_t& quantity, const uint64_t& q, const std::string& data_folder = "data");

		//Create a Linear Code with a target number of codewords (quantity). It requires the A matrix. Can set up a seed for the uniform sampling of codewords.
		CHHash(const uint64_t& quantity, const uint64_t& q, matrix<uint64_t>* A, const uint64_t& seed = 0);

		//Create a Linear Code with a target number of codewords and a target codeword length. Not well supported
		CHHash(const uint64_t& quantity, const uint64_t& code_length, const uint64_t& q, matrix<uint64_t>* A = nullptr, const uint64_t& seed = 0);

		//Obtain the minimun length of each codeword for the linear code to be a perfect hash family
		static uint64_t get_code_length(const uint64_t& q, const uint64_t& M, const uint64_t& S_A);

		~CHHash();

		uint64_t get(const std::string& key, _chhash_tree* branch = nullptr);

		void generate_tree(const uint64_t& number_keys, std::string* keys);

		void generate_tree(const uint64_t& number_keys, std::string keys_file);

		std::string to_c_code();

		qint* get_codeword(const std::string& key, _chhash_tree* branch = nullptr);
};