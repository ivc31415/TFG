#include "chhash.hpp"
#include "chrono2.hpp"
#include "lemmas.hpp"
#include "util.hpp"
#include <iostream>
#include <fstream>

struct _chhash_tree;
struct _chhash_tree {
	char c;
	uint64_t i;
	_chhash_tree* yes;
	_chhash_tree* no;
	uint64_t h;
};

//Returns true if it hasn't failed
bool CHHash::generate_linear_code(const uint64_t& quantity, const uint64_t& code_length, const uint64_t& q, matrix<uint64_t>* A, const uint64_t& seed) {
	if(!(q == 4 || q == 5 || q == 7)) {
		throw std::runtime_error("Invalid q parameter. Supported: 4, 5, 7, 8");
	}

	simple_clock c;
	c.reset();

	//Generate C2 code
	q_ary* C2;
	if(q == 5 || q == 7) {
		C2 = lemma_3_12(q);
	} else {
		C2 = lemma_3_14_case_1(q == 4 ? 2 : 3, 3);
	}

	//Generate C1 code and concatenate it with C2
	C = generate_qary_lemma_2_3(code_length, quantity, C2, A, seed);

	delete C2; //Not needed any more

	time_code = c.microseconds()/1000.0;

	return true;
}

_chhash_tree* init_tree(const uint64_t& number_keys, std::string* keys, const uint64_t& i, uint64_t* hash_number, bool* mask) {
	bool* mask_yes = make_array<bool>(number_keys, false);
	bool* mask_no = make_array<bool>(number_keys, false);

	//std::cout << "i: " << i << "\n";
	//std::cout << "hash_number: " << *hash_number << "\n";
	//if(mask != nullptr) {std::cout << "mask: "; _print_array(mask, number_keys); std::cout << "\n";}

	char c = 0;
	bool c_once = true;
	bool has_second = false;
	for(uint64_t j = 0; j < number_keys; j++) {
		if(mask != nullptr && !mask[j])	continue;

		std::string* key = &(keys[j]);

		if(c == 0 && key->length() <= i) { //Special case
			free(mask_no);
			mask_no = copy_array(number_keys, mask);
			mask_no[j] = false;
			_chhash_tree* t = new _chhash_tree{0, i, nullptr, init_tree(number_keys, keys, i, hash_number, mask_no), (*hash_number)++};
			free(mask_no);
			free(mask_yes);
			return t;
		}

		if(c == 0 && key->length() > i) {
			c = key->at(i);
			mask_yes[j] = true;
			//std::cout << "[debug 1] c = '" << c << "', j = " << j << "\n";
		} else if(c != 0 && key->length() > i) {
			if(c == key->at(i)) {
				c_once = false;
				mask_yes[j] = true;
				//std::cout << "[debug 2] c = '" << key->at(i) << "', j = " << j << "\n";
			} else {
				has_second = true;
				mask_no[j] = true;
				//std::cout << "[debug 3] c = '" << key->at(i) << "', j = " << j << "\n";
			}
		}
	}
	//std::cout << "c_once: " << c_once << "\n";
	//std::cout << "mask_yes: "; _print_array(mask_yes, number_keys); std::cout << "\n";
	//std::cout << "mask_no:  "; _print_array(mask_no, number_keys);  std::cout << "\n";

	_chhash_tree* t = new _chhash_tree{c, i, nullptr, nullptr, 0};
	if(c_once) {
		//std::cout << "Adding at i=" << i << " w/ character '" << c << "' to hash " << (*hash_number) << "...\n";
		t->h = (*hash_number)++;
	} else {
		if(has_second) {
			//std::cout << "Branching at i=" << i << " w/ character '" << c << "'...\n";
			t->yes = init_tree(number_keys, keys, i + 1, hash_number, mask_yes);
		} else {
			//std::cout << "Continuing at i=" << i << " w/ character '" << c << "' and no branching...\n";
			delete t;
			t = init_tree(number_keys, keys, i + 1, hash_number, mask_yes);
		}
	}

	if(has_second) {
		//std::cout << "Secondary branching at i=" << i << " w/ character '" << c << "'...\n";
		t->no = init_tree(number_keys, keys, i, hash_number, mask_no);
	}

	free(mask_yes);
	free(mask_no);

	return t;
}

CHHash::CHHash(const uint64_t& quantity, const uint64_t& code_length, const uint64_t& q, matrix<uint64_t>* A, const uint64_t& seed) {
	// Generate Code
	generate_linear_code(quantity * 3, (uint64_t)ceil(code_length / 4.0), q, A, seed);
}

CHHash::CHHash(const uint64_t& quantity, const uint64_t& q, matrix<uint64_t>* A, const uint64_t& seed) {
	//Get shortest code length to fulfil requirements
	uint64_t code_length = get_code_length(q, quantity * 3, A->width) * 4;

	// Generate Code
	generate_linear_code(quantity * 3, code_length, q, A, seed);
}

CHHash::CHHash(const uint64_t& quantity, const uint64_t& q, const std::string& data_folder) {
	//Get precalculated A Matrix
	matrix<uint64_t> A(data_folder + "/A_q_" + std::to_string(q) + ".dat");

	//Get shortest code length to fulfil requirements
	uint64_t code_length = get_code_length(q, quantity * 3, A.width) * 4;

	// Generate Code
	generate_linear_code(quantity * 3, code_length, q, &A, 0);
}

uint64_t CHHash::get_code_length(const uint64_t& q, const uint64_t& M, const uint64_t& S_A) {
	double a = (1 - (factorial(q) * (double)S_A)/(double)pow(q * q, q));
	double b = (double)M/(2 * q * combinatorialD(M, q));
	double c = log(b)/log(a);

	return (uint64_t)ceil(c);
}

void CHHash::delete_tree(_chhash_tree* t) {
	if(t == nullptr) {
		if(tree != nullptr) {
			delete_tree(tree);
			tree = nullptr;
		}
		return;
	}

	if(t->yes != nullptr)
		delete_tree(t->yes);

	if(t->no != nullptr)
		delete_tree(t->no);

	delete t;
}

CHHash::~CHHash() {
	delete_tree();
	delete C;
}

std::string CHHash::to_c_code(_chhash_tree* branch, const std::string& preamble) {
	std::string output = preamble + "if (key[" + std::to_string(branch->i) + "] == ";
	if(branch->c != 0)
		output += std::string("'") + branch->c + "') {\n";
	else
		output += "'\\0') {\n";

	if(branch->yes == nullptr) {
		output += preamble + "\treturn " + std::to_string(branch->h) + ";\n";
	} else {
		output += to_c_code(branch->yes, preamble + "\t");
	}
	
	if(branch->no != nullptr) {
		output += preamble + "} else {\n";
		if(branch->no->yes == nullptr && branch->no->no == nullptr) {
			output += preamble + "\treturn " + std::to_string(branch->no->h) + ";\n";
		} else {
			output += to_c_code(branch->no, preamble + "\t");
		}
	}

	output += preamble + "}\n";

	return output;
}

void CHHash::generate_tree(const uint64_t& number_keys, std::string* keys) {
	uint64_t number_hashes = 0;
	tree = init_tree(number_keys, keys, 0, &number_hashes, nullptr);
}

void CHHash::generate_tree(const uint64_t& number_keys, std::string keys_file) {
	std::ifstream f;
	f.open(keys_file);

	std::vector<std::string> keys;
	keys.reserve(number_keys);

	uint64_t i = 0;	
	while(i < number_keys && !f.eof()) {
		char buffer[128];
		f.getline(buffer, 1024, '\n');
		std::string k = std::string(buffer);
		keys.push_back(k);
		i++;
	}

	std::string* first = &keys[0];

	generate_tree(number_keys, first);
}

std::string CHHash::to_c_code() {
	if(this->tree == nullptr)
		return "int hash(const char* key) {\n\treturn 0; //No tree was generated\n}\n";
	return "int hash(const char* key) {\n" + to_c_code(this->tree, "\t") + "}\n";
}

uint64_t CHHash::get(const std::string& key, _chhash_tree* branch) {
	if(branch == nullptr) return get(key, tree);

	bool check = (key.length() > branch->i && key[branch->i] == branch->c) || (key[branch->i] == 0 && key.length() == branch->i);

	if(check) {
		if(branch->yes == nullptr) {
			return branch->h;
		} else {
			return get(key, tree->yes);
		}
	} else {
		if(branch->no->yes == nullptr && branch->no->no == nullptr) {
			return branch->no->h;
		} else {
			return get(key, tree->no);
		}
	}
}