#ifndef Q_ARY
#define Q_ARY

//Type of data saved in a q-ary's matrix
#define qint uint16_t

#include <stdint.h>
#include <matrix.hpp>

struct q_ary_subset;

/*
	Data type that holds information on a q-ary
*/
struct q_ary {
	//Hold Data of the q_ary, specifially array of n-tuples (codeword)
	matrix<qint>* C;
	//q-ary's q (Base of every element in each tuple)
	uint64_t q;
	//How big each tuple is
	uint64_t n;
	//How many n-tuples are in the q-ary
	uint64_t k;

	//Constructors
	//From matrix
	q_ary(matrix<qint>* C, uint64_t q);

	//From matrix in file
	q_ary(std::string file_path, uint64_t q) : q_ary(new matrix<qint>(file_path), q) {}

	//Empty constructor
	q_ary(uint64_t q, uint64_t n, uint64_t k);

	//Destructor (Free the memory)
	~q_ary();

	//Retrieve the i-th codeword (n-tuple)
	qint* codeword(const uint64_t& i);

	//Is the q-ary also a m-hash?
	bool check_mhash(const uint64_t& m);

	//Is the q-ary A-friendly? A must be a (m x a) matrix, a <= (q m)
	bool check_friendly(matrix<uint64_t>* A, const uint64_t& q_A);

	//Get A-unfriendly q-sized sets
	std::vector<bool*> get_unfriendly_subsets(const uint64_t& q_A, matrix<uint64_t>* A);

	//Remove codewords based on a mask. If true, remove
	void remove_codewords(bool* mask);

	//Remove codewords that lie in any A-unfriendly m-sized subset
	void remove_unfriendly_codewords(const uint64_t& m, matrix<uint64_t>* A, const uint64_t& max = 0);

	//Retrieve the weight of a codeword
	uint64_t codeword_weight(const uint64_t& i);

	//Calculate distance of q-ary (Assuming it's a linear code)
	uint64_t distance();

	//Hamming distance between two codewords
	uint64_t hamming_distance(const uint64_t& i, const uint64_t& j);

	//Create a subset from an array of bools containing what tuples are included
	q_ary_subset create_subset(bool* I);

	//Write a codeword as text
	void write_codeword(const uint64_t codeword_number, const char* data);

	//Write a codeword
	void write_codeword_raw(const uint64_t codeword_number, const qint* data);

	//Calculate the dual code of this q-ary, which is C⊥ = {u | u · v = 0, v ∈ C}
	q_ary* get_dual_code();

	//Calculate the distance of the dual code. Has to generate the dual code
	uint64_t dual_distance();

	//Get the rate
	double rate();
};

//Tree search to count separated subsets = |S(C)|
void count_separated_subsets(q_ary* C, uint64_t* output, bool* I = nullptr, const uint64_t step = 0, const uint64_t selected = 0);

//Tree search for every m-element subset
bool check_separated_aux(q_ary* ary, const uint64_t& m, bool* I, std::vector<bool*>* separated_subsets = nullptr, std::vector<bool*>* unseparated_subsets = nullptr, const uint64_t depth = 0, const uint64_t n_index_selected = 0);

/*
	Holds a subset of a q-ary's data
*/
struct q_ary_subset {
	bool* I; //Tuples of the original subset that are in this one. Size of the original q-ary's k value
	q_ary* q; //Original q-ary
	uint64_t n; //Number of elements in subset

	q_ary_subset(q_ary* q, bool* I) {
		this->I = I;
		this->q = q;
		this->n = 0;
		for(uint64_t i = 0; i < q->k; i++) {
			if(I[i])
				n++;
		}
	}
	
	//Are the tuples in the subset separated?
	bool separated();

	//Is the subset friendly to A? Numer of elements in subset must equal height of A matrix
	bool friendly(matrix<uint64_t>* A);

	//OLD METHOD. Is the subset friendly to A? Numer of elements in subset must equal height of A matrix
	bool friendly_old(matrix<uint64_t>* A);
};

//Print a q-ary to the screen
void print_qary(q_ary* q);
//Print a q-ary subset to the screen
void print_subset(q_ary* s);

//Sample an uniformly random matrix (Base q) with specific sizes
matrix<qint>* sample_uniform_matrix(const uint64_t& q, const uint64_t& height, const uint64_t& width);

#endif