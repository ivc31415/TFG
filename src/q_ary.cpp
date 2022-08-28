#include <iostream>
#include <algorithm>
#include <math.h>
#include <vector>

#include "q_ary.hpp"
#include "make_array.hpp"

/*-------------------------------------------------
	Q-Ary
-------------------------------------------------*/

q_ary::q_ary(matrix<qint>* C, uint64_t q) {
	this->C = C;
	this->q = q;
	this->n = C->height;
	this->k = C->width;
}

q_ary::q_ary(uint64_t q, uint64_t n, uint64_t k) {
	this->q = q;
	this->n = n;
	this->k = k;
	C = new matrix<qint>(n, k);
}

q_ary::~q_ary() {
	delete C;
	C = nullptr;
}

qint* q_ary::codeword(const uint64_t& i) {
	return C->data[i];
}

//Tree search to count separated q-subsets = |S(C)|
void count_separated_subsets(q_ary* C, uint64_t* output, bool* I, const uint64_t step, const uint64_t selected) {
	bool first = I == nullptr;
	if(first) {
		I = make_array<bool>(C->k, false);
		*output = 0;
	}

	if(selected == C->q) {
		//Check if separated
		q_ary_subset s = q_ary_subset(C, I);
		if(s.separated())
			(*output)++;
	} else if(step < C->k) {
		I[step] = true;
		count_separated_subsets(C, output, I, step + 1, selected + 1);
		I[step] = false;
		count_separated_subsets(C, output, I, step + 1, selected);
	}

	if(first) {
		free(I);
	}
}

//Tree search for every m-element subset
bool check_separated_aux(q_ary* ary, const uint64_t& m, bool* I, std::vector<bool*>* separated_subsets, std::vector<bool*>* unseparated_subsets, const uint64_t depth, const uint64_t n_index_selected) {
	bool output = true;

	if(n_index_selected == m) {
		q_ary_subset s = q_ary_subset(ary, I);
		output = s.separated();
		if(separated_subsets != nullptr && output) { //Push to friendly subsets list if user wants
			//std::cout << "Found separated subset\n";
			separated_subsets->push_back(copy_array(ary->k, I)); //Executable corrupted?
		}
		else if(unseparated_subsets != nullptr && !output) { //Push to unfriendly subsets list if user wants
			//std::cout << "Found unseparated subset\n";
			unseparated_subsets->push_back(copy_array(ary->k, I)); //Executable corrupted?
		}
	} else if((ary->k - depth) + n_index_selected >= m) {
		I[depth] = true;
		output = check_separated_aux(ary, m, I, separated_subsets, unseparated_subsets, depth + 1, n_index_selected + 1);

		I[depth] = false;
		output = check_separated_aux(ary, m, I, separated_subsets, unseparated_subsets, depth + 1, n_index_selected) || output;
	}

	/*std::cout << "I: (";
	for(uint64_t i = 0; i < ary->k; i++) {
		if(i < depth)
			std::cout << I[i] << ",";
		else
			std::cout << "?,";
	}
	std::cout << "\b), " << n_index_selected << "/" << depth << ", separated: " << output << "\n";*/

	return output;
}

bool q_ary::check_mhash(const uint64_t& m) {
	bool* I = make_array<bool>(k, false);
	bool output = check_separated_aux(this, m, I);
	free(I);

	return output;
}

//Tree search for every q-element subset
bool check_friendly_aux(q_ary* ary, const uint64_t& q, matrix<uint64_t>* A, bool* I, std::vector<bool*>* friendly_subsets = nullptr, std::vector<bool*>* unfriendly_subsets = nullptr, const uint64_t depth = 0, const uint64_t n_index_selected = 0) {
	bool output = true;
	uint64_t m = A->width;

	if(n_index_selected == q) {
		q_ary_subset s = q_ary_subset(ary, I);
		//std::cout << "I: "; _print_array(I, ary->k); std::cout << " (Depth: " << depth << ")\n";
		output = s.friendly(A); //New, faster way to check if subgroup is A friendly using divide and conquer. For the old one, see q_ary_subset::friendly_old()

		if(friendly_subsets != nullptr && output) { //Push to friendly subsets list if user wants
			//std::cout << "Found friendly subset\n";
			friendly_subsets->push_back(copy_array(ary->k, I)); //Executable corrupted?
		}
		else if(unfriendly_subsets != nullptr && !output) { //Push to unfriendly subsets list if user wants
			//std::cout << "Found unfriendly subset\n";
			unfriendly_subsets->push_back(copy_array(ary->k, I)); //Executable corrupted?
		}
	} else if((ary->k - depth) + n_index_selected >= q) {
		//std::cout << "ary->k - depth + n_index_selected >= q  ---> " << ary->k << " - " << depth << ") + " << n_index_selected << " >= " << q << "\n";

		I[depth] = true;
		output = check_friendly_aux(ary, q, A, I, friendly_subsets, unfriendly_subsets, depth + 1, n_index_selected + 1);
		I[depth] = false;
		if(output || n_index_selected == 0) { //A subset size m was found that was not separated
			output = check_friendly_aux(ary, q, A, I, friendly_subsets, unfriendly_subsets, depth + 1, n_index_selected);
		}
	}

	return output;
}

bool q_ary::check_friendly(matrix<uint64_t>* A, const uint64_t& q_A) {
	bool* I = make_array<bool>(k, false);
	bool output = check_friendly_aux(this, q_A, A, I);
	free(I);

	return output;
}

std::vector<bool*> q_ary::get_unfriendly_subsets(const uint64_t& q_A, matrix<uint64_t>* A) {
	std::vector<bool*> unfriendly_subsets_mask;
	bool* I = make_array<bool>(k, false);
	//Run binary search tree
	check_friendly_aux(this, q_A, A, I, nullptr, &unfriendly_subsets_mask);
	free(I);

	return unfriendly_subsets_mask;
}

void q_ary::remove_codewords(bool* mask) {
	bool* to_keep_mask = make_array<bool>(k, false);
	int new_k = 0;

	for(uint64_t i = 0; i < k; i++) {
		if(!mask[i]) {
			to_keep_mask[i] = true;
			new_k++;
		}
	}

	//Get new matrix
	matrix<qint>* C2 = new matrix(this->C, to_keep_mask);

	//Replace matrix
	delete C;
	C = C2;
	k = new_k;

	//Free data
	free(to_keep_mask);
}

void q_ary::remove_unfriendly_codewords(const uint64_t& m, matrix<uint64_t>* A, const uint64_t& max) {
	//Get unfriendly subsets
	std::vector<bool*> unfriendly_subsets = get_unfriendly_subsets(m, A);

	//Get a mask of codewords to keep
	uint64_t new_k = k;
	uint64_t counter = 0;
	bool* to_keep_mask = make_array<bool>(k, true);
	for(bool* single_mask : unfriendly_subsets) {
		for(uint64_t i = 0; i < k; i++) {
			if(single_mask[i] && to_keep_mask[i]) {
				to_keep_mask[i] = false;
				new_k--;
			}
		}
		free(single_mask);
		counter++;
		if(max != 0 && counter == max) {
			break;
		}
	}

	//Get new matrix
	matrix<qint>* C2 = new matrix(this->C, to_keep_mask);
	free(to_keep_mask);

	//Replace matrix
	delete C;
	C = C2;
	k = new_k;
}

uint64_t q_ary::codeword_weight(const uint64_t& i) {
	qint* cw = codeword(i);
	uint64_t o = 0;

	for(uint64_t _n = 0; _n < n; _n++) {
		o += (cw[_n] != 0);
	}

	return o;
}

uint64_t q_ary::distance() {
	uint64_t d = UINT64_MAX;
	
	for(uint64_t i = 1; i < k; i++) {
		for(uint64_t j = 0; j < i; j++) {
			uint64_t d_candidate = hamming_distance(i, j);
			d = d_candidate < d ? d_candidate : d;
		}
	}

	return d;
}

uint64_t q_ary::hamming_distance(const uint64_t& i, const uint64_t& j) {
	qint* cw_i = codeword(i);
	qint* cw_j = codeword(j);
	uint64_t o = 0;

	for(uint64_t _n = 0; _n < n; _n++) {
		o += (cw_i[_n] != cw_j[_n]);
	}

	return o;
}

//Get the rate
double q_ary::rate() {
	return log2((double)k)/((double)n);
}

q_ary_subset q_ary::create_subset(bool* I) {
	return q_ary_subset(this, I);
}

void q_ary::write_codeword(const uint64_t codeword_number, const char* data) {
	qint* cw = C->data[codeword_number];
	for(uint64_t i = 0; i < n; i++)
		cw[i] = ((qint)data[i])%q;
}

void q_ary::write_codeword_raw(const uint64_t codeword_number, const qint* data) {
	qint* cw = C->data[codeword_number];
	for(uint64_t i = 0; i < n; i++)
		cw[i] = ((qint)data[i])%q;
}

q_ary* q_ary::get_dual_code() {
	q_ary* output = new q_ary(q, n, k);

	//https://arxiv.org/pdf/cs/0506087.pdf
	//	The dual code C⊥ of a linear code C is defined as C⊥ = {u | u · v = 0 for all v ∈ C}
	//		Where u · v is the sum of the product of every i member of the codeword
	//	We have to perform Gauss Elimination

	for(uint64_t ik = 0; ik < k; ik++) {
		qint* cw_orig = codeword(ik);
		qint* cw_dual = output->codeword(ik);

		uint64_t t = 0;
		for(uint64_t in = 0; in < n; in++) {
			t += cw_orig[in];
		}
		uint64_t delta = q - (t % q);
		
	}

	return output;
}

uint64_t q_ary::dual_distance() {
	q_ary* dual_code = get_dual_code();
	uint64_t output = dual_code->distance();
	delete dual_code;
	return output;
}

/*-------------------------------------------------
	Q-Ary Subset
-------------------------------------------------*/

bool q_ary_subset::separated() {
	for(uint64_t i = 0; i < q->n; i++) {
		bool* seen = make_array<bool>(q->q, false);
		bool check = true;
		for(uint64_t _k = 0; _k < q->k; _k++) {
			if(I[_k]) { //Check only tuples on the subset
				uint64_t value = q->C->data[_k][i];
				if(seen[value]) {
					check = false;
					break;
				}
				seen[value] = true;
			}
		}
		free(seen);
		if(check)
			return true;
	}
	return false;
}

bool q_ary_subset::friendly(matrix<uint64_t>* A) {
	uint64_t* J = make_array<uint64_t>(n);
	uint64_t* row = make_array<uint64_t>(n);
	bool* seen = make_array<bool>(q->q);
	uint64_t start, end, mid, stage;

	//Transform the mask that defines the subgroup into a list with the position of the subgroup codewords in the linear code
	uint64_t l = 0;
	for(uint64_t i = 0; i < q->k; i++) {
		if(this->I[i]) {
			J[l++] = i;
		}
	}

	//For every index of each codeword
	for(uint64_t i = 0; i < q->n; i++) {
		//Get Row as Array
		bool abort = false;
		for(uint64_t j = 0; j < q->q; j++) seen[j] = false;
		for(uint64_t j = 0; j < n; j++) {
			uint64_t c = q->codeword(J[j])[i];
			row[j] = c;
			if(seen[c]) {
				abort = true;
				break;
			}
			seen[c] = true;
		}
		if(abort) { //Row has repeated numbers, SKIP
			continue;
		}

		//Divide and Conquer
		sort_array(row, n);
		start = 0;
		end = A->width - 1;
		stage = 0;

		while(stage < n) {
			mid = (start + end)/2;
			
			//Check how far both tuples are equal. Stage indicates the first element between the two that are not equal.
			stage = 0;
			while(stage < n && A->data[mid][stage] == row[stage]) {
				stage++;
			}

			//They are equal
			if(stage >= n) {
				free(J);
				free(seen);
				free(row);
				return true;
			}

			//Divide and Conquer ended, they are not equal
			if(start == end) {
				break;
			}

			//Divide and conquer's subdivision of search space
			if(A->data[mid][stage] > row[stage]) {
				end = mid;
			} else {
				start = mid + 1;
			}
		}
	}

	free(J);
	free(seen);
	free(row);

	return false;
}

bool q_ary_subset::friendly_old(matrix<uint64_t>* A) {
	/*
		Warning. This function will yield false negatives if the matrix A doesn't have all the permutations for every subset
		It's better to use the function q_ary_subset::friendly() instead, which doesn't have this problem
	*/
	for(uint64_t k = 0; k < A->width; k++) { //For each set in A
		auto a_set = A->data[k];
		for(uint64_t i = 0; i < q->n; i++) { //For each index i
			bool tuples_equal = true;
			uint64_t l = 0;
			for(uint64_t j = 0; j < q->k; j++) { //Compare with the set of the elements with index i
				if(this->I[j]) {
					auto q_set = this->q->codeword(j);

					if(a_set[l++] != q_set[i]) {
						tuples_equal = false;
						break;
					}
				}
			}
			if(tuples_equal) {
				return true;
			}
		}
	}
	return false;
}

/*-------------------------------------------------
	Misc. Functions
-------------------------------------------------*/

void print_qary_or_subset(q_ary* q, bool* I) {
	std::cout << "(k: " << q->k << ", n: " << q->n << ", q: " << q->q << ")\n";
	for(uint64_t _k = 0; _k < q->k; _k++) {
		if(I == nullptr || I[_k]) {
			std::cout << "(";
			for(uint64_t _n = 0; _n < q->n; _n++) {
				if(_n != 0) std::cout << ", ";
				std::cout << q->C->data[_k][_n];
			}
			std::cout << ")\n";
		}
	}
}

void print_qary(q_ary* q) {
	print_qary_or_subset(q, nullptr);
}

void print_subset(q_ary_subset* s) {
	print_qary_or_subset(s->q, s->I);
}
