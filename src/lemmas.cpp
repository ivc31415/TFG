#include <math.h>
#include "lemmas.hpp"
#include "util.hpp"

/*-------------------------------------------------
	Lemmas
-------------------------------------------------*/

uint64_t generate_qary_lemma_2_1_M_fast(const uint64_t& m, const uint64_t& q, const uint64_t& N, const uint64_t S_C2, const uint64_t MAX_M) {
	double left = ((double)(2*q))/factorial(q) * pow((1 - (factorial(q) * (double)S_C2)/(double)pow(m, q)), (double)N);
	uint64_t rangeMin = q;
	uint64_t rangeMax = MAX_M;
	uint64_t M = MAX_M;
	double right = 1;

	//First check if MAX is in range or not
	for(uint64_t i = 1; i < q; i++) {
		right *= (M - i);
	}
	right = 1/right;
	if(left <= right)
		return MAX_M;

	while(rangeMin != rangeMax) {
		//Calculate Midpoint
		M = rangeMin + (rangeMax - rangeMin) / 2;

		//std::cout << "Range: [" << rangeMin << ", " << rangeMax << "] (Midpoint: " << M << ")\n";

		//Calculate right side
		right = 1;
		for(uint64_t i = 1; i < q; i++) {
			right *= (M - i);
		}
		right = 1/right;

		//Check inequality
		if(left <= right) {
			//std::cout << "\t" << left << " <= " << right << "\n";

			//Calculate right side, Again, but for M + 1
			right *= (M - q + 1)/(double)M;

			if(left > right) {
				//std::cout << "Found M! (M = " << M << ")\n";
				return M;
			} else {
				rangeMin = M;
			}
		} else {
			//std::cout << "\t" << left << " > " << right << "\n";
			rangeMax = M;
		}
	}

	return 0;
}

uint64_t generate_qary_lemma_2_1_M(const uint64_t& m, const uint64_t& q, const uint64_t& N, const uint64_t S_C2, const uint64_t MAX_M) {
	uint64_t M = q;
	bool first = true;
	while(combinatorial(M, q) * pow((1 - (factorial(q) * (double)S_C2)/(double)pow(m, q)), (double)N) <= (double)M/(2*q) && M <= MAX_M) {
		M++;
		first = false;
	}
	if(first) return 0;
	return M - 1;
}

q_ary* generate_qary_lemma_2_1(const uint64_t& m, const uint64_t& q, const int64_t& N, const uint64_t& M_max, matrix<uint64_t>* A, const uint64_t& seed, const bool& check_for_correct_parameters) {
	/*-------------------------------------------------------------------------------------------------------------------
		Steps:
			1. Sample M codewords in [q]^N (M N-tuples of base q) with "replacement" (Replacement means elements can be selected more than once)
			2. Remove codewords that lie in any A-unfriendly q-sized subset
				2.1 Because A is going to be the collection of all q-sized subsets of [m], this step can be skipped
			3. Remove at most M/2 codewords
				3.1 We can simply sample M/2 codewords
			4. Remove the o(M) codewords with collisions
		Extra notes:
			1. |A|, because A is all the combinations of ([q] m), equals (q m)
	-------------------------------------------------------------------------------------------------------------------*/
	
	uint64_t M = M_max;

	if(check_for_correct_parameters) {
		//M = generate_qary_lemma_2_1_M_fast(m, q, N, A->width, M_max);

		//We check is the parameters are correct. The check is specified as figure 1 at Lemma 2.1
		if(M < q) throw std::invalid_argument("M has to be equal or greater than q");
		if(m < q) throw std::invalid_argument("m has to be equal or greater than q");
		double right_side = (double)M/(2*q);
		double left_side = (double)combinatorial(M, q) * pow((1 - (factorial(q) * (double)A->width)/(double)pow(m, q)), (double)N);
		if(left_side > right_side) {
			std::cout << "q: " << q << ", m: " << m << ", N: " << N << ", M: " << M << ", size(A): " << A->width << " (" << left_side << " > " << right_side << ")\n";
			std::cout << "part: " << pow((1 - (factorial(q) * (double)A->width)/(double)pow(m, q)), (double)N) << "\n";
			std::cout << "part: " << combinatorial(M, q) << "\n";
			throw std::invalid_argument("Called generate_qary_lemma_2_1 with invalid arguments. Check Lemma 2.1, Figure 1 for more details");
		}
	}

	//Obtain matrix for step 4
	q_ary* qary;
	if(A == nullptr) {
		// Uniform sampling of M codewords
		qary = new q_ary(sample_uniform_matrix(m, N, M, seed), m);
	} else {
		// Uniform sampling with replacement of M codewords
		qary = new q_ary(sample_uniform_matrix(m, N, M, seed), m);
		//Remove codewords in A-unfriendly subsets (Slow)
		qary->remove_unfriendly_codewords(q, A, M/2);
	}
	matrix<qint>* mtrx = qary->C;

	//o(M) in the paper	
	double oM = combinatorialD(M, 2) * pow((double)m, -(double)N);

	//Delete collisions
	uint64_t _possible_collision_counter = oM + 1;
	uint64_t k = qary->k; //Number of tuples in the final q-ary's matrix
	bool* codewords_to_keep = make_array<bool>(qary->k, true);

	//Detect repetitions
	for(uint64_t i1 = 1; i1 < qary->k; i1++) {
		if(_possible_collision_counter == 0)
			break;
		auto a1 = mtrx->data[i1];
		for(uint64_t i2 = 0; i2 < i1; i2++) {
			if(codewords_to_keep[i2]) {
				auto a2 = mtrx->data[i2];
				if(std::equal(a1, a1 + N, a2)) { //Are the arrays equal?
					codewords_to_keep[i1] = false;
					_possible_collision_counter--;
					k--;
					break;
				}
			}
		}
	}

	//Create the output q-ary
	q_ary* output = new q_ary(q, N, k);
	uint64_t j = 0;
	for(uint64_t i = 0; i < qary->k; i++) {
		if(codewords_to_keep[i]) {
			for(uint64_t l = 0; l < (uint64_t)N; l++) {
				output->C->data[j][l] = mtrx->data[i][l];
			}
			j++;
		}
	}

	delete qary;
	free(codewords_to_keep);

	return output;
}

//Construct C1 and a q-ary by concatination according to Lemma 2.3
q_ary* generate_qary_lemma_2_3(uint64_t N, uint64_t M, q_ary* C2, matrix<uint64_t>* A, const uint64_t& seed) {
	/*-------------------------------------------------------------------------------------------------------------------
		Steps:
			1. Obtain S(C2): Collection of all q-element subsets of C2 that are separated
			2. Make a biyection π from each tuple of C2 to [m]
			3. Define A as union of every subset in S(C2) with each tuple passing through the biyection
			4. Concatenate by Inverse-biyect every C1 tuple subset (Look at the article)
	-------------------------------------------------------------------------------------------------------------------*/

	//If we don't have the A matrix, calculate it
	bool calculateA = (A == nullptr);
	if(calculateA) {
		// Obtain S(C2)
		std::cout << "Calculating S(C2)...\n";
		bool* I = make_array<bool>(C2->k, false);
		std::vector<bool*> SC2;
		check_separated_aux(C2, C2->q, I, &SC2);

		// Make biyection
		// Codeword 0 in C2 goes to 0, codeword 1 to 1, and so on...

		// Define A

		std::cout << "Calculating A...\n";
		std::vector<uint64_t*> A_V; //This could have been initialized directly as a matrix
		for(bool* separated_subset : SC2) {
			uint64_t* A_V_s = make_array<uint64_t>(C2->q);
			uint64_t j = 0;
			for(uint64_t i = 0; i < C2->k; i++) {
				if(separated_subset[i]) {
					A_V_s[j] = i;
					j++;
				}
			}
			A_V.push_back(A_V_s);
		}

		A = new matrix<uint64_t>(C2->q, A_V.size());

		for(uint64_t i = 0; i < A_V.size(); i++) {
			A->data[i] = A_V[i];
		}
		
		free(I);
		for(bool* a : SC2) {
			free(a);
		}
	}
	
	//Generate C1
	q_ary* C1 = generate_qary_lemma_2_1(C2->k, C2->q, N, M, A, seed, true);
	
	//Concatenate
	q_ary* C = new q_ary(C2->q, C2->n * C1->n, C1->k);
	for(uint64_t c = 0; c < C->k; c++) { //For each codeword in C1
		qint* C1_cw = C1->codeword(c);
		for(uint64_t i = 0; i < C1->n; i++) { //For each element in that codeword
			qint* C2_cw = C2->codeword(C1_cw[i]);
			for(uint64_t j = 0; j < C2->n; j++) { //Copy array to specific part of C's codew.
				C->codeword(c)[j + i * C2->n] = C2_cw[j];
			}
		}
	}

	// Free Data
	delete C1;
	if(calculateA) delete A;
	return C;
}

//Fills output with a mask for M's columns which i-th row equals Gi
bool lemma_3_1_check_equal_row(const uint64_t q, const uint64_t k, const uint64_t i, 
								bool* output, matrix<qint>* M, qint* Gi, 
								const uint64_t depth = 0, const uint64_t n_index_selected = 0) {
 
	if(n_index_selected == q) {
		bool match = true;
		for(uint64_t j = 0; j < k; j++) {
			if(output[j] && M->data[i][j] != Gi[j]) {
				match = false;
				break;
			}
		}
		return match;
	} else if((k - depth) >= (q - n_index_selected)) {
		output[depth] = true;
		if(lemma_3_1_check_equal_row(q, k, i, output, M, Gi, depth + 1, n_index_selected + 1))
			return true;
		output[depth] = false;
		return lemma_3_1_check_equal_row(q, k, i, output, M, Gi, depth + 1, n_index_selected);
	}
	return false;
}

uint64_t* lemma_3_1_normal_A(const uint64_t n, const uint64_t k, matrix<bool>* FancyA,
						uint64_t* results = nullptr, bool* mask = nullptr,
						const uint64_t depth = 0, const uint64_t n_index_selected = 0) {

	bool generate_arrays = (results == nullptr);
	if(generate_arrays) {
		mask = make_array<bool>(n, false);
		results = make_array<uint64_t>(n, 0);
	}

	if(depth == n) {
		if(n_index_selected == 0)
			return nullptr;

		//Construct FancyA_T, and at the same time, calculate |FancyA_T| and add it to the results
		bool* AT_current = make_array<bool>(k, true);
		uint64_t counter = k;
		for(uint64_t j = 0; j < n; j++) {
			//For each FancyA_j that composes FancyA_T
			if(mask[j]) {
				for(uint64_t i = 0; i < k; i++) {
					if(AT_current[i] && !(FancyA->data[j][i])) {
						AT_current[i] = false;
						counter--;
					}
				}
			}
		}
		//Add to Ai (i = n_index_selected)
		results[n_index_selected] += counter;

	} else {
		mask[depth] = true;
		lemma_3_1_normal_A(n, k, FancyA, results, mask, depth + 1, n_index_selected + 1);
		mask[depth] = false;
		lemma_3_1_normal_A(n, k, FancyA, results, mask, depth + 1, n_index_selected);
	}

	if(generate_arrays) {
		free(mask);
		return results;
	}

	return nullptr;
}

void lemma_3_1(q_ary* C, matrix<qint>* G) {
	/*-------------------------------------------------------------------------------------------------------------------
	Notes:
		Distance in a linear code is minimun weight of its non-zero codewords
		Dual Distance is... ?????
	Steps:
		User inputed G
		User inputed C

		1- For each element i in G
			1.1- Check for q-sized subsets of C which row i equals element i of G
			1.2- You get Fancy Ai
		2- For every i in [n]
			2.1- Sum the size of every FancyA_T that has T of [n] of size i
	-------------------------------------------------------------------------------------------------------------------*/

	//if(C->q != G->height) {
	//	throw std::invalid_argument("Lemma 3.1 requires q (from the q-ary) and G's height (size of each element) to be the same");
	//}

	//if(C->n != G->width) {
	//	throw std::invalid_argument("Lemma 3.1 requires n (from the q-ary) and G's width (number of elements) to be the same");
	//}

	uint64_t q = 4; //= C->q;
	uint64_t k = 6; //= C->k;
	uint64_t n = 6; //= C->n;

	matrix<bool>* FancyA = new matrix<bool>(k, n, false);
	for(uint64_t i = 0; i < n; i++) {
		bool* FancyAi = FancyA->data[i];
		qint* Qi = G->data[i];
		lemma_3_1_check_equal_row(q, k, i, FancyAi, nullptr, nullptr);
	}
	
	uint64_t* Ai = lemma_3_1_normal_A(n, k, FancyA);
}

q_ary* lemma_3_3(const matrix<qint>& G, uint64_t q) {
	//The order 'q' of the abelian group 'G' is how many symbols it uses
	q_ary* C = new q_ary(q, G.height + 1, G.width);
	for(uint64_t j = 0; j < G.width; j++) {
		qint last_element = 0;
		for(uint64_t i = 0; i < G.height; i++) {
			qint to_copy = G.data[j][i];
			C->C->data[j][i] = G.data[j][i] % q;
			last_element += to_copy;
		}
		C->C->data[j][G.height] = last_element % q;
	}

	//Generate A_i for i = [0, n-1]
	uint64_t n = G.height + 1;
	std::vector<uint64_t> A;
	for(uint64_t i = 0; i < G.height; i++) {
		A.push_back(pow(q, q * (n - i)) * pow(factorial(q), i));
	}

	return nullptr;
}

//Lemma 3.14 Case 1
q_ary* lemma_3_14_case_1(uint64_t r, qint alpha) {
	uint64_t q = 1;
	for(uint64_t i = 0; i < r; i++)
		q *= 2;

	//Number of possible codewords
	uint64_t k = q * q;

	q_ary* C = new q_ary(q, 4, k);

	uint64_t counter = 0;
	for(uint64_t x = 0; x < q; x++) {
		for(uint64_t y = 0; y < q; y++) {
			C->C->data[counter][0] = x;
			C->C->data[counter][1] = y;
			C->C->data[counter][2] = (x + y) % q;
			C->C->data[counter][3] = (x + alpha * y) % q;

			counter++;
		}
	}

	return C;
}

//Lemma 3.12
q_ary* lemma_3_12(uint64_t q) {
	if(q%2 == 0 || q < 3) {
		throw std::invalid_argument("When using Lemma 3.12, q must be odd and greater than 3");
	}

	q_ary* C = new q_ary(q, 4, q * q);
	uint64_t counter = 0;
	for(uint64_t x = 0; x < q; x++) {
		for(uint64_t y = 0; y < q; y++) {
			C->C->data[counter][0] = x;
			C->C->data[counter][1] = y;
			C->C->data[counter][2] = (x + y) % q;
			C->C->data[counter][3] = (x - y + q) % q;

			counter++;
		}
	}

	return C;
}

q_ary* lemma_3_14_case_2(uint64_t r, qint alpha) {
	uint64_t q = 1;
	for(uint64_t i = 0; i < r; i++)
		q *= 2;

	//Number of possible codewords
	uint64_t k = q * q;

	q_ary* C = new q_ary(q, 4, k);

	uint64_t counter = 0;
	for(uint64_t x = 0; x < q; x++) {
		for(uint64_t y = 0; y < q; y++) {
			C->C->data[counter][0] = x;
			C->C->data[counter][1] = y;
			C->C->data[counter][2] = (x + y) % q;
			C->C->data[counter][3] = (x + alpha * y) % q;

			counter++;
		}
	}

	return C;
}