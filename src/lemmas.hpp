#pragma once

#include <stdint.h>
#include "q_ary.hpp"

uint64_t generate_qary_lemma_2_1_M_fast(const uint64_t& m, const uint64_t& q, const uint64_t& N, const uint64_t S_C2, const uint64_t MAX_M);

uint64_t generate_qary_lemma_2_1_M(const uint64_t& m, const uint64_t& q, const uint64_t& N, const uint64_t S_C2, const uint64_t MAX_M);

//Constructor of an q-ary according to Lemma 2.1, Remark 1
q_ary* generate_qary_lemma_2_1(const uint64_t& m, const uint64_t& q, const int64_t& N, const uint64_t& M_max, matrix<uint64_t>* A = nullptr, const uint64_t& seed = 0, const bool& check_for_correct_parameters = false);

//Construct a q-ary by concatination according to Lemma 2.3 - Leave A nullptr to calculate it
q_ary* generate_qary_lemma_2_3(uint64_t N, uint64_t M, q_ary* C2, matrix<uint64_t>* A = nullptr, const uint64_t& seed = 0);

//Lemma 3.1 - G is a finite set of q-tuples
void lemma_3_1(q_ary* C, matrix<qint>* G);

//Lemma 3.3
q_ary* lemma_3_3(const matrix<qint>& G, uint64_t q);

//Lemma 3.14 Case 1
q_ary* lemma_3_14_case_1(uint64_t r, qint alpha);

//Lemma 3.12
q_ary* lemma_3_12(uint64_t q);