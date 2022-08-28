#pragma once

#include <random>
#include <stdint.h>

#include "make_array.hpp"
#include "q_ary.hpp"
#include "matrix.hpp"

//Returns a vector 'u' for which it's dot product with 'v' yields 0. Base 'q', size 'n'
//'u' is output, 'sum' how much is left to sum
template <typename T>
T* zero_dot_product(T* v, const uint64_t& n, const uint64_t& q, T* u = nullptr, const uint64_t& i = 0, uint64_t sum = 0) {
	bool nonRecursive = u == nullptr;
	//Non-recursive call to function
	if(nonRecursive) {
		u = make_array<T>(n, 1);
		sum = 0;
		for(uint64_t i = 0; i < n; i++) {
			sum += v[i];
		}
		sum = sum % q;
	}

	//std::cout << "sum: " << sum << " / u: "; _print_array(u, n); std::cout << "\n";

	//Resursive bit
	if(sum == 0)
		return u;
	if(i < n) {
		//Case 1: +1 value of current index and check sum
		if(u[i] < (q - 1)) {
			u[i]++;
			T* u_valid = zero_dot_product(v, n, q, u, i, (sum + v[i]) % q);
			if(u_valid != nullptr)
				return u;
			u[i]--;
		}
		
		//Case 2: Skip to next index
		T* u_valid = zero_dot_product(v, n, q, u, i + 1, sum);
		if(u_valid != nullptr)
			return u;
	}

	//Non-resursive end to function
	if(nonRecursive) {
		free(u);
	}

	return nullptr;
}

//Obtain basis of vectors in columns
template <typename T>
matrix<T> matrix_basis(q_ary* C) {
	auto M = C->C;
}

/*-------------------------------------------------
	Global Functions
-------------------------------------------------*/

//Sample an uniformly random matrix (Base q) with specific sizes
matrix<qint>* sample_uniform_matrix(const uint64_t& q, const uint64_t& height, const uint64_t& width, const uint64_t& seed);

//Combinatorial numbers
double combinatorialD(const uint64_t& a, const uint64_t& b);

//Combinatorial numbers
uint64_t combinatorial(const uint64_t& a, const uint64_t& b);

//Factorial numbers
uint64_t factorial(const uint64_t& i);