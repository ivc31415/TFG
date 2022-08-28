#include "util.hpp"

//Sample an uniformly random matrix (Base q) with specific sizes
matrix<qint>* sample_uniform_matrix(const uint64_t& q, const uint64_t& height, const uint64_t& width, const uint64_t& seed) {
	//Setup matrix
	matrix<qint>* m = new matrix<qint>(height, width);

	//Setting up the Uniform Random Number Generator
	if(seed == 0) {
		std::random_device                  rand_dev;
		std::mt19937                        generator(rand_dev());
		std::uniform_int_distribution<qint> distr(0, q - 1);

		//Fill matrix with data
		for(uint64_t w = 0; w < width; w++) {
			for(uint64_t h = 0; h < height; h++) {
				m->data[w][h] = distr(generator);
			}
		}
	} else {
		#pragma GCC diagnostic push
		#pragma GCC diagnostic ignored "-Wnarrowing"
		std::mt19937                        generator{seed};
		#pragma GCC diagnostic pop
		std::uniform_int_distribution<qint> distr(0, q - 1);

		//Fill matrix with data
		for(uint64_t w = 0; w < width; w++) {
			for(uint64_t h = 0; h < height; h++) {
				m->data[w][h] = distr(generator);
			}
		}
	}

	return m;
}

//Combinatorial numbers
double combinatorialD(const uint64_t& a, const uint64_t& b) {
	if(a < b)
		throw std::invalid_argument("Combinatorial Number requires the top value (a) to be equal or greater than the lower value (b)");
	else if(a == b || b == 0)
		return 1;

	uint64_t b2 = (a < (2 * b)) ? (a - b) : b;
	uint64_t c = a - b2;

	double output = 1;
	for(uint64_t i = 2; i <= a; i++) {
		output *= i;
		if(i <= b2)
			output /= i;
		if(i <= c)
			output /= i;
	}

	return output;
}
uint64_t combinatorial(const uint64_t& a, const uint64_t& b) {
	return (uint64_t)combinatorialD(a,b);
}

//Factorial numbers
uint64_t factorial(const uint64_t& i) {
	uint64_t output = 1;
	for(uint64_t j = 2; j <= i; j++) {
		output *= j;
	}
	return output;
}