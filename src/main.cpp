#include "tests_hard.hpp"

using namespace std;

int main(int argc, char *argv[]) {

	simple_clock clk;
	clk.reset();

	//CHHash genera una familia de hashes perfectos con por lo menos x codewords y un valor específico de q (Tiene que ser q = 4, 5 o 7)
	uint64_t x = 10;
	uint64_t q = 4;
	CHHash c(x, q);

	//Pruebas
	//testFirstCodes();
	//generate_A_files();
	//test_run_many_2_auto_N_r(false, 20, 200, 10, 4);
	//test_run_many_2_auto_N(true, 10, 20, 1, 1);
	//test_size_N(1000, 200000, 1000);
	//test_size_M(10, 10000, 10);
	//test_run_many_2(100, 100, 1);

	//test_new_friendliness();
	//testArraySort();

	//testReadKeys();

	uint64_t t = clk.microseconds();
	std::cout << "Time of Test: " << (t/1000) << "ms\n";

	return 0;
}