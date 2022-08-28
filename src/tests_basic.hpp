#include <iostream>

#include "q_ary.hpp"
#include "util.hpp"
#include "lemmas.hpp"
#include "hashmap.hpp"
#include "chrono2.hpp"
#include "chhash.hpp"

using namespace std;

void test_new_friendliness() {
	matrix<uint64_t> A("data/A_q_4.dat");

	CHHash c(20, 4, &A, 2);
}

void testArraySort() {
	uint64_t n = 4;
	uint64_t* v = make_array<uint64_t>(n);
	v[0] = 6; v[1] = 2; v[2] = 1; v[3] = 10;

	_print_array(v, n);
	std::cout << "\n";

	sort_array(v, n);

	_print_array(v, n);
	std::cout << "\n";
}

void testReadKeys() {
	CHHash ch(10, 4, "data");
	ch.generate_tree(1000, "../kepler_stellar_17_names_1000.txt");
	std::cout << ch.to_c_code() << "\n";
}

void testFastMCalc() {
	generate_qary_lemma_2_1_M_fast(25, 5, 300, 11006, 300000);
}

//int64_t testFaster_2_1(bool useOldOne);

int __gen_hash(const char* key) {
	if (key[0] == 'h') {
		if (key[1] == 'e') {
			return 0;
		} else {
			if (key[2] == '\0') {
				return 2;
			} else {
				return 1;
			}
		}
	} else {
		if (key[0] == 'k') {
			return 3;
		} else {
			return 4;
		}
	}
}

void testFinal() {
	std::string keys[] = {"hello", "hi", "kepler 832", "hip 4582", "wasp 21"};
	CHHash h(55, 40, 5);
	h.generate_tree(5, keys);

	cout << "-----------------------\n";
	cout << h.to_c_code();
}

void testOne(int type, uint64_t value) {
	uint64_t N = 100, M_max = 8;

	std::cout << "Calculating...\n";

	uint64_t q;
	if(type == 0)
		q = pow(2, value);
	else
		q = value;

	cout << "q: " << q << ", N: " << N << "\n";

	simple_clock c2, cc;
	
	c2.reset();
	q_ary* C2;
	if(type == 0)
		C2 = lemma_3_14_case_1(value, 3); //Inner Code
	else
		C2 = lemma_3_12(value);
	auto time2 = c2.microseconds();
	
	cout << "Times:\n";
	cout << "\tC2 Generation: " << (time2/1000.0) << "ms\n";

	q_ary* C;
	cc.reset();
	C = generate_qary_lemma_2_3(N, M_max, C2);
	auto timeC = cc.microseconds();

	cout << "\tConcatenation: " << (timeC/1000.0) << "ms\n";
	cout << "Saving Codes C & C2 (q: " << q << ")...\n";
	//print_qary(C);

	C2->C->save("C2.code");
	C->C->save("C.code");

	delete C2;
}

void testOne(int argc, char *argv[]) {
	int type = 0; //0: 3.14, 1: 3.12
	uint64_t value = 0;

	if(argc != 3)
		return;
	
	if(argv[1][0] == '1')
		type = 1;
	
	value = (uint64_t)(argv[2][0] - '0');
	
	testOne(type, value);
}

void calculate_sc_3_14(uint64_t r) {
	uint64_t N = 10, M = 20;

	uint64_t q = pow(2, r);
	M = (q*q) > M ? (q*q) : M;
	cout << "q: " << pow(2, r) << ", N: " << N << ", M: " << M << "\n";

	simple_clock c2, cSC;
	
	c2.reset();
	q_ary* C2 = lemma_3_14_case_1(r, 3); //Inner Code
	auto time2 = c2.microseconds();

	cSC.reset();
	uint64_t sizeSC;
	count_separated_subsets(C2, &sizeSC);
	auto timeSC = cSC.microseconds();

	cout << "Times:\n";
	cout << "\tC2 Generation: " << (time2/1000.0) << "ms\n";
	cout << "\tS(C2) search: " << (timeSC/1000.0) << "ms\n";
	cout << "C2: ";
	print_qary(C2);
	cout << "|S(C)|: " << sizeSC << "\n";

	delete C2;
}
/*
void testOneNew_2_1(int argc, char *argv[]) {
	uint64_t N = 300, M_max = 8;
	uint64_t r = 2; //q = 2^r
	bool useNew = false;

	if(argv != nullptr && argc > 2) {
		if(argv[2][0] == '3')
			r = 3;
		if(argv[1][0] == 'n')
			useNew = true;
	}

	std::cout << "Calculating...\n";

	uint64_t sizeSC;
	if(r == 2)
		sizeSC = 856;
	if(r == 3)
		sizeSC = 66762112; //α = 3
	if(r > 3)
		return;

	uint64_t q = pow(2, r);
	cout << "q: " << pow(2, r) << ", N: " << N << "\n";

	simple_clock c2, cc;
	
	c2.reset();
	q_ary* C2 = lemma_3_14_case_1(r, 3); //Inner Code
	auto time2 = c2.microseconds();
	
	cout << "Times:\n";
	cout << "\tC2 Generation: " << (time2/1000.0) << "ms\n";

	q_ary* C;
	cc.reset();
	if(useNew)
		C = generate_qary_lemma_2_3_v2(N, M_max, C2, sizeSC);
	else
		C = generate_qary_lemma_2_3(N, M_max, C2);
	auto timeC = cc.microseconds();

	cout << "\tConcatenation: " << (timeC/1000.0) << "ms\n";
	cout << "C2: ";
	//print_qary(C2);
	cout << "|S(C)|: " << sizeSC << "\n";
	cout << "C: ";
	//print_qary(C);

	delete C2;
}*/
/*
void testFaster_2_1_loop() {
	double timeOld = 0, timeNew = 0;
	const int64_t I = 1;

	for(int i = 0; i < I; i++) {
		timeOld += testFaster_2_1(true);
		timeNew += testFaster_2_1(false);
	}

	timeOld /= (I * 1000.0);
	timeNew /= (I * 1000.0);

	cout << "Average time for old algorithm: " << timeOld << "\n";
	cout << "Average time for new algorithm: " << timeNew << "\n";
	cout << "Old to New Ratio: x" << (timeOld/timeNew) << " faster\n";
}*/
/*
int64_t testFaster_2_1(bool useOldOne) {
	const uint64_t N = 35, M_max = 20;
	simple_clock c;
	c.reset();

	q_ary* C2 = lemma_3_14_case_1(2, 3); //Inner Code
	
	//q_ary* C1 = generate_qary_lemma_2_1_v2(C2->k, C2->q, N, M, C2, true);

	q_ary* C;
	if(!useOldOne)
		C = generate_qary_lemma_2_3_v2(N, M_max, C2, 856);
	else
		C = generate_qary_lemma_2_3(N, M_max, C2);

	delete C2;

	auto time = c.microseconds();

	if(C->k < 1)
		cout << "?? (useOldOne: " << useOldOne << ")\n";

	delete C;

	return time;
}*/

void testFirstCodes() {
	const uint64_t Q[] = {4, 5, 7, 8};
	for(uint64_t i = 0; i < 3; i++) {
		uint64_t q = Q[i];

		q_ary* C2;//Inner Code
		if(q == 4)
			C2 = lemma_3_14_case_1(2, 3); 
		else if(q == 8)
			C2 = lemma_3_14_case_1(3, 3); 
		else
			C2 = lemma_3_12(q);

		print_qary(C2); std::cout << "\n";
	}
}

void testConcatenation() {
	q_ary* C2 = lemma_3_14_case_1(2, 3); //Inner Code
	//print_qary(C2);
	q_ary* C = generate_qary_lemma_2_3(10, 20, C2);
}

void testUtilZeroDot() {
	qint v[] = {1, 5, 3, 2};
	auto* u = zero_dot_product(v, 4, 8);
	_print_array(u, 4); cout << "\n";
}

void testHashMap() {
	LCHashMap<int> map(10);
	map.add("hey", 34345);
	printf("Value at 'hey': %i\n", *map.get("hey"));
}

template class LCHashMap<int>;

void weirdTest() {
	std::cout << "Hello!\n";
	
	//Make 256-ary of 7 5-tuples (2 codewords)
	q_ary* ary_file = new q_ary("data/m1.dat", 256);
	cout << ary_file->C << "\n";

	//Make 8-ary of 10 4-tuples

	//auto ary2 = generate_qary_lemma_2_1(8, 8, 4, 8, nullptr, true);
	q_ary* ary2 = new q_ary(8, 8, 8);
	for(uint64_t i = 0; i < 8; i++) {
		for(uint64_t j = 0; j < 8; j++) {
			ary2->C->data[i][j] = i;
		}
	}
	//auto a = generate_qary_lemma_2_3(nullptr, ary2);

	//lemma_3_1(nullptr, nullptr);

	/*std::cout << ary->C << "\n";
	bool* mask = make_array<bool>(7, true); mask[2] = false;
	matrix<qint> m2(ary->C, mask);
	std::cout << m2 << "\n";*/
}

