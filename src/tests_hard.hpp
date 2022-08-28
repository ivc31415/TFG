#pragma once

#include <stdio.h>

#include "tests_basic.hpp"
#include "chhash.hpp"

//Generate cached A files
void generate_A_files() {
	const uint64_t Q[] = {4, 5, 7, 8};

	simple_clock c_code, c_a;
	uint64_t t_code, t_a;

	for(uint64_t i = 0; i < 3; i++) {
		/*
			Generate Data
		*/
		uint64_t q = Q[i];

		c_code.reset();
		q_ary* C2; //Inner Code
		if(q == 4)
			C2 = lemma_3_14_case_1(2, 3);
		else if(q == 8)
			C2 = lemma_3_14_case_1(3, 3);
		else if(q == 5 || q == 7)
			C2 = lemma_3_12(q);
		t_code = c_code.microseconds();

		c_a.reset();
		// Obtain S(C2)
		std::cout << "Calculating S(C2)...\n";
		bool* I = make_array<bool>(C2->k, false);
		std::vector<bool*> SC2;
		check_separated_aux(C2, C2->q, I, &SC2);

		//std::cout << "S(C2) (Size: " << SC2.size() << ")\n";
		/*for(bool* a : SC2) {
			std::cout << "\t";
			_print_array(a, C2->k);
			std::cout << "\n";
		}*/

		// Make biyection
		// Codeword 0 in C2 goes to 0, codeword 1 to 1, and so on...

		// Define A

		std::cout << "Calculating A...\n";
		std::vector<uint64_t*> A_V;
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
		matrix<uint64_t> A(C2->q, A_V.size());
		for(uint64_t i = 0; i < A_V.size(); i++) {
			A.data[i] = A_V[i];
		}
		t_a = c_a.microseconds();

		/*
			Save Data
		*/
		A.save("data/A_q_" + std::to_string(q) + ".dat");
		std::cout << "q: " << q << " - Times:\n\tC2: " << (t_code/1000.0) << " ms\n\tA: " << (t_a/1000.0) << " ms\n";

		/*
			Free Memory
		*/
		free(I);
		for(bool* a : SC2) {
			free(a);
		}
	}
}

//rate
void test_run_many_2_auto_N_r(const bool& save_codes = false, const uint64_t& M_START = 10, const uint64_t& M_END = 100, const uint64_t& M_STEP = 10, const uint64_t& ITERATIONS = 1) {
	uint64_t M_NUMBER_STEPS = (M_END - M_START)/M_STEP + 1;
	const uint64_t qs[] = {4, 5, 7};
	const uint64_t qsize[] = {856, 11006, 2896222}; //Calculated with generate_A_files

	//CSV Header
	FILE* f_time;
	f_time = fopen(("output/tests/code_rate_" + std::to_string(M_START) + "_to_" + std::to_string(M_END) + "_at_" + std::to_string(M_STEP) + ".csv").c_str(), "w");
	fprintf_s(f_time, "q \\ M");
	for(uint64_t M = M_START; M <= M_END; M = M + M_STEP) {
		fprintf_s(f_time, ", %u", M);
	}
	fprintf_s(f_time, "\n");

	//For each q
	for(uint64_t Q = 0; Q < 2; Q++) {
		uint64_t q = qs[Q];
		uint64_t qq = q * q;

		simple_clock c;
		c.reset();
		matrix<uint64_t> A("data/A_q_" + std::to_string(q) + ".dat");
		auto time_load_A = c.microseconds();
		cout << "Time loading A matrix (q = " << q << "): " << (time_load_A/1000.0) << "ms\n";

		float* data = make_array<float>(M_NUMBER_STEPS, 0);

		fprintf_s(f_time, "%u", q);

		//Data Adquisition
		std::cout << "q: " << q << "\n";
		#pragma omp parallel for num_threads(22)
		for(uint64_t i = 0; i < M_NUMBER_STEPS; i++) {

			uint64_t M = M_START + M_STEP * i;

			for(uint64_t iter = 0; iter < ITERATIONS; iter++) {
				#pragma omp critical
				std::cout << "M: " << M << ", Iteration " << iter << "\n";

				CHHash ch(M, q, &A);
				data[i] += (float)(ch.C->rate());

				if(save_codes) {
					ch.C->C->save("output/tests/codes/q_" + std::to_string(q) + "_M_" + std::to_string(M) + "_N_" + std::to_string(CHHash::get_code_length(q, M, A.width)) + "_iteration_" + std::to_string(iter) + ".dat");
				}
			}

			data[i] /= ITERATIONS;
		}

		//Write data to files
		for(uint64_t i = 0; i < M_NUMBER_STEPS; i++) {
			fprintf_s(f_time, ", %f", data[i]);
		}
		fprintf_s(f_time, "\n");

		//Free memory
		free(data);
	}

	fclose(f_time);
}

//Calculates codes using optimal N for number of codewords
void test_run_many_2_auto_N(const bool& save_codes = false, const uint64_t& M_START = 10, const uint64_t& M_END = 100, const uint64_t& M_STEP = 10, const uint64_t& ITERATIONS = 1) {
	uint64_t M_NUMBER_STEPS = (M_END - M_START)/M_STEP + 1;
	const uint64_t qs[] = {4, 5, 7};
	const uint64_t qsize[] = {856, 11006, 2896222}; //Calculated with generate_A_files

	//CSV Header
	FILE* f_time;
	f_time = fopen(("output/tests/code_time_" + std::to_string(M_START) + "_to_" + std::to_string(M_END) + "_at_" + std::to_string(M_STEP) + ".csv").c_str(), "w");
	fprintf_s(f_time, "q \\ M");
	for(uint64_t M = M_START; M <= M_END; M = M + M_STEP) {
		fprintf_s(f_time, ", %u", M);
	}
	fprintf_s(f_time, "\n");

	//For each q
	for(uint64_t Q = 0; Q < 3; Q++) {
		uint64_t q = qs[Q];
		uint64_t qq = q * q;

		simple_clock c;
		c.reset();
		matrix<uint64_t> A("data/A_q_" + std::to_string(q) + ".dat");
		auto time_load_A = c.microseconds();
		cout << "Time loading A matrix (q = " << q << "): " << (time_load_A/1000.0) << "ms\n";

		float* times = make_array<float>(M_NUMBER_STEPS, 0);

		fprintf_s(f_time, "%u", q);

		//Data Adquisition
		std::cout << "q: " << q << "\n";
		#pragma omp parallel for num_threads(22)
		for(uint64_t i = 0; i < M_NUMBER_STEPS; i++) {

			uint64_t M = M_START + M_STEP * i;

			for(uint64_t iter = 0; iter < ITERATIONS; iter++) {
				#pragma omp critical
				std::cout << "M: " << M << ", Iteration " << iter << "\n";

				simple_clock c;
				c.reset();
				CHHash ch(M, q, &A);
				times[i] += c.microseconds();
			}

			times[i] /= ITERATIONS;

			//if(save_codes) {
			//	ch.C->C->save("output/tests/codes/q_" + std::to_string(q) + "_M_" + std::to_string(M) + "_N_" + std::to_string(CHHash::get_code_length(q, M, A.width)) + ".dat");
			//}
		}

		//Write data to files
		for(uint64_t i = 0; i < M_NUMBER_STEPS; i++) {
			fprintf_s(f_time, ", %f", times[i]/1000.0f);
		}
		fprintf_s(f_time, "\n");

		//Free memory
		free(times);
	}

	fclose(f_time);
}

//Calculates biggest N based on number of codewords (M is multiplied by 3 due to how C_1 works)
void test_size_N(const uint64_t& M_START = 10, const uint64_t& M_END = 1000000, const uint64_t& M_STEP = 1000) {
	uint64_t M_NUMBER_STEPS = (M_END - M_START)/M_STEP + 1;
	const uint64_t qs[] = {4, 5, 7};

	//CSV Header
	FILE* f_time;
	FILE* f_m;
	f_time = fopen("output/tests/N_time.csv", "w");
	f_m = fopen("output/tests/N.csv", "w");
	fprintf_s(f_time, "q \\ M");
	fprintf_s(f_m, "q \\ M");
	for(uint64_t M = M_START; M <= M_END; M = M + M_STEP) {
		fprintf_s(f_time, ", %u", M);
		fprintf_s(f_m, ", %u", M);
	}
	fprintf_s(f_time, "\n");
	fprintf_s(f_m, "\n");

	//For each q
	for(uint64_t Q = 0; Q < 3; Q++) {
		uint64_t q = qs[Q];
		uint64_t qq = q * q;

		simple_clock c;
		c.reset();
		matrix<uint64_t> A("data/A_q_" + std::to_string(q) + ".dat");
		auto time_load_A = c.microseconds();
		cout << "Time loading A matrix (q = " << q << "): " << (time_load_A/1000.0) << "ms\n";

		uint64_t* N = make_array<uint64_t>(M_NUMBER_STEPS);
		float* times = make_array<float>(M_NUMBER_STEPS);

		fprintf_s(f_time, "%u", q);
		fprintf_s(f_m, "%u", q);

		//Data Adquisition
		for(uint64_t i = 0; i < M_NUMBER_STEPS; i++) {
			//std::cout << "q: " << q << ", i: " << i << "/" << N_NUMBER_STEPS << "\n";

			uint64_t M = (M_START + M_STEP * i) * 3;

			simple_clock c;
			c.reset();
			N[i] = M < q ? 1 : CHHash::get_code_length(q, M, A.width);
			times[i] = c.microseconds();
		}

		//Write data to files
		for(uint64_t i = 0; i < M_NUMBER_STEPS; i++) {
			fprintf_s(f_m, ", %u", N[i]);
			fprintf_s(f_time, ", %f", times[i]/1000.0f);
		}
		fprintf_s(f_time, "\n");
		fprintf_s(f_m, "\n");

		//Free memory
		free(N);
		free(times);
	}

	fclose(f_time);
	fclose(f_m);
}

//Calculates biggest M based on a specific codeword length
void test_size_M(const uint64_t& N_START = 10, const uint64_t& N_END = 10000, const uint64_t& N_STEP = 10) {
	uint64_t N_NUMBER_STEPS = (N_END - N_START)/N_STEP + 1;
	const uint64_t qs[] = {4, 5, 7};

	//CSV Header
	FILE* f_time;
	FILE* f_m;
	f_time = fopen("output/tests/M_time.csv", "w");
	f_m = fopen("output/tests/M.csv", "w");
	fprintf_s(f_time, "q \\ N");
	fprintf_s(f_m, "q \\ N");
	for(uint64_t N = N_START; N <= N_END; N = N + N_STEP) {
		fprintf_s(f_time, ", %u", N);
		fprintf_s(f_m, ", %u", N);
	}
	fprintf_s(f_time, "\n");
	fprintf_s(f_m, "\n");

	//For each q
	for(uint64_t Q = 0; Q < 3; Q++) {
		uint64_t q = qs[Q];
		uint64_t qq = q * q;

		simple_clock c;
		c.reset();
		matrix<uint64_t> A("data/A_q_" + std::to_string(q) + ".dat");
		auto time_load_A = c.microseconds();
		cout << "Time loading A matrix (q = " << q << "): " << (time_load_A/1000.0) << "ms\n";

		uint64_t* M = make_array<uint64_t>(N_NUMBER_STEPS);
		float* times = make_array<float>(N_NUMBER_STEPS);

		fprintf_s(f_time, "%u", q);
		fprintf_s(f_m, "%u", q);

		//Data Adquisition
		for(uint64_t i = 0; i < N_NUMBER_STEPS; i++) {
			//std::cout << "q: " << q << ", i: " << i << "/" << N_NUMBER_STEPS << "\n";

			uint64_t N = N_START + N_STEP * i;

			simple_clock c;
			c.reset();
			M[i] = generate_qary_lemma_2_1_M_fast(qq, q, N, A.width, 3000000000);
			times[i] = c.microseconds();
		}

		//Write data to files
		for(uint64_t i = 0; i < N_NUMBER_STEPS; i++) {
			fprintf_s(f_m, ", %u", M[i]);
			fprintf_s(f_time, ", %f", times[i]/1000.0f);
		}
		fprintf_s(f_time, "\n");
		fprintf_s(f_m, "\n");

		//Free memory
		free(M);
		free(times);
	}

	fclose(f_time);
	fclose(f_m);
}

void test_run_many_2(const uint64_t& TARGET_N, const uint64_t& TARGET_M, uint64_t NUM_ITERATIONS = 100, bool try8 = false) {
	const uint64_t Q[] = {4, 5, 7, 8};

	for(uint64_t i = 0; i < (try8 ? 4 : 3); i++) {
		uint64_t q = Q[i];

		FILE* f_time;
		f_time = fopen(("output/output_" + std::to_string(NUM_ITERATIONS) + "_time_q_" + std::to_string(q) + ".csv").c_str(), "w");
		FILE* f_m;
		f_m = fopen(("output/output_" + std::to_string(NUM_ITERATIONS) + "_M_q_" + std::to_string(q) + ".csv").c_str(), "w");

		fprintf_s(f_time, "N \\ M");
		fprintf_s(f_m, "N \\ M");
		for(uint64_t M = 1; M < (TARGET_M + 1); M++) {
			fprintf_s(f_time, ", %u", M);
			fprintf_s(f_m, ", %u", M);
		}
		fprintf_s(f_time, "\n");
		fprintf_s(f_m, "\n");

		simple_clock c;
		c.reset();
		matrix<uint64_t> A("data/A_q_" + std::to_string(q) + ".dat");
		auto time_load_A = c.microseconds();

		cout << "Time loading A matrix (q = " << q << "): " << (time_load_A/1000.0) << "ms\n";

		matrix<float> M_all(TARGET_N, TARGET_M, -1); //[width][height]
		matrix<double> T_all(TARGET_N, TARGET_M, -1); //[width][height]

		//#pragma omp parallel for num_threads(20)
		for(uint64_t M = q; M < (TARGET_M + 1); M++) {
			for(uint64_t N = 1; N < (TARGET_N + 1); N++) {
				if(CHHash::get_code_length(q, M * 3, A.width) < (uint64_t)ceil(N / 4.0)) {
					continue;
				}

				#pragma omp critical
				std::cout << "M: " << M << ", N: " << N << ", optimal N: " << CHHash::get_code_length(q, M * 3, A.width) << "\n";

				float M_average = 0;
				double T_average = 0;

				for(uint64_t iter = 0; iter < NUM_ITERATIONS; iter++) {
					CHHash h(M, N, q, &A);

					if(h.C != nullptr) {
						M_average += h.C->k;
						T_average += h.time_code;
					}
				}

				//std::cout << (M-1) << ", " << (N-1) << " (" << M_all.width << ", " << M_all.height << ")\n";
				M_all.data[M - 1][N - 1] = M_average/NUM_ITERATIONS;
				T_all.data[M - 1][N - 1] = T_average/NUM_ITERATIONS;
			}
		}

	
		for(uint64_t N = 0; N < (TARGET_N); N++) {
			fprintf_s(f_m, "%u", (N + 1));
			fprintf_s(f_time, "%u", (N + 1));
			for(uint64_t M = 0; M < (TARGET_M); M++) {
				fprintf_s(f_m, ", %f", M_all.data[M][N]);
				fprintf_s(f_time, ", %f", T_all.data[M][N]);
			}
			fprintf_s(f_m, "\n");
			fprintf_s(f_time, "\n");
		}

		fclose(f_time);
		fclose(f_m);
	}
}

void test_run_many(const uint64_t& TARGET_N, const uint64_t& MAX_M, bool try8 = false) {
	const uint64_t Q[] = {4, 5, 7, 8};

	FILE* f;
	f = fopen("output.csv", "w");
	fprintf_s(f, "q, Iteration, Time for Concatenation (microseconds), Size of Codewords (N), Number of Codewords (M)\n");

	for(uint64_t i = 0; i < (try8 ? 4 : 3); i++) {
		uint64_t q = Q[i];

		simple_clock c;
		c.reset();
		matrix<uint64_t> A("data/A_q_" + std::to_string(q) + ".dat");
		auto time_load_A = c.microseconds();

		cout << "Time loading A matrix (q = " << q << ") (S(c2): " << A.width << "): " << (time_load_A/1000.0) << "ms\n";

		uint64_t N = TARGET_N, M_max = MAX_M;

		for(uint64_t iter = 0; iter < 100; iter++) {

			simple_clock c1, c2;

			//std::cout << "Calculating...\n";
			//cout << "q: " << q << ", N: " << N << "\n";
			
			c1.reset();
			q_ary* C2;//Inner Code
			if(q == 4)
				C2 = lemma_3_14_case_1(2, 3); 
			else if(q == 8)
				C2 = lemma_3_14_case_1(3, 3); 
			else
				C2 = lemma_3_12(q);
			auto time2 = c1.microseconds();
			
			//cout << "Data:\n";
			//cout << "\tC2 Generation: " << (time2/1000.0) << "ms\n";

			q_ary* C;
			c2.reset();
			C = generate_qary_lemma_2_3(N, M_max, C2, &A);
			auto timeC = c2.microseconds();

			//cout << "\tConcatenation: " << (timeC/1000.0) << "ms\n";
			//cout << "\tNumber of Codewords: " << C->k << "\n";

			fprintf_s(f, "%u, %u, %u, %u, %u\n", q, iter, timeC, C->n, C->k);

			delete C;
		}
	}

	fclose(f);
}