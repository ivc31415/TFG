#pragma once

#include <stdint.h>
#include <iostream>
#include <vector>

#include "make_array.hpp"

#define ARRAY_VERTICAL 0
#define ARRAY_HORIZONTAL 1

template <typename T>
struct matrix {
	T** data = nullptr; //[row][column]
	uint64_t height;
	uint64_t width;

	//Construct a new matrix, without or with a default value
	matrix(uint64_t _height, uint64_t _width) {
		height = _height;
		width = _width;
		data = make_array<T*>(_width);
		for(uint64_t i = 0; i < _width; i++) {
			data[i] = make_array<T>(_height);
		}
	}
	matrix(uint64_t _height, uint64_t _width, T default_value) {
		height = _height;
		width = _width;
		data = make_array<T*>(_width);
		for(uint64_t i = 0; i < _width; i++)
			data[i] = make_array(_height, default_value);
	}

	//Construct matrix from file
	matrix(std::string file_path) {
		FILE* file;
		file = fopen(file_path.c_str(), "r");

		fscanf_s(file, "%llu", &height);
		fscanf_s(file, "%llu", &width);

		data = make_array<T*>(width);
		for(uint64_t i = 0; i < width; i++) {
			data[i] = make_array<T>(height);
			for(uint64_t j = 0; j < height; j++) {
				uint64_t d = 0;
				fscanf_s(file, "%llu", &d);
				data[i][j] = (T)d;
			}
		}

		fclose(file);
		
		//std::cout << "Loaded matrix from file '" << file_path << "'\n";
	}

	//Construct a matrix from specific columns from a second matrix
	matrix(matrix<T>* original, bool* mask) {
		//Get size of new matrix
		height = original->height;
		width = 0;
		for(uint64_t i = 0; i < original->width; i++)
			if(mask[i])
				width++;
		
		//Allocate matrix
		data = make_array<T*>(width);
		//for(uint64_t i = 0; i < width; i++) {
		//	data[i] = make_array<T>(height);
		//}

		//Copy data
		uint64_t j = 0;
		for(uint64_t i = 0; i < original->width; i++) {
			if(mask[i]) {
				data[j] = copy_array(height, original->data[i]);
				//for(uint64_t x = 0; x < original->height; x++)
					//data[j][x] = rand()%3;
				j++;
			}
		}
	}

	void save(std::string file_path) {
		FILE* file;
		file = fopen(file_path.c_str(), "w");

		fprintf_s(file, "%llu %llu\n", height, width);
		for(uint64_t i = 0; i < width; i++) {
			for(uint64_t j = 0; j < height; j++) {
				if(j != 0)
					fprintf_s(file, " ");
				fprintf_s(file, "%llu", data[i][j]);
			}
			fprintf_s(file, "\n");
		}

		fclose(file);
	}

	~matrix() {
		for(uint64_t i = 0; i < width; i++) {
			free(data[i]);
		}
		free(data);
	}
	
	//Retrieve i-th column of the matrix
	/*
	T* operator [] (int i) {
		return data[i];
	}*/
};

//output stream overload
template <typename T>
std::ostream& operator<<(std::ostream& os, const matrix<T>* m) {
	os << "(" << m->height << " x " << m->width << ")\n";
	for(uint64_t j = 0; j < m->height; j++) {
		for(uint64_t i = 0; i < m->width; i++) {
			os << m->data[i][j] << ",\t";
		}
		os << "\b\b\n";
	}

	return os;
}
template <typename T>
std::ostream& operator<<(std::ostream& os, const matrix<T>& m) {
	os << &m;
	return os;
}