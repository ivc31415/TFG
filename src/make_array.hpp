#ifndef MAKE_ARRAY
#define MAKE_ARRAY

#include <stdlib.h>
#include <iostream>

/*
	Allocates space for an array of type 'T' and size 'size'
	Returns pointer to the array
*/
template <typename T>
T* make_array(unsigned int size) {
	T* p;
	p = (T*)malloc(size*sizeof(T));
	return p;
}

/*
	Allocates space for an array of type 'T' and size 'size'
	Returns pointer to the array
	Initializes all the array to the value 'initial_value'
*/
template <typename T>
T* make_array(unsigned int size, const T initial_value) {
	T* p;
	p = (T*)malloc(size*sizeof(T));
	for(unsigned int i = 0; i < size; i++)
		p[i] = initial_value;
	return p;
}

/* 
	Fills an array with a repeating value
*/
template <typename T>
void fill_array(T* data, unsigned int size, const T value) {
	for(unsigned int i = 0; i < size; i++) {
		data[i] = value;
	}
}

/*
	Copy an array from one to another
*/
template <typename T>
T* copy_array(unsigned int size, T* original) {
	T* p = make_array<T>(size);
	for(unsigned int i = 0; i < size; i++) {
		p[i] = original[i];
	}
	return p;
}

/*
	Print an Array
*/
template <typename T>
void _print_array(T* a, uint64_t size) {
	std::cout << "[";
	for(uint64_t i = 0; i < size; i++) {
		if(i != 0) std::cout << ", ";
		std::cout << a[i];
	}
	std::cout << "]";
}

/*
	Sort Array (Quicksort)
*/

template <typename T>
void sort_array(T* v, int64_t start, int64_t end) {
	if(start < end) {
		uint64_t p1 = v[end];
		uint64_t l = start - 1;
		for(uint64_t j = start; j < end; j++) {
			if(v[j] <= p1) {
				l++;
				uint64_t aux = v[l];
				v[l] = v[j];
				v[j] = aux;
			}
		}
		l++;
		uint64_t aux = v[l];
		v[l] = v[end];
		v[end] = aux;

		sort_array(v, start, l - 1);
		sort_array(v, l + 1, end);
	}
}

template <typename T>
void sort_array(T* v, uint64_t size) {
	sort_array(v, 0, (int64_t)(size - 1));
}

#endif
