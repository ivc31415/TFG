#pragma once

#include <stdint.h>
#include <iostream>

struct loading_bar {
	bool first = true;

	void update(uint64_t progress, uint64_t target) {
		if(!first) {
			std::cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
		}

		int p = (int)(20*((float)progress)/target);
		std::cout << "[";
		for(int i = 0; i < 20; i++)
			if(p > i)
				std::cout << "|";
			else
				std::cout << " ";

		std::cout << "]\n";

		first = false;
	}

	void reset() {
		first = true;
	}
};