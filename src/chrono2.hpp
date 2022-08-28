#pragma once

#include <chrono>

struct simple_clock {
	std::chrono::steady_clock::time_point begin;

	bool running;

	void reset() {
		begin = std::chrono::steady_clock::now();
	}

	int64_t microseconds() {
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		return std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
	}
};