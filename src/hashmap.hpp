#pragma once

#include <stdio.h>

#include <string>
#include "make_array.hpp"

template<typename T>
struct LCHashMapElement {
	std::string* key;
	T data;
	LCHashMapElement<T>* next_element = nullptr;
};

template<typename T>
class LCHashMap {
	private:
		LCHashMapElement<T>* data;
		uint64_t size;

		//Test basic hash function
		uint64_t hash_test(const std::string& key) {
			uint64_t h = 0;
			for(uint64_t i = 0; i < key.size(); i++) {
				h += key[i];
			}
			return h;
		}

		//Has function to use
		uint64_t hash(const std::string& key) {
			return this->hash_test(key);
		}

	public:
		LCHashMap(const uint32_t& size) {
			this->size = size;
			data = make_array<LCHashMapElement<T>>(size);
			for(uint32_t i = 0; i < size; i++) {
				data[i].key = nullptr;
				data[i].next_element = nullptr;
			}
		}

		//#warning remove temp counters in add() and get()
		//Return 0 if everything went correctly
		uint8_t add(const std::string& key, const T& value) {
			if(key.compare("") == 0) return 1;

			uint64_t i = this->hash(key)%size;
			auto element = &(data[i]);

			uint16_t TEMPCOUNTER = 0;
			if(element->key != nullptr) {
				do {
					element = element->next_element;
					TEMPCOUNTER++;
				} while(element->next_element != nullptr);
				element->next_element = make_array<LCHashMapElement<T>>(1);
				element = element->next_element;
			}

			element->key = new std::string(key);
			element->next_element = nullptr;
			element->data = value;

			//printf("Inserted key '%s' with value '%i' at index %u and subposition %u!\n", key.c_str(), element->data, i, TEMPCOUNTER);

			return 0;
		}

		T* get(const std::string& key) {
			uint64_t i = this->hash(key)%size;
			if(data[i].key == nullptr) return nullptr;

			auto element = &(data[i]);
			while(element->next_element != nullptr && element->key->compare(key)) { //While the key is not the same
				element = element->next_element;
			}

			return &(element->data);
		}
};