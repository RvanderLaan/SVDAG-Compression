
//=============================================================
//  This file is part of the SymVox (Symmetry Voxelization) software
//  Copyright (C) 2016 by CRS4 Visual Computing Group, Pula, Italy
//
//  For more information, visit the CRS4 Visual Computing Group 
//  web pages at http://vic.crs4.it
//
//  This file may be used under the terms of the GNU General Public
//  License as published by the Free Software Foundation and appearing
//  in the file LICENSE included in the packaging of this file.
//
//  CRS4 reserves all rights not expressly granted herein.
//  
//  This file is provided AS IS with NO WARRANTY OF ANY KIND, 
//  INCLUDING THE WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS 
//  FOR A PARTICULAR PURPOSE.
//=============================================================

#include "halton.hpp"

void Halton::testHalton(int n, int dim) {
	Halton h(dim);
	for (int i = 1; i <= n; ++ i) {
		std::cout << "halton " << i << " of dim " << dim << " = " << h.getNext() << std::endl;
	}
}


void Halton::testPrime() {
	int dim = 2;
	for (int i = 1; i < 10; ++ i) {
		std::cout << "prime for dim " << i << " " << Halton::findPrime(i) << std::endl;
	}
}


Halton::Halton(int dim): _index(0) {
	_base = findPrime(dim);
}


float Halton::getNext() {
	++ _index;

	float result = .0f;
	float fraction = 1.0f / (float)_base;
	int idx = _index;

	while (idx > 0) {
		int digit = idx % _base;
		result += fraction * (float)digit;

		idx  = (idx - digit) / _base;
		fraction /= (float)_base;
	}
	return result;
}


bool Halton::isPrime(int n) {
	bool isPrime = true;
	for (int i = 2; i < n; ++ i) {
		if (n % i == 0) {
			isPrime = false;
			break;
		}
	}
	return isPrime;
}
 

int Halton::findPrime(int idx) {
	if (idx < 1) {
		std::cerr << "error: cannot find " << idx << "th prime" << std::endl;
		return 0;
	}

	// only even prime numbers
	if (idx == 1) return 2;

	int number = 3;
	int primeFound = 1;

	while (1) {
		if (isPrime(number)) ++ primeFound;
		if (primeFound == idx) break;
		number += 2;
	}
	return number;
}


void HaltonSequence::testHalton(int n, int dim) {
	HaltonSequence h(dim);
	float *haltons = new float[dim];

	for (int i = 1; i <= n; ++ i) {
		h.getNext(haltons);
		std::cout << "halton " << i << " of dim " << dim << " = ";

		for (int j = 0; j < dim; ++ j)
			std::cout << haltons[j] << " ";

		std::cout << std::endl;
	}
}


HaltonSequence::HaltonSequence(int dim) {
	for (int i = 0; i < dim; ++ i) {
		_haltons.push_back(Halton(i + 1));
	}
}


void HaltonSequence::getNext(float *numbers) {
	const size_t dim = _haltons.size();
	for (size_t i = 0; i < dim; ++ i) {
		numbers[i] = _haltons[i].getNext();
	}
}


void HaltonSequence::getNext(float *samples, int n, int startIdx) {
	const int dim = (int)_haltons.size();
	for (int i = 0; i < n; ++ i) {
		getNext(samples + startIdx + i * dim);
	}
}