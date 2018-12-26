
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


#pragma once

#include <iostream>
#include <vector>

class Halton {
public:
	static void testHalton(int n, int dim);
	static void testPrime();
	Halton(int dim);
	float getNext();
	void restart() { _index = 0; }

protected:
	static bool isPrime(int n);
    static int findPrime(int idx);

	unsigned int _index;
	int _base;
};		


class HaltonSequence {
public:
	HaltonSequence(int dim);
	void getNext(float *numbers);
	// Returns n halton numbers and adds it to samples from startIdx
	void getNext(float *samples, int n, int startIdx);
	void restart() { for (size_t i = 0; i < _haltons.size(); ++ i) _haltons[i].restart(); }

	static void testHalton(int n, int dim);

protected:
	std::vector<Halton> _haltons;
};
