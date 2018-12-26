
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

#include <string>
#include <sl/cstdint.hpp>

static std::string byteToString(const sl::uint8_t byte) {
	std::string s;
	for (int i = 7; i >= 0; --i)
		s += ((1U << i) & byte) ? '1' : '0';
	return s;
}

static std::string byteToString(const sl::uint16_t byte) {
	std::string s;
	for (int i = 15; i >= 0; --i)
		s += ((1U << i) & byte) ? '1' : '0';
	return s;
}

static std::string byteToString(const sl::uint32_t byte) {
	std::string s;
	for (int i = 31; i >= 0; --i) {
		s += ((1U << i) & byte) ? '1' : '0';
		if (i == 8 || i == 16 || i == 24) s += '.';
	}
	return s;
}

static unsigned int bitCount(const sl::uint8_t &byte) { 
	return (byte * 01001001001ULL & 042104210421ULL) % 017;
}


static void writeBMP(const std::string filename, const int width, const int height, const sl::uint8_t * imgBuffer, bool vertInvert = false) {
	FILE * fp = fopen(filename.c_str(), "wb");
	if (!fp) {
		std::cout << "WARNING! writeBMP(): Couldn't open for write '" << filename << "'" << std::endl;
		return;
	}
	int filesize = 54 + 3 * width * height;

	unsigned char bmpfileheader[14] = { 'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0 };
	unsigned char bmpinfoheader[40] = { 40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0 };
	unsigned char bmppad[3] = { 0,0,0 };

	bmpfileheader[2] = (unsigned char)(filesize);
	bmpfileheader[3] = (unsigned char)(filesize >> 8);
	bmpfileheader[4] = (unsigned char)(filesize >> 16);
	bmpfileheader[5] = (unsigned char)(filesize >> 24);

	bmpinfoheader[4] = (unsigned char)(width);
	bmpinfoheader[5] = (unsigned char)(width >> 8);
	bmpinfoheader[6] = (unsigned char)(width >> 16);
	bmpinfoheader[7] = (unsigned char)(width >> 24);
	bmpinfoheader[8] = (unsigned char)(height);
	bmpinfoheader[9] = (unsigned char)(height >> 8);
	bmpinfoheader[10] = (unsigned char)(height >> 16);
	bmpinfoheader[11] = (unsigned char)(height >> 24);

	fwrite(bmpfileheader, 1, 14, fp);
	fwrite(bmpinfoheader, 1, 40, fp);
	for (int i = 0; i < height; i++)
	{
		if(vertInvert) 
			fwrite(imgBuffer + (width* i * 3), 3, width, fp);
		else
			fwrite(imgBuffer + (width*(height - i - 1) * 3), 3, width, fp);
		fwrite(bmppad, 1, (4 - (width * 3) % 4) % 4, fp);
	}
	fclose(fp);
}