
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

// Triangle - Box tester (by Tomas Akenine-M�ller
// see http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/)

#pragma once

#include <sl/fixed_size_vector.hpp>
#include <sl/fixed_size_point.hpp>

bool planeBoxOverlap(const sl::vector3d normal, const sl::point3d vert, const double maxbox);
bool testTriBox(const sl::point3d boxCenter, const double boxHalfSide, const sl::point3f *tri);
sl::point3f closestPointOnTri(const sl::point3d p, const sl::point3f *triangle);
