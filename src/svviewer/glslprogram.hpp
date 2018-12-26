
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

#include <GL/glew.h>
#include <map>
#include <string>
#include <iostream>

class GLSLProgram {
public:
	typedef enum {
		GEOMETRY_SHADER,
		VERTEX_SHADER,
		FRAGMENT_SHADER
	} TShaderType;
public:
	GLSLProgram(std::string programName = "");
	~GLSLProgram();
	bool setAndCompileShader(TShaderType shaderType, std::string shaderCode);
	bool loadAndCompileShader(TShaderType shaderType, std::string file);
	bool link();

	void setDefine(std::string token, std::string value);
	void setDefine(std::string token, int value);
	void setDefine(std::string token, unsigned int value);
	void setDefine(std::string token, float value);
	inline void setDefine(std::string token) { setDefine(token, ""); }
	int getNumDefines();

	inline GLuint getGLProgram() { return _programGLHandler; }
	inline GLint getUniformID(std::string name, bool isUniformBlock = false) {
		if (_uniforms.find(name) == _uniforms.end())
			initUniform(name, isUniformBlock);
		return _uniforms[name];
	}
	inline void setName(std::string name) { _name = name; }

	void clean();
	bool reBuild();

	// debug stuff
	void printUniforms();

private:
	bool initUniform(std::string name, bool isUniformBlock = false);
	bool loadCodeFile(std::string filename, TShaderType shaderType);
	bool compileShader(TShaderType shaderType);

private:
	std::string _name;
	std::string _vertCode, _geomCode, _fragCode;
	GLuint _vertShader, _geomShader, _fragShader;
	GLuint _programGLHandler;
	std::map< std::string, GLint > _uniforms;
	std::map< std::string, std::string > _defines;
};

