
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


#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

#include <GL/glew.h>

#include "glslprogram.hpp"


GLSLProgram::GLSLProgram(std::string programName) {
	_name = programName;
	_vertShader = _geomShader = _fragShader = _programGLHandler = 0;
}

GLSLProgram::~GLSLProgram() {
	if (_programGLHandler) {
		clean();
		glDeleteProgram(_programGLHandler);
	}
}

bool GLSLProgram::setAndCompileShader(TShaderType shaderType, std::string shaderCode)
{
	switch (shaderType) {
	case GEOMETRY_SHADER:
		_geomCode = shaderCode;
		break;
	case VERTEX_SHADER:
		_vertCode = shaderCode;
		break;
	case FRAGMENT_SHADER:
		_fragCode = shaderCode;
		break;
	}

	if (!compileShader(shaderType)) {
		return false;
	}

	return true;
}

bool GLSLProgram::loadAndCompileShader(TShaderType shaderType, std::string file) {
	
	if (!loadCodeFile(file, shaderType)) {
		std::cout << "WARNING GLSLProgram::loadAndCompileShader: Couldn't read " << file << " file." << std::endl;
		return false;
	}

	if (!compileShader(shaderType)) {
		return false;
	}

	return true;
}

bool GLSLProgram::loadCodeFile(std::string filename, TShaderType shaderType) {
	
	std::string code;
	std::ifstream ifs(filename, std::ios::in);

	if (ifs.is_open()) {
		std::string line = "";
		while (std::getline(ifs, line))
			code += line + " \n";
		ifs.close();
	}
	else return false;
	
	switch (shaderType) {
	case GEOMETRY_SHADER:
		_geomCode = code;
		break;
	case VERTEX_SHADER:
		_vertCode = code;
		break;
	case FRAGMENT_SHADER:
		_fragCode = code;
		break;
	}
	return true;
}

bool GLSLProgram::compileShader(TShaderType shaderType) {
	
	std::string code, shaderTypeStr;
	GLuint shaderHandler;
	switch (shaderType) {
	case GEOMETRY_SHADER:
		shaderHandler = glCreateShader(GL_GEOMETRY_SHADER);
		code = _geomCode;
		_geomShader = shaderHandler;
		shaderTypeStr = "Geometry Shader";
		break;
	case VERTEX_SHADER:
		shaderHandler = glCreateShader(GL_VERTEX_SHADER);
		_vertShader = shaderHandler;
		code = _vertCode;
		shaderTypeStr = "Vertex Shader";
		break;
	case FRAGMENT_SHADER:
		shaderHandler = glCreateShader(GL_FRAGMENT_SHADER);
		_fragShader = shaderHandler;
		code = _fragCode;
		shaderTypeStr = "Fragment Shader";
		break;
	}

	std::string definesCode = "\n";
	std::map<std::string, std::string>::iterator it;
	for (it = _defines.begin(); it != _defines.end(); it++) {
		std::string defineLine = "#define " + (*it).first + " " + (*it).second + "\n";
		definesCode += defineLine;
	}

	size_t versionPos = code.find("version");
	size_t insertDefinesPos = 0;
	if (versionPos != std::string::npos) {
		insertDefinesPos = code.find("\n", versionPos) + 1;
	}
	code.insert(insertDefinesPos, definesCode);

	char const * codePointer = code.c_str();
	//printf("%s\n", code.c_str());
	glShaderSource(shaderHandler, 1, &codePointer, NULL);
	glCompileShader(shaderHandler);

	GLint result = GL_FALSE;
	int  infoLogLenght = 0;

	glGetShaderiv(shaderHandler, GL_COMPILE_STATUS, &result);
	glGetShaderiv(shaderHandler, GL_INFO_LOG_LENGTH, &infoLogLenght);
	if (result == GL_FALSE)
	{
		if (infoLogLenght > 0)
		{
			std::cout << "GLSL ERROR (" << _name << ", " << shaderTypeStr << ")" << std::endl;
			GLchar* strInfoLog = new GLchar[infoLogLenght + 1];
			int dum;
			glGetShaderInfoLog(shaderHandler, infoLogLenght, &dum, strInfoLog);
			std::cout << strInfoLog << std::endl;
			std::cout << "------------------------------------------------" << std::endl;
			delete[] strInfoLog;
		}
	}

	return result == GL_TRUE ? true : false;
}

bool GLSLProgram::link() {

	// Link the program
	//std::cout << "* Loading GLSL program \"" << _name << "\"... ";
	if(!_programGLHandler) _programGLHandler = glCreateProgram();

	if(_geomShader) glAttachShader(_programGLHandler, _geomShader);
	if(_vertShader) glAttachShader(_programGLHandler, _vertShader);
	if(_fragShader) glAttachShader(_programGLHandler, _fragShader);

	glLinkProgram(_programGLHandler);

	GLint result;
	int infoLogLenght = 0;

	// Check the program

	glGetProgramiv(_programGLHandler, GL_LINK_STATUS, &result);
	glGetProgramiv(_programGLHandler, GL_INFO_LOG_LENGTH, &infoLogLenght);

	//if(result == GL_TRUE) std::cout << "OK!" << std::endl;
	if(result != GL_TRUE) {
		std::cout << "FAILED!!!" << std::endl;

		if ( infoLogLenght > 0 ){
			std::cout << "---------- GLSL Compiler Log -----------" << std::endl;
			std::vector<char> errorMessage(infoLogLenght+1);
			glGetProgramInfoLog(_programGLHandler, infoLogLenght, NULL, &errorMessage[0]);
			printf("%s\n", &errorMessage[0]);
			std::cout << "----------------------------------------" << std::endl;
		}
	}
	
	return result==GL_TRUE ? true : false;
}

void GLSLProgram::clean() {
	if (_geomShader) {
		glDetachShader(_programGLHandler, _geomShader);
		glDeleteShader(_geomShader);
	}
	if (_vertShader) {
		glDetachShader(_programGLHandler, _vertShader);
		glDeleteShader(_vertShader);
	}
	if (_fragShader) {
		glDetachShader(_programGLHandler, _fragShader);
		glDeleteShader(_fragShader);
	}
}

bool GLSLProgram::reBuild() {
	clean();
	if (!_geomCode.empty()) compileShader(GEOMETRY_SHADER);
	if (!_vertCode.empty()) compileShader(VERTEX_SHADER);
	if (!_fragCode.empty()) compileShader(FRAGMENT_SHADER);
	link();
	return true;
}

bool GLSLProgram::initUniform(std::string name, bool isUniformBlock) {
	GLint id;
	if (isUniformBlock)
		id = glGetUniformBlockIndex(_programGLHandler, name.c_str());
	else
		id = glGetUniformLocation(_programGLHandler, name.c_str());
	_uniforms[name] = id;
	if (id == -1) {
		std::cout << "WARNING GLSLProgram [" << _name << "] Uniform \"" << name << "\" does not exists or not used" << std::endl;
		glGetError();
		return false;
	}
	return true;
}

void GLSLProgram::setDefine(std::string token, std::string value) {
	_defines[token] = value;
}

void GLSLProgram::setDefine(std::string token, int value) {
	char str[256];
	sprintf(str, "%i", value);
	_defines[token] = std::string(str);
}

void GLSLProgram::setDefine(std::string token, unsigned int value) {
	char str[256];
	sprintf(str, "%u", value);
	_defines[token] = std::string(str);
}

void GLSLProgram::setDefine(std::string token, float value) {
	char str[256];
	sprintf(str, "%.20ff", value);
	_defines[token] = std::string(str);
}

int GLSLProgram::getNumDefines()
{
	return _defines.size();
}

void GLSLProgram::printUniforms() {
	printf("UNIFORMS [%s]  --------------\n", _name.c_str());
	for (std::map< std::string, GLint >::iterator i = _uniforms.begin(); i != _uniforms.end(); i++) {
		std::cout << (*i).first << " : " << (*i).second << std::endl;
	}
}
