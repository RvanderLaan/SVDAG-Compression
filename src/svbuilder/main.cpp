
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


#include <fstream>
#include <string.h>

#include <symvox/scene.hpp>
#include <symvox/geom_octree.hpp>
#include <symvox/encoded_svo.hpp>
#include <symvox/encoded_svdag.hpp>
#include <symvox/encoded_ssvdag.hpp>
#include <symvox/encoded_ussvdag.hpp>

#include "../symvox/cluster.hpp"

void printUsage() {
	printf(
		"\n"
		"Usage:\n"
		"      svbuilder input_model.obj numLevels numBuildSteps [--lossy or -l] [--cross-level-merging or -c] [<output.[svo | svdag | ussvdag | ssvdag]>\n"
		"Where:\n"
		"      input_model.obj: a 3D model in ASCII OBJ format\n"
		"      numLevels: levels of the octree to build (i.e. 10->1K^3, 13->8K^3...)\n"
		"      numBuildSteps:\n"
		"          if   0, raw build of the full octree, only 1 thread\n"
		"          if > 0, levels of the 'base octree' to build. Then the subtrees\n"
		"                  of the full children will be computed in parallel, 1 thread per CPU\n"
		"          TIP: usually 0 when less than 10 levels,\n"
		"               then 1 o 2 for 11-14 levels\n"
		"               then 3,4 or even 5 for higher levels\n"
		"               (depending on PC capability and complexity of the model)\n"
		"      <output.[svo | svdag | ussvdag | ssvdag]>: (optional) output files to be written\n"
		"\n"
	);
}

int main(int argc, char ** argv) {

//	cluster::TestSubgraph2("edges11-orig.txt");
//	return 0;

	printf(
		"\n"
		"===============================================================================\n"
		"==========   SymVox Builder (c)2016 ViC / CRS4   http://vic.crs4.it   =========\n"
		"===============================================================================\n"
	);

	if (argc < 4) {
		printUsage();
		exit(1);
	}

	std::string inputFile(argv[1]);
	int nLevels = atoi(argv[2]);
	int levelStep = atoi(argv[3]);

	printf(" MODEL: '%s'   [ %d levels, step %i ]   (%.0fK^3)\n", inputFile.c_str(), nLevels, levelStep, pow(2, nLevels) / 1024.f);
	printf("===============================================================================\n\n");

	sl::real_time_clock ck;
	ck.restart();

	Scene scene;

    bool isLas = false;

	if (strstr(inputFile.c_str(), ".obj") || strstr(inputFile.c_str(), ".OBJ")) {
		scene.loadObj(inputFile, true, false, false, false, true);
    } else if (strstr(inputFile.c_str(), ".las") || strstr(inputFile.c_str(), ".LAS")) {
        scene.loadLas(inputFile);
        isLas = true;
    } else if (strstr(inputFile.c_str(), ".laz") || strstr(inputFile.c_str(), ".LAZ")) {
        scene.loadLas(inputFile);
        isLas = true;
	} else {
		printf("Can't read input file '%s'. Only supported ASCII Obj files.\n", inputFile.c_str());
		exit(1);
	}

    bool lossy = false;
	bool multiLevel = false;
	bool exploitHiddenGeom = false;

	float lossyInflation = 2.0;
	float allowedLossyDiffFactor = 0.99;
	int includedNodeRefCount = 1;

	for (int i = 4; i < argc; ++i) {
		std::string arg = argv[i];
		if (arg == "--lossy" || arg == "-l") {
		    try {
				if (i + 1 < argc) {
					// Clamped between 1.2 and 10 (https://micans.org/mcl/man/mclfaq.html#toc-granularity)
		        	lossyInflation = std::stof(argv[i + 1]);
					lossyInflation = std::min(10.0f, std::max(1.2f, lossyInflation));
				}
		        if (i + 2 < argc) {
					// Clamped between 1 and 8, as >8 diff would mean ???
					allowedLossyDiffFactor = std::stof(argv[i + 2]);
					allowedLossyDiffFactor = std::min(1.0f, std::max(0.1f, allowedLossyDiffFactor));
				}
		        if (i + 3 < argc) {
					includedNodeRefCount = std::stoi(argv[i + 3]);
					includedNodeRefCount = std::min(1000000, std::max(1, includedNodeRefCount));
				}
		    } catch(std::invalid_argument&) {
		        printf(" --- WARNING: Could not parse lossy compression parameters, using defaults. ---\n");
		    }
			printf("Lossy parameters: inflation %.1f, allowed diff factor: %.2f, included node ref count: %i\n", lossyInflation, allowedLossyDiffFactor, includedNodeRefCount);
			lossy = true;
		}
		else if (arg == "--cross-level-merging" || arg == "-c")
			multiLevel = true;
        if (arg == "--hidden-geometry" || arg == "-h")
            exploitHiddenGeom = true;
	}

	printf("Lossy: %d, Cross-level: %d, Hidden geom: %d\n", lossy, multiLevel, exploitHiddenGeom);

	GeomOctree octree(&scene);

	//sl::aabox3d sceneBBoxD = sl::conv_to<sl::aabox3d>::from(scene.getAABB());
	auto minF = scene.getAABB()[0],
		 maxF = scene.getAABB()[1];

	sl::point3d minD = sl::point3d(minF[0], minF[1], minF[2]);
	sl::point3d maxD = sl::point3d(maxF[0], maxF[1], maxF[2]);
	sl::aabox3d sceneBBoxD = sl::aabox3d(minD, maxD);

	EncodedSVO svo;
	if (levelStep == 0) {
        if (!isLas) {
            octree.buildSVO(nLevels, sceneBBoxD, false, NULL, false);
        } else {
            octree.buildSVOFromPoints(inputFile, nLevels, sceneBBoxD, false, NULL);
        }
        svo.encode(octree);
		if (exploitHiddenGeom) {
            octree.hiddenGeometryFloodfill();
//            octree.toHiddenGeometryDAG();
		}

		octree.toDAG();
        
	}
	else {
		octree.buildDAG(nLevels, levelStep, sceneBBoxD, true);
	}

	// For OBJs with very large bboxes, rendering is bugging out. Fix: Rescale bbox
	float bboxSize = scene.getAABB()[0].distance_to(scene.getAABB()[1]);
	float maxBboxSize = 100000.f;
	if (bboxSize > maxBboxSize || bboxSize < 0.1) {
		std::cout << "Normalizing bbox; too small/large: " + std::to_string(bboxSize) + ". Initially: " << scene.getAABB() << std::endl;
		sl::vector3f newBboxMax = (scene.getAABB()[1] - scene.getAABB()[0]).ok_normalized();
		newBboxMax = maxBboxSize * newBboxMax; // Normalize * 1000
		std::cout << "Corrected BBOX: " << newBboxMax << std::endl;
		scene.setAABB(
			sl::aabox3f(
				sl::point3f(0, 0, 0),
				sl::point3f(newBboxMax[0], newBboxMax[1], newBboxMax[2])
			)
		);
		octree.resizeSceneBbox(scene.getAABB());
	}

#if 0 // DEBUG TEST
	if (!octree.checkIntegrity()) {
		printf("Something wrong!\n");
		return 1;
	}
#endif

    octree.initChildLevels();

//	octree.DebugHashing();

	// Prepare for saving file to disk
	std::string path = sl::pathname_directory(inputFile);
	std::string baseName = sl::pathname_base(sl::pathname_without_extension(inputFile));
	std::string basePath = path + "/" + baseName + "_" + std::to_string(nLevels);
	std::string sep = sl::pathname_directory_separators();

	int precisionVal = 1;
	std::string trimmedInfl = std::to_string(lossyInflation).substr(0, std::to_string(lossyInflation).find(".") + precisionVal + 1);
	std::string trimmedLossFact = std::to_string(allowedLossyDiffFactor).substr(0, std::to_string(allowedLossyDiffFactor).find(".") + precisionVal + 3);
		
	std::string paramStr = "-l_" + trimmedInfl + "_" + trimmedLossFact + "_" + std::to_string(includedNodeRefCount);

	// Save base SVDAG and ESVDAG in any case
	EncodedSVDAG svdag2;
	svdag2.encode(octree);
	svdag2.save(basePath + ".svdag");

	EncodedSSVDAG esvdag2;
	esvdag2.encode(octree);
	esvdag2.save(basePath + ".esvdag");

	if (lossy) {
		auto origDagTime = octree.getStats().toDAGTime; // since toDag is also called on the lossy dag, restore original dag time
        octree.toLossyDag(lossyInflation, allowedLossyDiffFactor, includedNodeRefCount);
		octree.getStats().toDAGTime = origDagTime;
	}

	// Encode conventional SVDAG
	EncodedSVDAG svdag;
	svdag.encode(octree);

	if (multiLevel) {
		// TODO: Clean this up later
		// This block will save both the single and multi level merged SVDAG
		svdag.save(basePath + "-single.svdag");

		octree.mergeAcrossAllLevels();
		EncodedSVDAG svdag2;
		svdag2.encode(octree);
		svdag2.save(basePath + "-multi.svdag");
	}

	EncodedSSVDAG esvdag;
	if (!multiLevel) esvdag.encode(octree);

    if (!multiLevel && !exploitHiddenGeom) {
        octree.toSDAG(false, false);
    }

	EncodedUSSVDAG ussvdag;
	if (!multiLevel) ussvdag.encode(octree);

	EncodedSSVDAG ssvdag;
	if (!multiLevel) ssvdag.encode(octree);

    EncodedSSVDAG psvdag;
    // psvdag.encode(octreeCopy);

    bool saveAll = true;

    if (saveAll) {
		std::string infix = "";
		if (lossy)
			infix += paramStr;

        svdag.save(basePath + infix + ".svdag");
        ussvdag.save(basePath + infix + ".ussvdag");
        ssvdag.save(basePath + infix + ".ssvdag");
        esvdag.save(basePath + infix + ".esvdag");
//        psvdag.save(basePath + ".psvdag");
    } else {
        for (int i = 4; i < argc; ++i) {
            std::string outputFile = argv[i];
            if     (strstr(outputFile.c_str(), ".svo") || strstr(outputFile.c_str(), ".SVO")) {
                if (levelStep == 0) svo.save(outputFile);
            }
            else if (strstr(outputFile.c_str(), ".svdag") || strstr(outputFile.c_str(), ".SVDAG")) {
                svdag.save(outputFile);
            }
            else if (strstr(outputFile.c_str(), ".ussvdag") || strstr(outputFile.c_str(), ".USSVDAG")) {
                ussvdag.save(outputFile);
            }
            else if (strstr(outputFile.c_str(), ".ssvdag") || strstr(outputFile.c_str(), ".SSVDAG")) {
                ssvdag.save(outputFile);
            }
            else if (strstr(outputFile.c_str(), ".psvdag") || strstr(outputFile.c_str(), ".PSVDAG")) {
                psvdag.save(outputFile);
            }
        }
    }

	sl::time_duration totalTime = ck.elapsed();

	GeomOctree::Stats stats = octree.getStats();

	// Pointless SVO stats are THEORICAL
	// meaning that, is calculated as 1 byte / node (ESVO paper codification)

	printf("\n========= RESULTS '%s' [%d levels] (%.0fK^3) =========\n", inputFile.c_str(), nLevels, pow(2, nLevels) / 1024.f);
	printf("Voxels:     %s\t(%zu)\n", sl::human_readable_quantity(stats.nTotalVoxels).c_str(), stats.nTotalVoxels);
	printf("SVO Nodes:  %s\t(%zu)\n", sl::human_readable_quantity(stats.nNodesSVO).c_str(), stats.nNodesSVO);
	printf("DAG Nodes:  %s\t(%zu)\n", sl::human_readable_quantity(stats.nNodesDAG).c_str(), stats.nNodesDAG);
	printf("SDAG Nodes: %s\t(%zu)\n", sl::human_readable_quantity(stats.nNodesSDAG).c_str(), stats.nNodesSDAG);
	printf("Pointerless SVO    : %s\t(%.3f bits/vox)\n", sl::human_readable_size(stats.nNodesSVO).c_str(), (8 * stats.nNodesSVO) / float(octree.getNVoxels()));
	if(levelStep == 0) printf("Encoded SVO       : %s\t(%.3f bits/vox)\n", sl::human_readable_size(svo.getDataSize()).c_str(), (8 * svo.getDataSize()) / float(octree.getNVoxels()));
	printf("Encoded SVDAG      : %s\t(%.3f bits/vox)\n", sl::human_readable_size(svdag.getDataSize()).c_str(),  (8 * svdag.getDataSize())   / float(octree.getNVoxels()));
	printf("Encoded ESVDAG     : %s\t(%.3f bits/vox)\n", sl::human_readable_size(esvdag.getDataSize()).c_str(), (8 * esvdag.getDataSize()) / float(octree.getNVoxels()));
	printf("Encoded USSVDAG    : %s\t(%.3f bits/vox)\n", sl::human_readable_size(ussvdag.getDataSize()).c_str(), (8 * ussvdag.getDataSize()) / float(octree.getNVoxels()));
	printf("Encoded SSVDAG     : %s\t(%.3f bits/vox)\n", sl::human_readable_size(ssvdag.getDataSize()).c_str(), (8 * ssvdag.getDataSize()) / float(octree.getNVoxels()));
	// printf("Encoded PSVDAG     : %s\t(%.3f bits/vox)\n", sl::human_readable_size(psvdag.getDataSize()).c_str(), (8 * psvdag.getDataSize()) / float(octree.getNVoxels()));
	printf("SSVDAG / DAG       : %.1f %%\n", 100.f * ssvdag.getDataSize() / (float)svdag.getDataSize());
	printf("SSVDAG / SVO       : %.1f %%\n", 100.f * ssvdag.getDataSize() / (float)stats.nNodesSVO);
	printf("SVO->SVDAG time    : %s\n", sl::human_readable_duration(stats.toDAGTime).c_str());
	printf("SVDAG->LSVDAG time : %s\n", sl::human_readable_duration(stats.toLSVDAGTime).c_str());
	printf("SVDAG->SSVDAG time : %s\n", sl::human_readable_duration(stats.toSDAGTime).c_str());
	printf("Total time         : %s\n", sl::human_readable_duration(totalTime).c_str());
	printf("===================================================================\n\n");

	FILE *baseStdout = stdout;
	// Print to console and to stats file
	for (int i = 0; i < 2; i++) {
		if (i == 1) {
			stdout = fopen("stats.txt", "a");
		}
		printf("%s, %d, lossy: %d\n", baseName.c_str(), nLevels, lossy);
		printf("#Voxels, %zu\n", stats.nTotalVoxels);
		printf(", SVDAG, ESVDAG, SSVDAG, SVO\n");
		printf("#nodes, %zu, '', %zu, %zu\n", stats.nNodesDAG, stats.nNodesSDAG, stats.nNodesSVO);
		printf("memory (bytes), %zu, %zu, %zu, %zu\n", svdag.getDataSize(), esvdag.getDataSize(), ssvdag.getDataSize(), stats.nNodesSVO);
		
		if (lossy) {
			printf("Construction times:\n");
			printf(",SVDAG, LSVDAG, Total,Hashing, SimNodes, Clustering\n");
			printf("time (ms), %zu, %zu, %zu, %zu, %zu, %zu\n", stats.toDAGTime.as_milliseconds(), stats.toLSVDAGTime.as_milliseconds(), totalTime.as_milliseconds(), stats.lHashing.as_milliseconds(), stats.lSimNodes.as_milliseconds(), stats.lClustering.as_milliseconds());

			printf("Lossy stats, %.2f, %.3f, %i\n", lossyInflation, allowedLossyDiffFactor, includedNodeRefCount);
			printf("TotalVoxDifference, #NodesIn, #ClustersOut, #edges\n");
			printf("%zu, %zu, %zu, %zu\n", stats.totalLossyVoxelDifference, stats.nClusteredNodes, stats.nClusters, stats.nEdges);
		}
		printf("\n\n");
	}
	fclose(stdout);
	stdout = baseStdout;

	return 0;
}
