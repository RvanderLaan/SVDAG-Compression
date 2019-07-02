
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

void printUsage() {
	printf(
		"\n"
		"Usage:\n"
		"      svbuilder input_model.obj numLevels numBuildSteps <output.[svo | svdag | ussvdag | ssvdag]>\n"
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

	GeomOctree octree(&scene);

	sl::aabox3d sceneBBoxD = sl::conv_to<sl::aabox3d>::from(scene.getAABB());

	EncodedSVO svo;
	if (levelStep == 0) {
        if (!isLas) {
            octree.buildSVO(nLevels, sceneBBoxD, false, NULL, false);
        } else {
            octree.buildSVOFromPoints(inputFile, nLevels, sceneBBoxD, false, NULL);
        }
        svo.encode(octree);
        if (lossy) {
            octree.toLossyDAG();
        } else {
            octree.toDAG();
        }
	}
	else {
		octree.buildDAG(nLevels, levelStep, sceneBBoxD, true);
	}


#if 0 // DEBUG TEST
	if (!octree.checkIntegrity()) {
		printf("Something wrong!\n");
		return 1;
	}
#endif

    octree.initChildLevels();


	// TODO: Clean this up later
	// This block will save both the single and multi level merged SVDAG
    std::string path = sl::pathname_directory(inputFile);
    std::string baseName = sl::pathname_base(sl::pathname_without_extension(inputFile));
    std::string basePath = path + "/" + baseName + "_" + std::to_string(nLevels);

    // Save single-level merged
	EncodedSVDAG svdag;
	svdag.encode(octree);
    svdag.save(basePath + "-single.svdag");

    octree.mergeAcrossAllLevels();
    EncodedSVDAG svdag2;
    svdag2.encode(octree);
    svdag2.save(basePath + "-multi.svdag");





    if (!lossy) {
//        octree.toSDAG(false, false);
//        octreeCopy.toSDAG(false, true);
    }

//    octree.findAllSymDuplicateSubtrees();

	EncodedUSSVDAG ussvdag;
//	ussvdag.encode(octree);

	EncodedSSVDAG ssvdag;
//	ssvdag.encode(octree);

    EncodedSSVDAG psvdag;
//    psvdag.encode(octreeCopy);

    bool saveAll = false;

    if (saveAll) {
        std::string path = sl::pathname_directory(inputFile);
        std::string baseName = sl::pathname_base(sl::pathname_without_extension(inputFile));

        std::string basePath = path + "/" + baseName + "_" + std::to_string(nLevels);

        svo.save(basePath + ".svo");
        svdag.save(basePath + ".svdag");
        ussvdag.save(basePath + ".ussvdag");
        ssvdag.save(basePath + ".ssvdag");
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
	printf("Encoded SVDAG      : %s\t(%.3f bits/vox)\n", sl::human_readable_size(svdag.getDataSize()).c_str(),   (8 * svdag.getDataSize())   / float(octree.getNVoxels()));
	printf("Encoded USSVDAG    : %s\t(%.3f bits/vox)\n", sl::human_readable_size(ussvdag.getDataSize()).c_str(), (8 * ussvdag.getDataSize()) / float(octree.getNVoxels()));
	printf("Encoded SSVDAG     : %s\t(%.3f bits/vox)\n", sl::human_readable_size(ssvdag.getDataSize()).c_str(), (8 * ssvdag.getDataSize()) / float(octree.getNVoxels()));
	printf("Encoded PSVDAG     : %s\t(%.3f bits/vox)\n", sl::human_readable_size(psvdag.getDataSize()).c_str(), (8 * psvdag.getDataSize()) / float(octree.getNVoxels()));
	printf("SSVDAG / DAG       : %.1f %%\n", 100.f * ssvdag.getDataSize() / (float)svdag.getDataSize());
	printf("SSVDAG / SVO       : %.1f %%\n", 100.f * ssvdag.getDataSize() / (float)stats.nNodesSVO);
	printf("SVO->SVDAG time    : %s\n", sl::human_readable_duration(stats.toDAGTime).c_str());
	printf("SVDAG->SSVDAG time : %s\n", sl::human_readable_duration(stats.toSDAGTime).c_str());
	printf("Total time         : %s\n", sl::human_readable_duration(totalTime).c_str());
	printf("===================================================================\n\n");




	return 0;
}
