#include "octree_parser.hpp"

// Parser for .octree file format: https://github.com/Forceflow/cpu_voxel_raycaster/blob/master/src/octree_build.cpp

octree_parser::octree_parser()
{

}

void octree_parser::readOctreeData(OctreeInfo const &octree_info, VoxelData** data) {
    string filename = octree_info.base_filename+string(".octreedata");
    FILE* file = fopen(filename.c_str(), "rb");

    *data = new VoxelData[octree_info.n_data];

    // read data
    for(size_t i = 0; i< octree_info.n_data; i++){
        (*data)[i] = VoxelData();
        readVoxelData(file,(*data)[i]);
    }
    fclose(file);
}

void octree_parser::readOctreeNodes(OctreeInfo const &octree_info, std::vector<OctNode> &nodes){
    string filename = octree_info.base_filename+string(".octreenodes");
    FILE* file = fopen(filename.c_str(), "rb");

    nodes.reserve(octree_info.n_nodes);

    for(size_t i = 0; i< octree_info.n_nodes; i++){
        OctNode n;
        n.read(file);
        nodes.push_back(n);
    }
    fclose(file);
}

// Parse a given octree header, store info in OctreeInfo struct
int octree_parser::parseOctreeHeader(const std::string &filename, OctreeInfo &i){
    cout << "  reading octree header from " << filename << " ... " << endl;
    ifstream headerfile;
    headerfile.open(filename.c_str(), ios::in);

    i.base_filename = filename.substr(0,filename.find_last_of("."));

    string line; bool done = false;
    headerfile >> line >> i.version;
    if (line.compare("#octreeheader") != 0) {cout << "    Error: first line reads [" << line << "] instead of [#octreeheader]" << endl; return 0;}

    while(headerfile.good() && !done) {
        headerfile >> line;
        if (line.compare("END") == 0) done = true; // when we encounter data keyword, we're at the end of the ASCII header
        else if (line.compare("gridlength") == 0) {headerfile >> i.gridlength;}
        else if (line.compare("n_nodes") == 0) {headerfile >> i.n_nodes;}
        else if (line.compare("n_data") == 0) {headerfile >> i.n_data;}
        else { cout << "  unrecognized keyword [" << line << "], skipping" << endl;
        char c; do { c = headerfile.get(); } while(headerfile.good() && (c != '\n'));
        }
    }

    headerfile.close();
    return 1;
}

unsigned int BinaryToGray(sl::uint8_t num) {
    return num ^ (num >> 1);
}

int octree_parser::buildAttrOctree(std::string basefilename){
    cout << "Reading octree from cache..." << endl;

    // compute inputfile
    size_t splitpoint = basefilename.find_last_of(".");
    stringstream octreecachefile, nodefile, datafile;
    octreecachefile << basefilename.substr(0,splitpoint) << ".octree";
    nodefile << basefilename.substr(0,splitpoint) << ".octreenodes";
    datafile << basefilename.substr(0,splitpoint) << ".octreedata";

    OctreeInfo info;
    parseOctreeHeader(basefilename, info);

    unsigned int _levels = sl::Log2_E(info.gridlength);

    vector<OctNode> octNodes;
    VoxelData* octData;
    readOctreeNodes(info, octNodes);
    readOctreeData(info, &octData);
    OctNode octRoot = octNodes[0];

    // Convert .octree format into GeomOctree
    //////////////////////////////////////////

    // TOdo: Create struct of GeomOctree data to pass into GeomOctree::BuildAttrSVO

    // TOdo: For attributes, increment _levels by one and start lev at 1. Insert 8 nodes instead of 1 (for each bit)
    GeomOctree::NodeData _data;
    _data.resize(_levels);

    // To keep track of the indices of children of OctNodes in the next level
    vector<size_t> childIndices, nextChildIndices;
    unsigned int lev = 0;

    // Initialize GeomOctree root node
    GeomOctree::Node root;
    _data[0].push_back(root);
    childIndices.push_back(0);

    // Convert the flat list of nodes into nodes per level as in GeomOctree._data
    for (unsigned int lev = 0; lev < _levels; ++lev) {
        int nodeIndex = 0;
        // For every node in the previous level, add its children to the current level
        for (auto node : _data[lev]) {
            // The corresponding node of the .octree format
            size_t childIndex = childIndices[nodeIndex];
            OctNode octNode = octNodes[childIndex];

            // Set its children in the GeomOctree
            for (int i = 0; i < 8; ++i) {
                char offset = octNode.children_offset[i];
                if (offset != -1) { // -1 means no child, otherwise 0 through to 7
                    node.setChildBit(i);
                    // Init children of this child - but no children in the final level
                    if (lev != _levels - 1) {
                        node.children[i] = _data[lev + 1].size();
                        _data[lev + 1].emplace_back();

                        size_t childIndex = octNode.children_base + offset;
                        nextChildIndices.push_back(childIndex);
                    }
                }
            }
        }
        childIndices.clear();
        childIndices.swap(nextChildIndices);
    }

    // compute NVoxels
    unsigned int _nVoxels = 0;
    for (id_t i = 0; i < _data[_levels - 1].size(); ++i) {
        _nVoxels += _data[_levels - 1][i].getNChildren();
    }

    // Update stats
//    _stats.nNodesSVO = info.n_nodes;
//    _stats.nNodesLastLevSVO = _data[_levels - 1].size();
//	_stats.simulatedEncodedSVOSize = (_stats.nNodesSVO - _stats.nNodesLastLevSVO) * 4;
//	_stats.nTotalVoxels = _nVoxels;
//    _stats.nTotalVoxels = 0;

//    _state = GeomOctree::S_SVO;

    // &(octree.getLevels()) = _levels;
    


//    // start reading octree
    // octree = new Octree(); // create empty octree
//    octree->gridlength = info.gridlength;
//    octree->n_data = info.n_data;
//    octree->n_nodes = info.n_nodes;

//    // read header and print statistics


    // Todo: Read octree
    // For one channel (e.g. red);
    // Construct GeomOctree for each bit separately
    // Join them together with a new root
    // Call toDAG
    // ?
    // Profit!

    for (int i = 0; i < 8; ++i) {
        
    }

    return 1;
}
