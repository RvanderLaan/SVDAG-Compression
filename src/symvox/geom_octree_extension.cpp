//
// Created by remi on 01-08-19.
// My own extensions to GeomOctree, to keep it separate from the original SSVDAG codebase
//

#include <memory>
#include <queue>
#include <stack>
#include <fstream>
#include <iomanip>      // std::setprecision

#include "geom_octree.hpp"
#include "cluster.hpp"

////////////////////////////////////////////////////////////////////
/////////////////////// HELPER FUNCTIONS ///////////////////////////
////////////////////////////////////////////////////////////////////

// Initializes the level for each child pointer to be in the next level
void GeomOctree::initChildLevels() {
    // Correcting child levels, in case multiple build steps were used
    for (unsigned int lev = 0; lev < _levels; ++lev) {
        for (id_t nodeIndex = 0; nodeIndex < _data[lev].size(); ++nodeIndex) {
            for (unsigned int c = 0; c < 8; ++c) {
                _data[lev][nodeIndex].childLevels[c] = lev + 1;
            }
        }
    }
}

/** Computes a uint64_t hash based on the child bitmasks of a node's children **/
std::hash <uint64_t> hash64;
uint64_t GeomOctree::computeNodeHash(const GeomOctree::Node &node, unsigned int depth) {
    uint64_t key = 0;

    // At depth 0, just return the child mask
    if (depth == 0) {
        key += node.childrenBitmask;
        return key;
    }
    // At depth 1, concat all 8-bit child masks into a 64 bit int
    else if (depth == 1) {
        for (unsigned int c = 0; c < 8; ++c) {
            // shift the bit mask of children into the key
            key = key << 8u;
            if (node.existsChildPointer(c)) {
                key += _data[node.childLevels[c]][node.children[c]].childrenBitmask;
            }
        }
    }
    // For higher depths, combine all children keys with XOR
    else {
        for (unsigned int c = 0; c < 8; ++c) {
            const auto &child = _data[node.childLevels[c]][node.children[c]];
            uint64_t origKeyHash = hash64(key);
            uint64_t childHash = origKeyHash << (c + 1u); // add some pseudo randomness when no child is at index c

            if (node.existsChildPointer(c)) {
                childHash = hash64(computeNodeHash(child, depth - 1));
                // TODO: reuse previous hash
//                childHash = matchMaps[lev][node.children[c]];
            }

            // bit shift so that identical child hashes don't cancel each other out
            key = (origKeyHash << 1u) ^ childHash; // xor of hash of current key and child key
            key = key ^ origKeyHash; // Hash with original key to make the hash truly asymmetrical
        }
    }
    return key;
};

uint64_t GeomOctree::computeNodeHashBotUp(const GeomOctree::Node &node, const unsigned int lev,
        const std::vector<std::vector<uint64_t>> &hashes, uint8_t childMaskMask) {
    uint64_t key = 0;
    if (lev == _levels - 1) {
        key += node.childrenBitmask & childMaskMask;
        return key;
    } else if (lev == _levels - 2) {
        for (unsigned int c = 0; c < 8; ++c) {
            // shift the bit mask of children into the key
            key = key << 8u;
            if (node.existsChildPointer(c) && (childMaskMask & (1U << c)) != 0) {
                key += _data[lev + 1][node.children[c]].childrenBitmask;
            }
        }
    } else {
        for (unsigned int c = 0; c < 8; ++c) {
            if (node.existsChildPointer(c) && (childMaskMask & (1U << c)) != 0) {
                uint64_t childKey = hashes[lev + 1][node.children[c]];

                key = key << 1u; // bit shift so that identical child hashes don't cancel each other out
                key = hash64(key) ^ hash64(childKey); // xor of hash of current key and child key
            }
        }
    }
    return key;
};

#define WRITE_KEYS_TO_FILE 0
/** Builds a multi-map of NodeHash -> Nodes with same hash */
void GeomOctree::buildMultiMap(unsigned int depth, std::vector<std::multimap<uint64_t, id_t>> &matchMaps, unsigned int levStart, unsigned int levEnd) {
    printf("[Match maps @D%u...", depth); fflush(stdout);

#if WRITE_KEYS_TO_FILE
    std::ofstream myfile("keys-" + std::to_string(depth) + ".txt", std::ios::out | std::ios::trunc);
#endif
    for (unsigned int lev = levStart; lev < levEnd - depth; ++lev) {
        matchMaps[lev].clear();
        for (id_t nodeIndex = 0; nodeIndex < _data[lev].size(); ++nodeIndex) {
            uint64_t key = computeNodeHash(_data[lev][nodeIndex], depth);
            matchMaps[lev].insert(std::make_pair(key, nodeIndex));
#if WRITE_KEYS_TO_FILE
                myfile << std::to_string(lev) + ", " + std::to_string(_data[lev][nodeIndex].childrenBitmask) + " ->\t " + std::to_string(key) + "\n";stdout
#endif
        }
    }
#if WRITE_KEYS_TO_FILE
    myfile.close();
#endif
    printf("]"); fflush(stdout);
};

/** Bottom-up version: Will compute hashes for all nodes in a bottom-up fashion, which should be faster, instead of specific levels as in the previous func */
void GeomOctree::buildMultiMapBotUp(std::vector<std::multimap<uint64_t, id_t>> &matchMaps,
        std::vector<std::vector<uint64_t>> &hashes, uint8_t childMaskMask) {
//    printf("[Match maps (BU)..."); fflush(stdout);
    for (int lev = _levels - 1; lev >= 0; --lev) {
        matchMaps[lev].clear();
        hashes.clear();
        hashes[lev].reserve(_data[lev].size());
        for (id_t nodeIndex = 0; nodeIndex < _data[lev].size(); ++nodeIndex) {
            // Only store hash if the childMaskMask completely overlaps with the inside part of the outsideMask
            const auto& node = _data[lev][nodeIndex];
            if ((node.outsideMask & childMaskMask) != 0) continue;

            int64_t hash = computeNodeHashBotUp(_data[lev][nodeIndex], lev, hashes, childMaskMask);
            matchMaps[lev].insert(std::make_pair(hash, nodeIndex));
            hashes[lev].push_back(hash);
        }
    }
//    printf("]"); fflush(stdout);
}

void GeomOctree::DebugHashing() {
    printf("DEBUG: Check how often hash collisions happen]\n");
    unsigned int histSize = 10; // num of collision counts to keep track of

    std::vector<std::multimap<uint64_t, id_t>> matchMaps(_levels);

    // csv header
    printf("Level, Depth");
    for (unsigned int i = 1; i < histSize; ++i) {
        printf(", %u", i);
    }
    printf(", Total hashes, Total nodes\n");

    for (unsigned int lev = _levels - 2; lev > 0; --lev) {
        std::vector<unsigned int> hist(histSize);

        printf("Level %u: ,", lev); fflush(stdout);

        // Compute node hashes for this node and its children up to the level above the leaves
        // Modify _levels - 1 to -2 to see how hashes collide when ignoring the leaf level
        unsigned int currentMatchDepth = std::max(_levels - lev - 2, 1u);
        buildMultiMap(currentMatchDepth, matchMaps, lev, lev + currentMatchDepth + 1);

        // Loop over every unique hash
        const auto& mm = matchMaps[lev];
        for(auto it = mm.begin(), end = mm.end(); it != end; it = mm.upper_bound(it->first)) {
            unsigned int count = mm.count(it->first);
            hist[std::min(count, (unsigned int) (histSize - 1))] += 1;
        }
//        printf("%u", lev);
        for (unsigned int i = 1; i < histSize; ++i) {
            printf(", %u", hist[i]);
        }
        printf(", %zu, %zu\n", mm.size(), _data[lev].size());
    }
}


/** Deep recursive subtree comparison **/
bool GeomOctree::compareSubtrees(
        unsigned int levA,
        unsigned int levB,
        Node &nA, // Should be in a lower level than b
        Node &nB, // Should be in a higher level than a
        std::vector<std::map<id_t, std::pair<unsigned int, id_t>>> &nodesInSubtree
) {
    // If B terminates before A, they are not equal since one or more levels are lost
    // Note: Subtree B is always higher up in the tree compared to subtree A (so a numerically lower level)
    if (nA.hasChildren() && !nB.hasChildren()) {
        return false;
    }
    else if (levA == _levels - 1) {
        // If nA is a leaf node, simply compare their child masks, since it doesn't matter what happens further down in B
        return nA.childrenBitmask == nB.childrenBitmask;
    }

    // Todo: Could maybe abort early on by looking up previously computed similar nodes to the child nodes of nA
    // If the the child of nB is not in that set of similar nodes to the child of nA, nA is not similar either (?)

    unsigned int childLevA = levA + 1;
    unsigned int childLevB = levB + 1;

    // Else, compare individual children: For every child
    for (int i = 0; i < 8; ++i) {

        // If the child bits don't match, they are not equal
        if (nA.existsChild(i) != nB.existsChild(i)) {
            return false;
        }
        // Otherwise the child mask bits are equal.
        // If there is no child node, they can be seen as equal since there is no subtree below this node
        // Note: This only holds for the leaf level.

        else if (!nA.existsChild(i) && !nB.existsChild(i)) {
            continue; // if they both have no children, continue looking through other children
        }
        // Else, compare child subtrees

        // Retrieve the child nodes from their respective levels
        Node &cA = _data[childLevA][nA.children[i]];
        Node &cB = _data[childLevB][nB.children[i]];

        // Add child node to the set of unique nodes in the subtree of node A
        nodesInSubtree[childLevA][nA.children[i]] = std::make_pair(childLevB, nB.children[i]); // stores the correspondence between these two nodes

        // Now compare the subtrees of these children - only if they are not equal, we can immediately return
        if (!this->compareSubtrees(childLevA, childLevB, cA, cB, nodesInSubtree)) {
            return false;
        }
    }
    return true;
}


/** Deep recursive subtree difference counter: Counts the difference in leaf nodes between two subtrees **/
void GeomOctree::diffSubtrees(
        unsigned int levA,
        unsigned int levB,
        const Node &nA, // Should be in a lower level than b
        const Node &nB, // Should be in a higher level than a,
        const unsigned int abortThreshold, // return before full comparison is finished if diff exceeds this number
        unsigned int &accumulator
) {
    // If nA is a leaf node, simply compare their child masks, since it doesn't matter what happens further down in B
    if (levA == _levels - 1) {
        for (int i = 0; i < 8; ++i) {
            if (nA.existsChild(i) != nB.existsChild(i))
                accumulator++;
        }
        return;
    }

    unsigned int childLevA = levA + 1;
    unsigned int childLevB = levB + 1;

    // Else, compare individual children: For every child
    for (int i = 0; i < 8; ++i) {

        // If the child bits don't match, compare...
        if (nA.existsChild(i) != nB.existsChild(i)) {
            if (levA == _levels - 2) {
                // If we're at the level above the leaves, simply check how different they are
                for (int j = 0; j < 8; ++j) {
                    if (nA.existsChild(i)
                        ? _data[levA + 1][nA.children[i]].existsChild(j)
                        : _data[levB + 1][nB.children[i]].existsChild(j)) {
                        accumulator++;
                    }
                }
                continue;
            } else {
                // If it's at a higher level, abort
                // Todo: Could still be a match in some weird edge cases, e.g. if a parent contains 1 child which contains 1 voxel
                accumulator += abortThreshold;
                return;
            }
        } else if (!nA.existsChild(i) && !nB.existsChild(i)) {
            // if they both have no children, continue looking through other children
            continue;
        }

        // Else, both nodes have children at this index: -> compare their subtrees
        const Node &cA = _data[childLevA][nA.children[i]];
        const Node &cB = _data[childLevB][nB.children[i]];

        this->diffSubtrees(childLevA, childLevB, cA, cB, abortThreshold, accumulator);
    }
}

/** Checks whether two subtrees at different levels are equal under a specific symmetry similarity */
bool GeomOctree::compareSymSubtrees(unsigned int levA, unsigned int levB, Node &nA_in, Node &nB, bool sX, bool sY, bool sZ) {
    // Still a WIP
    // Mirror node A
    Node nA = nA_in.mirror(sX, sY, sZ);
    invertInvs(nA, levA, sX, sY, sZ);

    // Same as normal subtree comparison
    if (nA.hasChildren() && !nB.hasChildren()) {
        return false;
    } else if (!nA.hasChildren() || !nB.hasChildren()) {
        // If they both don't have children, simply compare their child masks
        return nA.childrenBitmask == nB.childrenBitmask;
    }

    unsigned int childLevA = levA + 1;
    unsigned int childLevB = levB + 1;

    // If the end of the graph is reached, they are equal (?)
    if (childLevB >= _levels || childLevA >= _levels) {
        return true;
    }

    for (int i = 0; i < 8; ++i) {
        // If they child bits don't match, they are not equal
        if (nA.existsChild(i) != nB.existsChild(i)) {
            return false;
        }
            // Otherwise they are equal.
            // If either one is not set, we cannot compare them further. Therefore they can be seen as equal
        else if (!(nA.existsChild(i) || nB.existsChild(i))) {
            continue;
        }

        Node &cA = _data[childLevA][nA.children[i]];
        Node &cB = _data[childLevB][nB.children[i]];

        // Compare their child masks...
        if (cA.childrenBitmask != cB.childrenBitmask) {
            return false;
        } else {
            // If both nodes have children...
            if (cA.hasChildren() && cB.hasChildren()) {

                // Then compare the subtrees of these children
                if (!this->compareSymSubtrees(childLevA, childLevB, cA, cB, sX, sY, sZ)) {
                    return false;
                }
            }
            // If only 1 node has children, they can be seen as equal. One has more detail than the other
        }
    }
    return true;
}

/** Sorts nodes so that the node with the most refences is at index 0. Also updates pointers accordingly */
std::vector<std::vector<unsigned int>> GeomOctree::sortByRefCount() {
    /** Stolen from encoded_ssvdag.cpp */
    // histogram vector with pairs of idx, number of pointers to the nodes
    std::vector<std::pair<id_t, unsigned int>> hist;

    // indirection vector to map between one level indices (position in the vector) and their children offsets
    std::vector<GeomOctree::id_t> indirection;

    std::vector<std::vector<unsigned int>> refCounts(_levels);

    printf("Sorting nodes based on ref count: Level "); fflush(stdout);

    for (int lev = _levels - 1; lev >= 0; --lev) {
        printf(" %u..", lev); fflush(stdout);

        hist.resize(_data[lev].size());
        hist.shrink_to_fit();

        // refs initialization to ordered indices and zeros
        for (GeomOctree::id_t i = 0; i < _data[lev].size(); ++i) {
            hist[i].first = i; // idx
            hist[i].second = 0; // num of refs
        }

        if (lev > 0) { // Don't do it for the root node
            // count num of references in the superior level
            for (const GeomOctree::Node &n : _data[lev - 1])
                for (int c = 0; c < 8; ++c)
                    if (n.existsChildPointer(c)) hist[n.children[c]].second++;

            // sort by number of references (c++11's lambdas r00lez) ;D
            std::sort(hist.begin(), hist.end(),
                      [](std::pair< id_t, id_t > a, std::pair< id_t, id_t> b) { return a.second > b.second; }
            );
        }

        std::vector<id_t> newIndirection(_data[lev].size());
        std::vector<GeomOctree::Node> sortedNodes(_data[lev].size());

        refCounts[lev].reserve(_data[lev].size());

        for (id_t i = 0; i < hist.size(); ++i) {

            if (lev != _levels - 1) {
                GeomOctree::Node &n = _data[lev][hist[i].first];
                for (int c = 7; c >= 0; --c) {
                    if (n.existsChild(c)) {
                        n.children[c] = indirection[n.children[c]];
                    }
                }
            }

            newIndirection[hist[i].first] = i;
            sortedNodes[i] = _data[lev][hist[i].first];
            refCounts[lev][i] = hist[i].second;
        }

        indirection = newIndirection;

        _data[lev].clear();
        _data[lev].shrink_to_fit();
        _data[lev] = sortedNodes;
    }
    printf(" Done!\n");
    return refCounts;
}


/** Sorts nodes so that the node with the most refences is at index 0. Also updates pointers accordingly */
std::vector<std::vector<unsigned int>> GeomOctree::sortByEffectiveRefCount() {
    printf("Sorting nodes based on effective ref count: Level "); fflush(stdout);
    std::vector<std::vector<unsigned int>> refCounts(_levels);
    refCounts[0].emplace_back(1);

    // Pre-count effective ref counts
    for (int lev = 1; lev < _levels; ++lev) {
        refCounts[lev].assign(_data[lev].size(), 0);

        // count num of references in the superior level
        for (id_t nodeIndex = 0; nodeIndex < _data[lev-1].size(); ++nodeIndex) {
            const Node &n = _data[lev - 1][nodeIndex];
            // Add the ref counts from the parent to each child. -> If same child is referenced multiple times, add multiple times
            for (int c = 0; c < 8; ++c)
                if (n.existsChildPointer(c)) refCounts[lev][n.children[c]] += refCounts[lev - 1][nodeIndex];
        }
    }

    /** Stolen from encoded_ssvdag.cpp */
    // histogram vector with pairs of idx, number of pointers to the nodes
    std::vector<std::pair<id_t, unsigned int>> hist;

    // indirection vector to map between one level indices (position in the vector) and their children offsets
    std::vector<GeomOctree::id_t> indirection;

    // Sort
    for (int lev = _levels - 1; lev >= 0; --lev) {
        printf(" %u..", lev); fflush(stdout);

        hist.resize(_data[lev].size());
        hist.shrink_to_fit();

        // refs initialization to ordered indices and zeros
        for (GeomOctree::id_t i = 0; i < _data[lev].size(); ++i) {
            hist[i].first = i; // idx
            hist[i].second = refCounts[lev][i]; // num of refs
        }

        if (lev > 0) { // Don't do it for the root node
            // sort by number of references (c++11's lambdas r00lez) ;D
            std::sort(hist.begin(), hist.end(),
                      [](std::pair< GeomOctree::id_t, GeomOctree::id_t > a, std::pair< GeomOctree::id_t, GeomOctree::id_t> b) { return a.second > b.second; }
            );
        }

        std::vector<id_t> newIndirection(_data[lev].size());
        std::vector<GeomOctree::Node> sortedNodes(_data[lev].size());

        for (id_t i = 0; i < hist.size(); ++i) {

            if (lev != _levels - 1) {
                GeomOctree::Node &n = _data[lev][hist[i].first];
                for (int c = 7; c >= 0; --c) {
                    if (n.existsChild(c)) {
                        n.children[c] = indirection[n.children[c]];
                    }
                }
            }

            newIndirection[hist[i].first] = i;
            sortedNodes[i] = _data[lev][hist[i].first];
            refCounts[lev][i] = hist[i].second;
        }

        indirection = newIndirection;

        _data[lev].clear();
        _data[lev].shrink_to_fit();
        _data[lev] = sortedNodes;
    }
    printf(" Done!\n");


#if 1
    printf("DEBUG: Check how often nodes are referenced [effective ref count]\n");
    unsigned int histSize = 10; // num of ref counts to keep track of

    // csv header
    printf("Level");
    for (unsigned int i = 1; i < histSize; ++i) {
        printf(", %u", i);
    }
    printf(", Total nodes\n");

    for (unsigned int lev = 1; lev < _levels; ++lev) {
//        printf("- LEVEL %u, total nodes: %zu\n", lev, _data[lev].size());
        std::vector<unsigned int> hist(histSize);
        for (id_t nodeIndex = 0; nodeIndex < _data[lev].size(); ++nodeIndex) {
            hist[std::min(refCounts[lev][nodeIndex], (unsigned int) (histSize - 1))] += 1;
        }
        printf("%u", lev);
        for (unsigned int i = 1; i < histSize; ++i) {
            printf(", %u", hist[i]);
        }
        printf(", %zu\n", _data[lev].size());
    }
#endif

    return refCounts;
}

/**
 * Note: This is quite old and unusable. Should be merged into compareSymSubtrees
 * @brief GeomOctree::findAllSymDuplicateSubtrees Same as findAllDuplicateSubtrees, but with
 * symmetry as well: Find symmetrically identical subtrees over all levels
 * @return
 */
unsigned int GeomOctree::findAllSymDuplicateSubtrees() {
    // Still a WIP

    struct MirroredNode{
        MirroredNode() : mirrorX(false), mirrorY(false), mirrorZ(false), id(0) {}
        MirroredNode(bool mX, bool mY, bool mZ, id_t id) :
                mirrorX(mX), mirrorY(mY), mirrorZ(mZ), id(id) {}
        bool mirrorX, mirrorY, mirrorZ;
        id_t id;
    };

    // Note: Comparing ALL combinations of subtree pairs is not what we want
    // If a subtree high-up is equal to another subtree, all of its subtrees will equal it as well

    // Solution(?): Keep track of children that are NOT equal to another subtree
    // Only compare those instead of all nodes in the next level
    std::vector<Node> currentNodesToCheck;
    std::vector<Node> nextNodesToCheck;

    // Initialize all nodes of the highest level for comparison to subtrees at other levels
    unsigned int levA = 1;
    for (id_t i = 0; i < _data[levA].size(); i++) {
        Node &n = _data[levA][i];
        currentNodesToCheck.push_back(n);
    }

    unsigned int numEqualSubtrees = 0;
    unsigned int numTotalComparisons = 0;

    // Most efficient methinks:
    // Start at level 1 (A), compare to all nodes at levels above it (B), repeat for following levels
    // This way, the matches are found as early as possible

    // For every node to check...
    // Check for any level if there is an equal subtree
    // Note: If node B terminates before node A, they are NOT equal. The other way around is fine
    // For nodes that are not equal, add their children in the next level to the nextNodesToCheck list!

    std::vector<int> nodesEqualPerLevel(_levels, 0);

    for (; levA < _levels; ++levA) {
        printf("Comparing level %u\n", levA);
        // For all nodes to be checked...
        for (auto& nA : currentNodesToCheck) {

            bool foundMatch = false;
            for (unsigned int levB = 0; levB < levA; ++levB) {
                // For all nodes in the other level above A...
                for (id_t j = 0; j < _data[levB].size(); j++) {
                    numTotalComparisons++;

                    Node &nB = _data[levB][j];

                    // TODO: Compare all mirrored nodes of subtree nA to original nodes of subtree nB
                    // * Compute mirrored subtrees of nA
                    // * Compare each of them as before with nB
                    // * Repeat!

                    unsigned int nodesInSubtree = 1;

                    // Brute force check equality of subtrees of nA and nB
                    if (compareSymSubtrees(levA, levB, nA, nB, false, false, false)
                        || compareSymSubtrees(levA, levB, nA, nB, true, false, false)
                        || compareSymSubtrees(levA, levB, nA, nB, false, true, false)
                        || compareSymSubtrees(levA, levB, nA, nB, false, false, true)
                        || compareSymSubtrees(levA, levB, nA, nB, true, true, false)
                        || compareSymSubtrees(levA, levB, nA, nB, true, false, true)
                        || compareSymSubtrees(levA, levB, nA, nB, false, true, true)
                        || compareSymSubtrees(levA, levB, nA, nB, true, true, true)
                            ) {
                        foundMatch = true;
                        numEqualSubtrees++;
                        nodesEqualPerLevel[levA]++;
                        break;
                    }
                }
                if (foundMatch) {
                    // If a match is found, we are done for this subtree
                    // It is equal to a subtree higher up in the graph!
                    break;
                }
            }
            if (!foundMatch) {
                // If no match was found, try later for all child nodes
                for (int i = 0; i < 8; ++i) {
                    if (nA.existsChild(i)) {
                        nextNodesToCheck.push_back(_data[levA + 1][nA.children[i]]);
                    }
                }
            }
        }

        // Since children of a node may also be children of other nodes in a DAG,
        // we need to ensure children are only present once to the nextNodesToCheck vector
        std::sort(nextNodesToCheck.begin(), nextNodesToCheck.end() );
        nextNodesToCheck.erase(std::unique(nextNodesToCheck.begin(), nextNodesToCheck.end() ), nextNodesToCheck.end());

        // After all nodes for this level have been checked, swap the current and next nodes to check
        currentNodesToCheck.clear();
        currentNodesToCheck.swap(nextNodesToCheck);
    }

    printf("Nodes equal to another one: %u. Total #nodes %zu. Total #comparisons: %u\n", numEqualSubtrees, _nNodes, numTotalComparisons);

    for (unsigned int i = 0; i < _levels; ++i) {
        double pct = 100 * (nodesEqualPerLevel[i] / double(_data[i].size()));
        printf("- Level %u:   \t%i / %zu (%.2f%%) subtrees are equal to a subtree higher up\n", i, nodesEqualPerLevel[i], _data[i].size(), pct);
    }
    printf("Total equal: %u / %zu (%.2f%%)\n", numEqualSubtrees, _nNodes, (100 * numEqualSubtrees / double(_nNodes)));

    return numEqualSubtrees;
}

////////////////////////////////////////////////////////////////////
/////////////////////// LOSSY COMPRESSION //////////////////////////
////////////////////////////////////////////////////////////////////

/**
 * ToLossyDag1 idea:
 * - Merge nodes that are similar (e.g. only have 1 child difference)
 * - Most effective to perform on most common nodes
 *
 * Other implementation ideas
 * - Represent every child bit as a float, how close it is to the true geometry
 *
 * Proof of concept, greedy: Merge the largest sets of nodes that have diff 1
 * - Just per level for now, not earlier levels
 *
 * Looking at #nodes compressed per level w/ lossless compression, lossy compression should be avoided on:
 * - the leaves, there are very little benefits with large visual side effects
 * - the first half of levels, since this will have the largest side effects
 * -> therefore, focus on the deepest levels above the leaves
 */
/**
 * What needs to happen
 * - Convert to DAG as usual
 *   - This gives us a list of correspondences, which points to the new index in uniqueNodes given an index of _data[lev]
 * - Sort uniqueNodes by ref count
 * - Create new list of lossyUniqueNodes, of nodes that have correspondence to itself
 * - Find matches with a low diff (loop over every node for every node (N^2)
 *   - Replace list of correspondences, which points to nodes in lossyUniqueNodes given an index in _data[lev]
 * - Update pointers in level above
 */


/**
 * ToLossyDag2 Idea;
 * - Given a target percentage, like jpeg compression, merge as many closest matches until that target is reached
 * - Similar matches are found by comparing subtrees up to one or more levels before their maximum depth,
 *   and then the best match is found by performing deep comparisons
 * - Most effective on least common nodes, so it will act like outlier elimination
 *   Also had the idea of noisiness detection, but that's likely correlated with being rarely referenced, so not needed for now
 *
 * - The algorithm should probably start at the top of the graph, else it could cause a cascade of lossy errors
 *
 * - The target could be relative to the amount of nodes that is merged
 * - Start at the top of the graph, going down, merge all subtrees with diff 1. Then repeat for diff 2, 3, etc. until target is reached
 * - This means that up until diff 8, only one level should be ignored in finding match candidates
 *   Maybe also take into account whether all those diffs occur in the same node, or all in different ones?
 *
 * Other option
 *  -
 *
 * Other ideas:
 * - Generating new in-between node for two or more outliers, as a more accurate lossy representative
 *
 */


void GeomOctree::toLossyDag3() {
    // Reimplementation of toLossyDag2
    // Before starting anew:
    /*
     * Why doesn't toLossyDag2 work?
     * Problem: Too many detail is lost, even though diff is 1: E.g. large holes appear
     * - Hypothesis 1: Cascade of lossy errors
     * - How to (dis)prove this: Only perform on a single level
     * - Result: [experiment]
     *      - In some cases, it looks like the diff is larger than 1
     *      - Looks like the diff checking is happening at one level too high
     *      - No, I'm certain: There are nodes appearing where it was completely empty before.
     *        Even in the level above, the diff is at least 2 in this particular case, while the max was 1...
     *        Something smells fishy...
     *        Oh no.. Could it be that the pointers are updated one level higher than they're supposed to?
     *        Hmm, in other cases if works as intended: Only 1 diff at correct level
     * - How to resolve this: "Freeze" subtrees that are a match target.
     *      - Risk: Not taking into account ref counts of nodes in the subtree, only the root
     *
     * - Hypothesis 2: Pointer update mistake
     * -
     */

    ///////////////////////////////////////
    // From scratch:
    ///////////////////////////////////////
    /**
     * ToLossyDag3 idea: Find nodes in the DAG that are similar, merge them together, until a target number is reached
     * - Nodes are similar if the children in the leaf level are almost identical
     * - Removing nodes that are referenced frequently will result in more loss than those that are infrequently referenced
     * - - The reference count must take into account how often their parents are referenced as well, not just the individual node
     * - Merging two nodes in this manner might result in some child-nodes ending up unreferenced -> indirect gains
     * - - Q: can this be taken into account in the merging process?
     * - A node of which its parent was merged should be merged, as that will cause a cascade of lossy errors (more loss than intended)
     * - - Therefore, prefer to merge nodes with a node that has many references over one that only has a few
     */

    // Sort nodes based on effective number of references

    /*
     * Approach for finding similar nodes:
     * - Top-down/bottom-up: Doesn't matter, as long as it's infrequently referenced
     * - - Easier to do top-down: Then you can keep track of which children can be skipped, when already merged
     * - - Though, might not be needed when using effective ref count: Nodes with > 1 ref are skipped
     * - Comparison method: Compute hash based on all children, except the leaf level
     * - - For the candidate matches, compare leaf child masks to find match with least difference
     * - Preventing lossy error cascade
     * - - Check whether a node is a match before attempting to merge
     */

    /*
     * Pointer updating:
     * Option 1: Filter out nodes with 0 references
     * Option 2: Create new list, put all children in the list but check whether they already are in there first
     * Either case: Keep track of correspondences
     */

    // Maybe perform similarly to toDAG():
    // - For each level, merge and update pointers immediately afterwards (per level)

    /**
     * NEW new idea:
     * - For each level, starting at bottom:
     *   - Find all candidates for nodes with 1 ref to other nodes with 1 ref
     *      (note: could also be nodes with 2, 3, 4 refs)
     *   - Compute node clusters (only for nodes with few references)
     *   - Make list of unique nodes (center of clusters)
     *   - If a cluster center is similar to a node with more references, use that one
     *   - Replace duplicate nodes
     *   - Update parent pointers
     */

    if (_state == S_DAG) {
        printf("* Transforming DAG -> Lossy DAG3 ... \n");
    } else {
        printf("ERROR! This is not a DAG!\n");
        return;
    }

    auto refCounts = this->sortByEffectiveRefCount();

    _nNodes = 1;
    /** Every index denotes the index of the first duplicate of that node in uniqueNodes. Reset for each level. */
    std::vector<id_t> correspondences;
    std::vector<Node> uniqueNodes;
    std::vector<std::multimap<uint64_t, id_t>> matchMaps(_levels);
    std::map<id_t, bool> uniqueNodesChecker;

    sl::time_point ts = _clock.now();

    unsigned int nMatchesTotal = 0;

    // For every level, starting at two levels above the leaves...
    for (unsigned int lev = _levels - 2; lev > 0; --lev) {
		_clock.restart();
		sl::time_duration timeStamp = _clock.elapsed();
		int stepLogger = (int)round(_data[lev].size() / 10.0f);

		// Todo: Should depend on how many leaf nodes there are (missing) on the source node
        const int lossyDiff = (int) (_levels - lev) * 2; // should be outside variable: maximum allowed lossy diff

        std::vector<cluster::Edge> edges;

        // Clear the lists used to keep track of correspondences etc
        size_t oldLevSize = _data[lev].size();
        uniqueNodes.clear();
        uniqueNodes.shrink_to_fit();
        correspondences.clear();
        correspondences.resize(oldLevSize, (id_t)-1);
        uniqueNodesChecker.clear();

        printf("Level %u: ", lev); fflush(stdout);

        // Compute node hashes for this node and its children up to the level above the leaves
        unsigned int currentMatchDepth = std::max(_levels - lev - 2, 1u);
        buildMultiMap(currentMatchDepth, matchMaps, lev, lev + currentMatchDepth + 1);

        //////// FINDING EDGES FOR CLUSTERING ////////
        // For all nodes in this level, in reverse order (starting with least referenced)
        printf("\n\tFinding similar nodes [%i threads]... ", (int)omp_get_max_threads()); fflush(stdout);

        // Todo: Try out #pragma omp declare reduction (merge : std::vector<int> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp parallel for schedule(dynamic,2)
        for (int idA = (int) _data[lev].size() - 1; idA >= 0; --idA) {

			if ((idA % stepLogger == 0)) {
				printf("%.0f%%..", round(100.f * ((float) idA / (float)_data[lev].size())));
				fflush(stdout);
			}
            // only cluster nodes with 1 ref
            if (refCounts[lev][idA] > 1) continue; // should break but not sure if that affects omp stuff

            const Node &n = _data[lev][idA];
            // Todo: Should use pre-computed hash to avoid cascade of lossy error
            uint64_t nAKey = computeNodeHash(n, currentMatchDepth);

            // Todo: don't merge nodes that are a parent of a cluster representative

            // Break when ref count is greater than 1
            if (refCounts[lev][idA] == 1) {

                if (lev == _levels - 2) {
                    for (int i = 0; i < 64; ++i) {
                        uint64_t nAKeyMod = nAKey ^(1UL << i);
                        auto candidates = matchMaps[lev].equal_range(nAKeyMod);
                        for (auto it = candidates.first; it != candidates.second; ++it) {
                            id_t idB = it->second;
                            if (idA <= idB) continue; // don't match with itself or with previous nodes (since idB will already have been compared to idA)
                            if (refCounts[lev][idB] != 1) continue; // Only compare to other 1 ref nodes

                            // To avoid two identical edges, only add edges from low id to high id
                            // (No need to compare idB to idA and adding another edge, as it would be the same edge)
#pragma omp critical
                            {
                                edges.emplace_back(cluster::Edge{(unsigned int) idA, idB, 1});
                                uniqueNodesChecker[idA] = false;
                                uniqueNodesChecker[idB] = false;
                            };
                        }
                    }
                } else {
                    auto candidates = matchMaps[lev].equal_range(nAKey);

                    for (auto it = candidates.first; it != candidates.second; ++it) {
                        id_t idB = it->second;
                        if (idA <= idB) continue; // don't match with itself or with previous nodes (since idB will already have been compared to idA)
                        if (refCounts[lev][idB] != 1) continue; // Only compare to other 1 ref nodes

                        // Todo: Check how similar this node actually is: Deep comparison (expensive!!)
                        const Node &nB = _data[lev][idB];
                        unsigned int diff = 0;
                        this->diffSubtrees(lev, lev, n, nB, lossyDiff + 1, diff);

#pragma omp critical
                        if (diff <= lossyDiff) {
                            // Weights are similarity, so the inverse of the difference. Difference = dissimilarity
                            // Todo: Difference is currently linearly converted to similarity. Could also try 1 / diff
                            float sim = 1.0f - ((float) diff / ((float) lossyDiff + 1.0f));
                            edges.emplace_back(cluster::Edge{(unsigned int) idA, idB, sim});
                            uniqueNodesChecker[idA] = false;
                            uniqueNodesChecker[idB] = false;
                        }
                    }
                }
            }
        }

        // Clustering stage:
        // - Any node might have one or more potential matches, however, we can only pick one
        // - We prefer to match nodes with a node that is referenced more than once
        // - Then, find which nodes are potential matches the most frequently
        printf("\n\tClustering... "); fflush(stdout);
        const std::vector<std::vector<unsigned int>> clusters = cluster::clusterSubgraphs(edges, lev);
        printf(" Done! \n"); fflush(stdout);

        //////// COMPARING CLUSTERS TO OTHER NODES ////////
        for (id_t idA = 0; idA < _data[lev].size(); ++idA) {
            // Add all nodes with > 1 ref to uniqueNodes
            if (refCounts[lev][idA] == 1 && uniqueNodesChecker.find(idA) != uniqueNodesChecker.end()) continue;
            correspondences[idA] = (id_t) uniqueNodes.size(); // the correspondence is this node itself
            uniqueNodes.push_back(_data[lev][idA]);
        }

        // Replacing nodes from clusters
        for (unsigned int c = 0; c < clusters.size(); ++c) {
            // Representative node is put at index 0
            id_t repId = clusters[c][0];
            id_t repIdNew = uniqueNodes.size();
            /*
            const Node &n = _data[lev][repId];
            uint64_t repKey = computeNodeHash(n, currentMatchDepth);

            // Check if there is an node with > 1 ref count similar to the representative
            auto candidates = matchMaps[lev].equal_range(repKey);
            for (auto it = candidates.first; it != candidates.second; ++it) {
                id_t idB = it->second;
                if (repId == idB) continue; // don't match with itself
                if (refCounts[lev][idB] == 1) continue; // Only compare to nodes with more than 1 references

                // Check how similar this node actually is
                const Node &nB = _data[lev][idB];
                unsigned int diff = 0;
                this->diffSubtrees(lev, lev, n, nB, lossyDiff + 1, diff);

                if (diff <= lossyDiff) { // Todo: Allowed loss should be lower here
                    // Weights are similarity, so the inverse of the difference. Difference = dissimilarity
                    repId = idB;
                    repIdNew = correspondences[idB];
                    // Todo: Find best match instead of first
                    break;
                }
            }
            */

            // When representative is not replaced, add it as a unique node
            if (repIdNew == uniqueNodes.size()) {
                uniqueNodes.push_back(_data[lev][repId]);
                correspondences[repId] = repIdNew;
            }

            // The correspondence of all nodes in this cluster is the representative node
            for (const id_t corId : clusters[c]) {
                correspondences[corId] = repIdNew;
            }
        }

        /////////////////////////////////
        //// Replace previous level data
        _data[lev].clear();
        _data[lev].shrink_to_fit();
        uniqueNodes.shrink_to_fit();
        _data[lev] = uniqueNodes; // Replace all SVO nodes with the unique DAG nodes in this level
        _data[lev].shrink_to_fit();
        _nNodes += _data[lev].size();

        /////////////////////////////////
        //// Update all pointers in the level above
        for (id_t i = 0; i < _data[lev-1].size(); i++) {
            Node * bn = &_data[lev-1][i];
            // For all children...
            for (int j = 0; j < 8; j++) {
                // If this child exists...
                if (bn->existsChild(j)) {
                    // Set the child pointer to the unique node that replaced this child
                    bn->children[j] = correspondences[bn->children[j]];
                }
            }
        }
    }

    _stats.toDAGTime = _clock.now() - ts;

    _stats.nNodesDAG = _nNodes;

    printf("OK! [%s]\n", sl::human_readable_duration(_stats.toDAGTime).c_str());

    this->sortByEffectiveRefCount(); // only for getting new ref counts

    this->toDAG(false); // filter out new identical nodes (small probability, but might as well check)

}


////////////////////////////////////////////////////////////////////
////////////////////// CROSS LEVEL MERGING /////////////////////////
////////////////////////////////////////////////////////////////////


/**
 * @brief mergeAcrossAllLevels Brute force search over all levels, looking for equal subtrees and merging them
 * Goal: Research to see what the benefit of multi-level merging would be.
 * Just for SVDAGs now, but SSVDAGs would likely perform better - harder to check though
 */
unsigned int GeomOctree::mergeAcrossAllLevels() {
    /*
    std::function<std::string (const GeomOctree::Node &node, unsigned int depth)> computeNodeString;
    computeNodeString = [&](const GeomOctree::Node &node, unsigned int depth) {
        std::string key = "";
        // At depth 0, just return the child mask
        if (depth == 0) {
            key += std::to_string(node.childrenBitmask);
            return key;
        }
            // At depth 1, concat all 8-bit child masks into a 64 bit int
        else if (depth == 1) {
            for (unsigned int c = 0; c < 8; ++c) {
                // shift the bit mask of children into the key

                if (node.existsChildPointer(c)) {
                    key += std::to_string(_data[node.childLevels[c]][node.children[c]].childrenBitmask);
                } else {
                    key += "x";
                }
            }
        }
            // For higher depths, combine all children keys with XOR
        else {
            for (unsigned int c = 0; c < 8; ++c) {
                if (node.existsChildPointer(c)) {
                    const auto &child = _data[node.childLevels[c]][node.children[c]];
                    std::string childKey = computeNodeString(child, depth - 1);
                    key += "-" + childKey;
                } else {
                    key += "-x";
                }
            }
        }
        return key;
    };
     */

    /*
    // Prints all childmasks of two subtrees
	std::function<int (unsigned int levA, unsigned int levB,
		const GeomOctree::Node &nA, const GeomOctree::Node &nB)> getMasks;
	getMasks = [&](unsigned int levA, unsigned int levB, const GeomOctree::Node &nA, const GeomOctree::Node &nB) {
		printf("%u %u\n", nA.childrenBitmask, nB.childrenBitmask);
		if (levA == _levels - 1) {
			return 0;
		}
		for (int c = 0; c < 8; ++c) {
			printf("C%u: ", c);
			if (nA.existsChild(c)) {
				getMasks(levA + 1, levB + 1, _data[levA + 1][nA.children[c]], _data[levB + 1][nB.children[c]]);
			} else if (nB.existsChild(c)) {
				printf("X %u !!! !!! !!!\n", _data[levB + 1][nB.children[c]].childrenBitmask);
			} else {
				printf("X X\n");
			}
		}
		return 0;
	};
    */

    ///////////////////////////////////////////////////////////////
    /// Building multi-maps for finding potential matches faster //
    ///////////////////////////////////////////////////////////////
    // try out unordered map for performance improvements --->>> nothing changed
    std::vector<std::multimap<uint64_t, id_t>> matchMaps(_levels);

    // Initial depth should depend on total # of levels, 1 seems enough for ~8K, 2 or higher for more
    unsigned int currentMatchDepth = _levels / 2;
    buildMultiMap(currentMatchDepth, matchMaps);

    /////////////////////////////////
    /// Finding identical subtrees //
    /////////////////////////////////

    // Comparing ALL combinations of subtree pairs is not what we want
    // If a subtree high-up is equal to another subtree, all of its sub-subtrees will equal it as well
    // So, only keep track of children that are NOT already equal to another subtree
    // Only compare those instead of all nodes in the next level
    // Main algorithm:
    // For every node to currently check...
    // Check for levels higher-up whether there is an identical subtree
    // Note: If node B terminates before node A, they are NOT equal. This will not happen since we only compare to subtrees higher-up.
    // For nodes that are not equal, add their children in the next level to the nextNodesToCheck list.
    std::set<id_t> currentNodesToCheck;
    std::set<id_t> nextNodesToCheck;

    // Initialize all nodes of the highest level - 1 for comparison to subtrees at higher levels
    // This way, the matches are found as early as possible
    unsigned int levA = 1;
    for (id_t i = 0; i < _data[levA].size(); i++) {
        currentNodesToCheck.insert(i);
    }

    // Some variables for analytics
    unsigned int numTotalComparisons = 0;
    size_t prevNNodes = _nNodes;
    _nNodes = 1;

    // We need to store the correspondences of one subtree to another:
    // Contains for each level, a map of node IDs that point to level and index of an identical subtree higher-up in the graph.
    std::vector<std::map<id_t, std::pair<unsigned int, id_t>>> multiLevelCorrespondences(_levels); // store correspondences across different levels

    // Store all nodes that are visited in a subtree, so they can be potentially removed if a correspondence is found
    std::vector<std::map<id_t, std::pair<unsigned int, id_t>>> nodesInCurSubtree(_levels);

    // The level of each child pointer is set in node.childLevels[k]
    // then at the end, replace the node data in each level in those subtrees

    // The larger subtree is always found first, since duplicate subtrees are found top-down

    // Algorithm for removing duplicates so that indices don't get messed up:
    // - Idea: Same way as toDAG: Keep track of correspondences, replace data of whole level at one time bottom-up

    printf("Starting multi-level subtree comparisons:\n");

    int totHashCollisions = 0;

    for (; levA < _levels; ++levA) {
        _clock.restart();
        printf(" - L%u (%zu / %zu to check)... \t", levA, currentNodesToCheck.size(), _data[levA].size());
        fflush(stdout);

        // Build new match maps for the lowest levels with lower depths, when those nodes do not have subtrees of that depth
        unsigned int maxMatchDepth = _levels - levA - 1;
        if (maxMatchDepth < currentMatchDepth) {
            currentMatchDepth = maxMatchDepth; // currentMatchDepth / 2;
            buildMultiMap(currentMatchDepth, matchMaps);
        }

        int stepLogger = (int)round((currentNodesToCheck.size() + 1) / 10.f);

        unsigned long totalMatchCount = 0;

        // For all nodes to be checked...
        unsigned int curNodeIndex = 0;
        for (const auto& idA : currentNodesToCheck) {
            if ((curNodeIndex % stepLogger == 0)) {
                printf("%.0f%%..", round(100.f * (curNodeIndex / (float)currentNodesToCheck.size())));
                fflush(stdout);
            }
            curNodeIndex++;

            Node &nA = _data[levA][idA];
            uint64_t nAKey = computeNodeHash(nA, currentMatchDepth);

            bool foundMatch = false;

            int hashCollisions = 0;

            // also loop over levA itself, to find matches in the same level? No, that is already done in DAG
            for (unsigned int levB = 0; levB < levA; ++levB) {
                // idea: instead of checking in 'chronological' order, start looking at commonly chosen subtrees first
                // would probably not be worth it, since only a small portion of all nodes gets matches. Current optimization (hash) seems good enough

                // Instead of looping over EVERY node, just loop over potential matches!
                auto matchResult = matchMaps[levB].equal_range(nAKey);
                totalMatchCount += matchMaps[levB].count(nAKey);
                for (auto it = matchResult.first; it != matchResult.second; ++it) {
                    id_t j = it->second;
                    numTotalComparisons++;
                    Node &nB = _data[levB][j];

                    for (int i = 0; i < _levels; ++i) {
                        nodesInCurSubtree[i].clear();
                    }

                    // Brute force check equality of subtrees of nA and nB
                    bool areSubtreesEqual = compareSubtrees(levA, levB, nA, nB, nodesInCurSubtree);

                    if (areSubtreesEqual) {
                        // printf("\nMatch! LA %u LB %u, IDA: %u IDB %u\n", levA, levB, idA, j);
                        // getMasks(levA, levB, nA, nB);
                        foundMatch = true;
                        // Store that the subtree under the root of nodeA is identical to the subtree under nodeB
                        multiLevelCorrespondences[levA][idA] = std::make_pair(levB, j);

                        // Append all nodes in this subtree to the nodes that can be removed
                        for (int levRem = 0; levRem < _levels; ++levRem) { // levels below the root
                            multiLevelCorrespondences[levRem].insert(nodesInCurSubtree[levRem].begin(), nodesInCurSubtree[levRem].end());
                        }
                        //printf("-OK MATCH-");
                        break;
                    } else if (maxMatchDepth == currentMatchDepth) {
//                        std::cout << std::endl << computeNodeString(nA, currentMatchDepth) << std::endl;
//                        std::cout << computeNodeString(nB, currentMatchDepth) << std::endl;
////                        getMasks(levA, levB, nA, nB);
//                        printf("Matches correspond at same depth level, but no equal?!?!\n");
                        hashCollisions++;
                    }
                }
                if (foundMatch) {
                    break; // If a match is found, we are done for this subtree. It is equal to a subtree higher up in the graph!
                }
            }
            if (hashCollisions > 0) {
//                printf("Hash collisions: %i\n", hashCollisions);
                totHashCollisions += hashCollisions;
            }
            if (!foundMatch) {
                // If no match was found, try to find duplicate subtrees for all child nodes in the next iteration
                for (int i = 0; i < 8; ++i) {
                    if (nA.existsChild(i)) {
                        nextNodesToCheck.insert(nA.children[i]);
                    }
                }
            }
        }

        printf("Avg matches found per node: %.0f\n", totalMatchCount / float(currentNodesToCheck.size()));

        auto time = _clock.elapsed();

        printf(" -> %lu (%.0f%%) [%s]\n", _data[levA].size() - multiLevelCorrespondences[levA].size(), 100 * multiLevelCorrespondences[levA].size() / (float) _data[levA].size(), sl::human_readable_duration(time).c_str());

        // After all nodes for this level have been checked, swap the current and next nodes to check
        currentNodesToCheck.clear();
        currentNodesToCheck.swap(nextNodesToCheck);
    }


    printf("Total hash collisions: %i\n", totHashCollisions);

    ///////////////////////////////////
    /// Removing identical subtrees ///
    ///////////////////////////////////

    // Now that all identical subtrees have been identified, the duplicate subtrees can be removed and the pointers to them should be updated.
    for (unsigned int lev = _levels - 1; lev > 0; --lev) {
        std::vector<Node> uniqueNodes;
        uniqueNodes.reserve(_data[lev].size() - multiLevelCorrespondences[lev].size());

        std::vector<id_t> correspondences(_data[lev].size(), -1); // normal correspondences for this level (old index -> new index) for those that are not removed

        for (id_t nodeIndex = 0; nodeIndex < _data[lev].size(); ++nodeIndex) {
            // insert all nodes in uniqueNodes that do not have a correspondence in a higher level
            if (multiLevelCorrespondences[lev].count(nodeIndex) == 0) {
                correspondences[nodeIndex] = uniqueNodes.size();
                uniqueNodes.push_back(_data[lev][nodeIndex]);
            }
        }

        // Replace node data for this level
        _data[lev].clear();
        _data[lev].shrink_to_fit();
        _data[lev] = uniqueNodes;
        _data[lev].shrink_to_fit();

        _nNodes += _data[lev].size();

        int numReplaced = 0;
        int numRemained = 0;

        // Update all pointers in the level above
        for (id_t nodeIndex = 0; nodeIndex < _data[lev-1].size(); ++nodeIndex) {
            Node *node = &_data[lev - 1][nodeIndex];
            // For all children...
            for (int j = 0; j < 8; j++) {
                // If this child exists...
                if (node->existsChild(j)) {
                    // If it was replaced by a subtree higher up
                    auto it = multiLevelCorrespondences[lev].find(node->children[j]);
                    if (it != multiLevelCorrespondences[lev].end()) {
                        node->childLevels[j] = it->second.first;
                        node->children[j] = it->second.second;
                        // Node order in a higher level will change in future iteration...
                        // Therefore, the next loop updates pointers of nodes in lower level that point to nodes in this level
                        numReplaced++;
                    } else {
                        // Else, update the index from the normal list of correspondences
                        node->children[j] = correspondences[node->children[j]];
                        numRemained++;
                    }
                }
            }
        }


        // Update pointers from lower levels to nodes in this level
        for (unsigned int levLow = _levels - 2; levLow >= lev; --levLow) {
            for (id_t nodeIndex = 0; nodeIndex < _data[levLow].size(); ++nodeIndex) {
                Node *node = &_data[levLow][nodeIndex];
                // For all children...
                for (int j = 0; j < 8; j++) {
                    // If this node from a lower level points a node in the current level...
                    if (node->existsChild(j) && node->childLevels[j] == lev) {
                        if (correspondences[node->children[j]] == (id_t) -1) {
                            printf("\t\t- Missing correspondence on lev %u: Node %u, child %u\n", levLow, nodeIndex, j);
                        }

                        // Update where that pointer has been moved to
                        node->children[j] = correspondences[node->children[j]];
                    }
                }
            }
        }

        printf(" - L %i Replaced/remained: %i / %i\n", lev, numReplaced, numRemained);
    }

    _stats.nNodesDAG = _nNodes;


    /////////////////////////////////
    /// Done: Print the results!   //
    /////////////////////////////////
    printf("Total #nodes %zu / %zu. Total #comparisons: %u\n", _nNodes, prevNNodes, numTotalComparisons);

    int totalElimNodes = 0;
    for (unsigned int i = 0; i < _levels; ++i) {
        totalElimNodes += multiLevelCorrespondences[i].size();
        id_t origSize = _data[i].size() + multiLevelCorrespondences[i].size();
        double pct = 100 * (multiLevelCorrespondences[i].size() / double(origSize));
        printf(" - Level %u:   \t%zu subtrees are equal to a subtree higher up. %zu / %i (%.2f%%) nodes of this level have been removed\n", i, multiLevelCorrespondences[i].size(), multiLevelCorrespondences[i].size(), origSize, pct);
    }
    printf("Total number of nodes that was removed: %u / %zu (%.2f%%)\n ", totalElimNodes, prevNNodes, (100 * totalElimNodes / double(prevNNodes)));


//    printf("Indirect subtree feasibility: How many unique pointers there are to other levels, per level\n");
//    Todo: Should use this to check how limiting it is to only have pointers to 1 or 2 levels instead of all
//    for (unsigned int lev = 0; lev < _levels; ++lev) {
//        std::vector<std::set<id_t>> uniquePointers(_levels);
//        std::vector<unsigned int> numPointers(_levels);
//        for (const auto &node : _data[lev]) {
//            for (int c = 0; c < 8; ++c) {
//                if (node.childLevels[c] != lev + 1) {
//                    uniquePointers[node.childLevels[c]].insert(node.children[c]);
//                    numPointers[node.childLevels[c]]++;
//                }
//            }
//        }
//
//        unsigned int totalNumUniquePointers = 0;
//        unsigned int totalNumPointers = 0;
//        for (unsigned int lev2 = 0; lev2 <= lev; ++lev2) {
//            printf("    L%u -> L%u: Unique / total = %zu / %u\n", lev, lev2, uniquePointers[lev2].size(), numPointers[lev2]);
//            totalNumUniquePointers += uniquePointers[lev2].size();
//            totalNumPointers += numPointers[lev2];
//        }
//        printf("  L%u total: Unique / total = %u / %u\n", lev, totalNumUniquePointers, totalNumPointers);
//    }


    return totalElimNodes;
}



void GeomOctree::symMergeAcrossAllLevels() {
    // Needs to be in SDAG state

    if (_state == S_SDAG) {
        printf("* Transforming SDAG -> CSDAG ... "); fflush(stdout);
    } else {
        printf("ERROR! This is not a SDAG!\n");
        return;
    }

    // Build multi-map of MirrorNodes
    struct MirroredNode{
        MirroredNode() : mirrorX(false), mirrorY(false), mirrorZ(false), id(0) {}
        MirroredNode(bool mX, bool mY, bool mZ, id_t id) :
                mirrorX(mX), mirrorY(mY), mirrorZ(mZ), id(id) {}
        bool mirrorX, mirrorY, mirrorZ;
        id_t id;
    };

    // For each of the 8 symmetry options, for each level, a multi-map of hash -> index for each node
    // None, X, Y, Z, XY, XZ, YZ, XYZ
    std::vector<std::vector<std::multimap<uint64_t, id_t, uint8_t>>> matchMaps(9);
    for (int i = 0; i < 8; ++i) {
        matchMaps[i].resize(_levels);
    }

    auto buildMultiMap = [&](unsigned int depth) {
        printf("[Building match maps @ depth %u... ", depth);
        fflush(stdout);
        for (unsigned int lev = 0; lev < _levels - depth; ++lev) {
            matchMaps[lev].clear();
            for (id_t nodeIndex = 0; nodeIndex < _data[lev].size(); ++nodeIndex) {

                // Todo: Insert all symmetry options in matchMaps
                for (int symIndex = 0; symIndex < 8; ++symIndex) {
                    // Todo: Check invariance
                    uint64_t key = 0; // computeNodeKey(_data[lev][nodeIndex], depth);
//                    matchMaps[symIndex][lev].insert(std::make_pair(key, nodeIndex));
                }
            }
        }
        printf("Done!]");
    };


    unsigned int currentMatchDepth = _levels / 2;

    // Store correspondences across different levels: For every level, a map of a node and its match (level, index, symmetry index)
    std::vector<std::map<id_t, std::pair<unsigned int, std::pair<id_t, int>>>> multiLevelCorrespondences(_levels);
    // Store all nodes that are visited in a subtree, so they can be removed when a correspondence is found
    std::vector<std::map<id_t, std::pair<unsigned int, id_t>>> nodesInCurSubtree(_levels);

    // Todo: Use structs instead of pairs in pairs in pairs

    //
    unsigned int levA = 1;
    for (; levA < _levels; ++levA) {

        unsigned int maxMatchDepth = _levels - levA - 1;
        if (maxMatchDepth < currentMatchDepth) {
            currentMatchDepth = maxMatchDepth;
            buildMultiMap(currentMatchDepth);
        }

        for (unsigned int idA = 0; idA < _data[levA].size(); ++idA) {
            Node &nA = _data[levA][idA];
            // Todo: Check if already matched through a parent node

            bool foundMatch = false;

            uint64_t nAKey = 0; // computeNodeKey(nA, currentMatchDepth);
            for (unsigned int levB = 0; levB < levA; ++levB) {
                for (int symIndex = 0; symIndex < 8; ++symIndex) {
//                    auto matchResult = matchMaps[symIndex][levB].equal_range(nAKey);
//                    for (auto it = matchResult.first; it != matchResult.second; ++it) {
//                        id_t idB = it->second;
//                        bool areSubtreesEqual = false; // compareSymSubtrees(levA, levB, nA, nB, nodesInCurSubtree);
//                        if (areSubtreesEqual) {
//                            foundMatch = true;
//                            multiLevelCorrespondences[levA][idA] = std::make_pair(levB, std::make_pair(idB, symIndex));
//                            // Append all nodes in this subtree to the nodes that can be removed
//                            for (int levRem = 0; levRem < _levels; ++levRem) { // levels below the root
//                                // Todo: take into account mirror bits of child pointers
////                                multiLevelCorrespondences[levRem].insert(nodesInCurSubtree[levRem].begin(), nodesInCurSubtree[levRem].end());
//                            }
//                            break;
//                        }
//                    }
//                    if (foundMatch) break;
                }
                if (foundMatch) break;
            }
        }
    }
}



////////////////////////////////////////////////////////////////////
////////////////// EXPLOITING HIDDEN GEOMETRY //////////////////////
////////////////////////////////////////////////////////////////////


void GeomOctree::hiddenGeometryFloodfill() {
    printf("Starting flood fill...\n");
    if (_state != S_SVO) {
        printf("Not in SVO format, aborting!\n");
        return;
    }

    ///// Idea //////
    // From the most NXNYNZ node, start a flood fill:
    // For every one of its neighbours, if it hasn't been visited yet, check if it intersects with geometry
    // If it doesn't intersect, then continue the floodfill for its neighbours

    // Todo: This doesn't work if there is blocking geometry, e.g. a box across the whole XY axis

    ////// Helper structs/functions ////////

    // Contains all info to find a node and a reference to its parent
    struct TravNode {
        TravNode(std::shared_ptr<TravNode> prnt, id_t pIdx, ChildrenIdx i, unsigned int l)
                : parent(std::move(std::move(prnt))), parentIdx(pIdx), childIdx(i), level(l) {}
        std::shared_ptr<TravNode> parent;
        id_t parentIdx;
        ChildrenIdx childIdx;
        unsigned int level;

        bool operator<(const TravNode &other) const {
            return level < other.level
                   || parentIdx < other.parentIdx
                   || childIdx < other.childIdx;
        }
    };

    TravNode nullTravNode(nullptr, 0, NXNYNZ, -1);

    auto getNode = [&](const TravNode &tn) {
        if (tn.level == 0) return _data[0][0]; // Root node
//        assert((tn.level < 0 || tn.level >= _levels));
        Node& p = _data[tn.level - 1][tn.parentIdx];
        return (Node&) _data[tn.level][p.children[tn.childIdx]];
    };

    // A leaf is not represented with a node, it is a child of a node with a child bitmask of 0 or at max level (_levels)
    auto isLeaf = [&](const TravNode &tn) {
        if (tn.level == _levels) return true;
        // if tn points to an actual node, and that node has no children, then it's a leaf
        // if there is a child pointer in the parent, it's a leaf it it has no children. Else, it's a leaf
        if (getNode(*(tn.parent)).existsChildPointer(tn.childIdx)) return !getNode(tn).hasChildren();
        return true;
    };

    // Child indices with same dir (e.g. NX -> NXNYNZ, NXPYNZ, NXNYPZ, NXPYPZ)
    // Returns one of the 4 child indices in the given direction of a node (i between 0 < 4)
    auto getChildIndexForDirection = [](DirectionIdx dir, unsigned int i) {
        int     /* NX */    childIdx = i; // NX
        if      (dir % 3 == NY) childIdx = i + 2 * (i / 2);
        else if (dir % 3 == NZ) childIdx = i * 2;
        if      (dir == PX) childIdx += 4;
        else if (dir == PY) childIdx += 2;
        else if (dir == PZ) childIdx += 1;
        return childIdx;
    };


    // Neighbour finding algorithm: https://geidav.wordpress.com/2017/12/02/advanced-octrees-4-finding-neighbor-nodes/

    std::function<TravNode (const TravNode &tn, DirectionIdx dir)> getNeighbourGrtrEqSz;
    getNeighbourGrtrEqSz = [&](const TravNode &tn, const DirectionIdx dir) {
        if (tn.level == 0) return nullTravNode; // Root node, cannot go further up

        // For each of the 4 children on a side
        for (int i = 0; i < 4; ++i) {
            // If tn is a child at positive side, return a neighbour within the same parent at the opposite side, vice versa
            auto oppositeSideChildIdx = (ChildrenIdx) getChildIndexForDirection((DirectionIdx) ((dir + 3) % 6), i);
            if (tn.childIdx == oppositeSideChildIdx) {
                auto curSideChildIdx = (ChildrenIdx) getChildIndexForDirection((DirectionIdx) dir, i);
                return TravNode(tn.parent, tn.parentIdx, curSideChildIdx, tn.level);
            }
        }

        // Else, try to find the neighbour of the parent in this direction
        TravNode pNb = getNeighbourGrtrEqSz(*(tn.parent), dir);
        // Return it if it's the nullNode, or it's a leaf
        if (pNb.parent == nullptr || isLeaf(pNb)) return pNb;

        auto pNbPtr = std::make_shared<TravNode>(pNb);
        id_t pNbIdx = _data[pNb.level - 1][pNb.parentIdx].children[pNb.childIdx];
        // If the neighbour of the parent is not a leaf node, return its child that is the closest to the given node
        // For each of the 4 children on this side of the node...
        for (int i = 0; i < 4; ++i) {
            auto curSideChildIdx = (ChildrenIdx) getChildIndexForDirection((DirectionIdx) dir, i);
            if (tn.childIdx == curSideChildIdx) {
                auto oppositeSideChildIdx = (ChildrenIdx) getChildIndexForDirection((DirectionIdx) ((dir + 3) % 6), i);
                return TravNode(pNbPtr, pNbIdx, oppositeSideChildIdx, pNb.level + 1);
            }
        }

        return nullTravNode;
    };

    static std::stack<TravNode> candidates;
    auto getNeighboursSmSz = [&](const TravNode &tn, const TravNode &tnNb, const DirectionIdx dir, std::vector<TravNode> &neighbours) {
        if (tnNb.parent != nullptr) candidates.push(tnNb);

        while (!candidates.empty()) {
            const TravNode &can = candidates.top(); candidates.pop();

            if (isLeaf(can)) {
                // if it's a leaf, add as a neighbour
                neighbours.emplace_back(can);
            } else {
                // else, add children in on the opposite side of the given direction as candidates
                auto canPtr = std::make_shared<TravNode>(can);
                id_t canIdx = getNode(*(can.parent)).children[can.childIdx];
                unsigned int newLev = can.level + 1;

                // For each of the 4 children on a side of the node, add candidates
                for (int i = 0; i < 4; ++i) {
                    auto oppositeSideChildIdx = (ChildrenIdx) getChildIndexForDirection((DirectionIdx) ((dir + 3) % 6), i);
                    candidates.push(TravNode(canPtr, canIdx, oppositeSideChildIdx, newLev));
                }
            }
        }
    };

    auto getNeighbours = [&](const TravNode &tn, std::vector<TravNode> &neighbours) {
        // Get neighbours for all directions
        for (int dir = 0; dir <= PZ; ++dir) {
            // Get the neighbour directly next to it in this direction
            const TravNode &nb = getNeighbourGrtrEqSz(tn, (DirectionIdx) dir);
            // find the leaf nodes on the side of the input node inside the neighbour
            getNeighboursSmSz(tn, nb, (DirectionIdx) dir, neighbours);
        }
    };

    ///// Finding node to start the floodfill /////
    TravNode rootTrav(NULL, 0, NXNYNZ, 0);
    TravNode startTrav = nullTravNode;
    std::shared_ptr<TravNode> prevParent = std::make_shared<TravNode>(rootTrav);

    // Find the deepest node at the most negative corner of the scene
    ChildrenIdx startIdx = PXPYPZ; // Start flood fill at the top right of the scene
    id_t curNodeIdx = 0;
    for (unsigned int lev = 0; lev < _levels; ++lev) {
        Node &node = _data[lev][curNodeIdx];

        // Create new TravNode containing parent info
        TravNode curNode(prevParent, curNodeIdx, startIdx, lev + 1);
        prevParent = std::make_shared<TravNode>(curNode);
        curNodeIdx = _data[lev][curNodeIdx].children[startIdx];
        // Replace the start node with its child
        startTrav = curNode;


        // Stop if no child is to be found
        if (!node.existsChildPointer(startIdx)) {
            break;
        }
    }

    // Check if the starting node intersects with geometry: Then the initial node is inside, so not a good starting point
    if (_data[startTrav.level - 1][startTrav.parentIdx].existsChild(startTrav.childIdx)) {
        printf("Deepest corner node intersects with geometry, not dealing with this now... Aborting flood fill\n");
        exit(1);
    }

    // Flood fill, using a queue of nodes that need to be visited
    static std::vector<TravNode> neighbours;

    static std::stack<TravNode> queue;
    queue.push(startTrav);
    _data[0][0].outsideMask = 0;
    _data[startTrav.level - 1][startTrav.parentIdx].setChildOutsideBit(startTrav.childIdx);

    unsigned int numOutsideLeafs = 0;

    while (!queue.empty()) {
        const TravNode travNode = queue.top(); queue.pop();

        // For all neighbours
        neighbours.clear();
        getNeighbours(travNode, neighbours);
        for (const auto &nbTravNode : neighbours) {

            const Node &nbP = getNode(*(nbTravNode.parent)); // neighbour parent
            // Node should be checked if 1. not intersects with geometry and 2. not already marked as outside
            if (!nbP.getChildOutsideBit(nbTravNode.childIdx)) {
                // Any neighbour of an outside node is also outside, as nodes on a surface also count as outside
                _data[nbTravNode.level - 1][nbTravNode.parentIdx].setChildOutsideBit(nbTravNode.childIdx);
                if (!nbP.existsChild(nbTravNode.childIdx)) {
                    // Add to queue to visit next
                    queue.push(nbTravNode);
                    // mark as being outside of any geometry in parent
    //                nbP.setChildOutsideBit(nbTravNode.childIdx);
                    numOutsideLeafs++;
                }
            }
        }
    }
    // Todo: Also consider what to happens with multiple build steps!

    printf("Propagating inside/outside node status up the graph...\n");
    // After all of the deepest nodes have been marked as inside/outside, propagate those values to their parents
    for (int lev = _levels - 2; lev >= 0; --lev) {
        int outsideNodes = 0;
        for (id_t nodeId = 0; nodeId < _data[lev].size(); ++nodeId) {
            Node &n = _data[lev][nodeId];
            for (int c = 0; c < 8; ++c) {
                if (n.existsChildPointer(c) && !_data[lev + 1][n.children[c]].isInside()) {
                    _data[lev][nodeId].setChildOutsideBit(c);
                    outsideNodes++;
                }
            }
        }
        printf("Level %i: %.2f%% of SVO nodes outside\n", lev + 1, 100.0 * outsideNodes / (float) _data[lev + 1].size());
    }


    // Check to see if floodfill works: Set all inside nodes to bitmask with 1
#if 0
    printf("DEBUG: Setting inside nodes to first node of level...\n");
    for (int lev = _levels - 1; lev >= 0; --lev) {
        for (id_t nodeId = 0; nodeId < _data[lev].size(); ++nodeId) {
            Node &parent = _data[lev][nodeId];
            for (int c = 0; c < 8; ++c) {
                if (!parent.existsChild(c) && !parent.getChildOutsideBit(c)) {
                    parent.setChildBit(c);
                    parent.children[c] = 0;
                }
//                if (parent.existsChild(c) && parent.getChildOutsideBit(c)) {
//                    _data[lev][nodeId].unsetChildBit(c);
//                    _data[lev][nodeId].children[c] = 0;
//                }
            }
        }
    }

#endif


    printf("Done with flood fill! %u leaves outside of geometry\n", numOutsideLeafs);
}

void GeomOctree::toHiddenGeometryDAG() {

    // Same as normal DAG, but ignore the nodes that are inside of geometry
    // Also compare with most common referenced nodes first

    printf("* Performing hidden geometry exploitation... "); fflush(stdout);

    // Subtree comparison
    // Option 1.
    // - Precompute hashes for all 256 posibilities of the inside mask
    // - Computational speedup possible from current implementation: Perform bottom-up instead of from scratch
    // - Note: Needs to be computed before toDAG, so there will be a large amount of nodes...

    // Todo: For memory issues, try to merge to a DAG while also taking into account the outsideMask
    // High likelihood that identical nodes have an identical outsideMask, at least in lower levels

    printf(" * - Preparing hash maps:"); fflush(stdout);
    // For every combination of hidden child masks (256) > For every level > For every node
    std::vector<std::vector<std::multimap<uint64_t, id_t>>> hidChildMatchMaps;
    std::vector<std::vector<std::vector<uint64_t>>> hidChildHashes;
    for (int i = 0; i < 256; i++) {
        if (i % 16 == 0) {
            printf("%.2f%% ", 100.0 * i / 256.0); fflush(stdout);
        }
        hidChildMatchMaps.emplace_back(std::vector<std::multimap<uint64_t, id_t>>(_levels));
        hidChildHashes.emplace_back(std::vector<std::vector<uint64_t>>(_levels));

        // Compute all candidate match maps for all levels where the hidden child mask equals i
        buildMultiMapBotUp(hidChildMatchMaps[i], hidChildHashes[i], i);
    }
    printf("\n");

    ///////////////////////////////////////////////////////
    // Todo: Other option: Modify SVO to an optimal state, so that we can call toDAG as usual

    printf(" * - Converting to DAG:\n");
    // Now perform toDAG and check for matches 256 times
    _nNodes = 1;
    std::vector<id_t> correspondences;
    std::vector<Node> uniqueNodes;

    // For every level, starting at the leaves...
    for (unsigned int lev = _levels - 2; lev > 1; --lev) {
        // Clear the lists used to keep track of correspondences etc
        size_t oldLevSize = _data[lev].size();
        uniqueNodes.clear();
        uniqueNodes.shrink_to_fit();
        correspondences.clear();
        correspondences.resize(oldLevSize);

        bool foundMatch = false;

        // For all nodes in this level...
        // Todo: compare to most referenced first?
        for (id_t i = 0; i < _data[lev].size(); i++) {
            Node n = _data[lev][i];
            if (!n.hasChildren()) continue; // skip empty nodes
//            if (lev == 1) {
//                printf("Node %i has hash %u", i, hidChildHashes[i][lev][0])
//            }

            // Todo: Should apply clustering here, or at least find the best candidate
            // If this node is already used as a unique node, skip it
            if (std::find(correspondences.begin(), correspondences.end(), i) != correspondences.end()) {
                continue;
            }

            uint64_t hash = hidChildHashes[n.outsideMask][lev][i];

            for (int h = 0; h < 256; h++) {
                auto candidates = hidChildMatchMaps[h][lev].equal_range(hash);
                for (auto it = candidates.first; it != candidates.second; ++it) {
                    id_t matchId = it->second;
                    if (matchId != i) {
                        // Todo: Deep comparison?
                        foundMatch = true;
                        correspondences[i] = matchId; // store duplicate node
                        break;
                    }
                }
                if (foundMatch) break;
            }


            if (!foundMatch) { // !found
                correspondences[i] = (id_t)uniqueNodes.size(); // the correspondence is this node itself
                uniqueNodes.push_back(n);
            }
        }

//        if (!iternalCall)
            printf("Reduced level %u from %lu to %lu nodes\n", lev, _data[lev].size(), uniqueNodes.size());

        _data[lev].clear();
        _data[lev].shrink_to_fit();
        uniqueNodes.shrink_to_fit();
        _data[lev] = uniqueNodes; // Replace all SVO nodes with the unique DAG nodes in this level
        _data[lev].shrink_to_fit();
        _nNodes += _data[lev].size();

        // Update all pointers in the level above
        for (id_t i = 0; i < _data[lev-1].size(); i++) {
            Node * bn = &_data[lev-1][i];
            // For all children...
            for (int j = 0; j < 8; j++) {
                // If this child exists...
                if (bn->existsChild(j)) {
                    // Set the child pointer to the unique node that replaced this child
                    bn->children[j] = correspondences[bn->children[j]];
                }
            }
        }
    }

//    _stats.toDAGTime = _clock.now() - ts;

    _state = S_DAG;
    _stats.nNodesDAG = _nNodes;


    // It's been a while, rewriting pseudoce as a fresh view
    /**
     * Input: SVO with an outsideMask for each node, specifying whether the node is visible from the outside
     * Goal: Merging nodes that appear the same from the outside
     * - This is where the visible side of a node is identical to that same part of another node, and the invisible part is arbitrary
     * - So, for each node, look for a node that has the same outside mask, or where there are more parts on the outside
     *
     * Nodes are put into multimaps to find them quickly, based on their children and outsideMask
     * - We have a node A, and we look for another node B that we can merge it with
     * - We compute a hash for node A, and variations for each of the bits in the outsideMask that are inside
     * - We only need 1 multi-map as usual, which integrates the outsideMask for each hash, since we look for each of the variations by computing different hashs for A
     *
     * Seems like it should work
     * * For every node
     *   - If it is partially inside
     *     - Create hashes for all possible candidates (each hash of the node with one of the inside bits off)
     *     - Find best candidate
     * *
     *
     * For what is left over:
     *   - If it completely inside, remove all children
     *
     * Should be down top-down:
     * *
     *
     * SO. Either
     * 1) create a single hash for each node, taking into consideration their insideMasks, or:
     *      Then for candidate lookups, compute each hash possibility for the node in question
     * 2) Create several hashes, for each combination of the children that are inside as if they are outside
     *      Then for candidate lookups, compute one hash possibility
     *
     * In other words, for each node that has one or more children that are not exposed to the outside world,
     * we can compute what the hash of a node would be that could replace it, and we can quickly look up
     * whether such a node exists. Such a hash value can be computed for each of the children that is inside.
     *
     * Similar problem as lossy: If a node is marked as a correspondence, and itself can also correspond
     * to another node, which is the better option? Compute all possibilities again?
     *
     * Also, should the initial hash values take into account the insideMask?
     * Yes: Many nodes will not have a correspondence, so a lot of potential matches are lost, as the hidden geometry
     *      can match to outside geometry of other nodes. When taking into account the mask, it is effectively set to empty space
     * No:  Then it doesn't work at all. Can't generate hashes for nodes when take can have any possible geometry.
     * Conclusion: Creating the initial hash map for all options seems like the best.
     *      Essentially: This node can serve as any of these options, this node as any of these other options, etc.
     *
     * -> Create hashes bottom-up. For each of the combinations of x hidden children (2^x), create a hash
     * -> Loop over every node, for each of the combinations of hidden children, look up candidates
     *
     * Note: probably needs a pre-processing step (mentioned earlier) where identical nodes with identical outsideMask are merged
     *
     *
     */


    // Option 2.
    // - Use set intersections to find which nodes to compare to
}

