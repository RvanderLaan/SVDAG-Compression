//
// Created by remi on 19-09-19.
//

#pragma once

#include <vector>
#include <map>
#include <set>
#include <stack>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <omp.h>
#include <assert.h>

class cluster {

public:
    struct Edge {
        unsigned int sourceId;
        unsigned int targetId;
        float weight;

        bool operator < (const Edge &other) const {
            if (sourceId != other.sourceId) return sourceId < other.sourceId;
            else                            return targetId < other.targetId;
        }
        bool operator == (const Edge &other) const {
            return sourceId == other.sourceId && targetId == other.targetId;
        }
    };

    static inline std::vector<std::vector<unsigned int>> clusterSubgraphs(
            const std::vector<Edge> &edges, int id = 0) {
        const auto subGraphs = findSubGraphs(edges);
        printf("Found %zu subgraphs of %zu edges\n", subGraphs.size(), edges.size());
        std::vector<std::vector<unsigned int>> clusters;

        // Threshold of cluster size for when to run MCL or when to just create a cluster of the whole subgraph
        double graphSizeThresholdToCluster = sqrt(edges.size()) / 4.0;

        for (const auto& subGraph : subGraphs) {
            if (subGraph.size() <= graphSizeThresholdToCluster) {
                // For extremly small clusters, the subgraph is one cluster
                std::set<unsigned int> clusterSet;
                for (const auto& e : subGraph) {
                    clusterSet.emplace(e.sourceId);
                    clusterSet.emplace(e.targetId);
                }
                std::vector<unsigned int> clusterVec(clusterSet.begin(), clusterSet.end());

                // Find cluster representative, with highest total edge weight
                std::map<unsigned int, float> summedEdgeWeights;
                for (const auto &e : subGraph) {
                    if (summedEdgeWeights.count(e.sourceId) == 0) {
                        summedEdgeWeights[e.sourceId] = e.weight;
                    } else {
                        summedEdgeWeights[e.sourceId] += e.weight;
                    }
                    if (summedEdgeWeights.count(e.targetId) == 0) {
                        summedEdgeWeights[e.targetId] = e.weight;
                    } else {
                        summedEdgeWeights[e.targetId] += e.weight;
                    }
                }

                // Swap two elements so that node with highest weight is at the index 0
                auto max = std::max_element(summedEdgeWeights.begin(), summedEdgeWeights.end());
                auto maxIt = std::find(clusterVec.begin(), clusterVec.end(), max->first);
                unsigned int maxIndex = std::distance(clusterVec.begin(), maxIt);
                unsigned int tmp = clusterVec[0];
                clusterVec[0] = max->first;
                clusterVec[maxIndex] = tmp;

                clusters.emplace_back(clusterVec);
            } else {
                // Else, perform MCL
                auto subClusters = MCL(subGraph, id);
                clusters.insert(clusters.end(), subClusters.begin(), subClusters.end());

                // Todo: Set cluster representative at index 0
            }
        }
        printf("Total of %zu clusters\n", clusters.size());
        return clusters;
    }

    static inline std::vector<std::vector<unsigned int>> MCL(
            const std::vector<Edge> &edges, int id = 0) {
        std::vector<std::vector<unsigned int>> clusters;
        if (edges.empty()) return clusters;

        std::string fnIn = "edges" + std::to_string(id) + ".txt";
        std::string fnOut = "clusters" + std::to_string(id) + ".txt";

        std::ofstream edgesFile(fnIn, std::ios::out | std::ios::trunc);

        for (auto const& edge : edges) {
            edgesFile << std::to_string(edge.sourceId) << " "
                   << std::to_string(edge.targetId) << " "
                   << edge.weight << "\n"; stdout;
        }
        edgesFile.close();

        std::string nT = std::to_string((int) omp_get_max_threads());
        std::string cmd = "$HOME/local/bin/mcl " + fnIn + " --abc -I 2.0 -q x -V all -te " + nT + " -o " + fnOut; // -scheme 1
//        std::string cmd = "export OMP_NUM_THREADS=" + nT + "; ../../hipmcl/bin/hipmcl -M " + fnIn + " -I 2.0 -per-process-mem 0.1 -o " + fnOut;

        int res = std::system(cmd.c_str());

        if (res != 0) {
            printf("Could not perform MCL\n");
            exit(res);
        }

        std::ifstream clusterFile(fnOut);
        std::string line;
        while(std::getline(clusterFile, line)) {
            std::stringstream linestream(line);
            std::string substr;

            clusters.emplace_back(std::vector<unsigned int>());

            // mcl uses tab, hipmcl uses space
            const char delim = '\t'; // ' '

            while (std::getline(linestream, substr, delim)) {
                clusters[clusters.size() - 1].emplace_back(std::stoul(substr));
            }
        }
        clusterFile.close();

        // Todo: Remove files
#if 0
        std::remove(fnIn);
        std::remove(fnOut);
#endif

        return clusters;
    }

    static std::vector<std::vector<Edge>> findSubGraphs(const std::vector<Edge> &edges) {
        std::map<unsigned int, std::vector<Edge>> edgeMap;
        std::vector<std::vector<Edge>> subGraphs;

        // Create map of each node with all of its edges
        // Add each to both nodes since it's an undirected graph
        for (const Edge &e : edges) {
            if (edgeMap.count(e.sourceId) == 0) {
                edgeMap[e.sourceId] = std::vector<Edge>{e};
            } else {
                edgeMap[e.sourceId].emplace_back(e);
            }
            if (edgeMap.count(e.targetId) == 0) {
                edgeMap[e.targetId] = std::vector<Edge>{e};
            } else {
                edgeMap[e.targetId].emplace_back(e);
            }
        }

        // Flood fill through all edges to find subgraphs
        std::stack<unsigned int> queue; // Queue of nodes to visit
        while (!edgeMap.empty()) { // Until all edges have been visitied
            unsigned int n = edgeMap.begin()->first; // Start with the first item left in the edge map
            const auto nEdges = edgeMap[n];
            edgeMap.erase(n);

            std::set<Edge> g(nEdges.begin(), nEdges.end());

            for (const Edge &e : nEdges) {
                queue.emplace(n == e.sourceId ? e.targetId : e.sourceId);
            }

            while (!queue.empty()) {
                unsigned int n2 = queue.top(); queue.pop();
                const auto &n2EdgesIt = edgeMap.find(n2);

                if (n2EdgesIt != edgeMap.end()) {
                    const std::vector<Edge> n2Edges = n2EdgesIt->second;
                    edgeMap.erase(n2);

                    // Add edges to the subgraph
                    // g.insert(n2Edges.begin(), n2Edges.end()); // WHY DOES THIS NOT WORK C++?? THAT COST ME A COUPLE OF HOURS!!12jsdfoguasdfhgas
                    std::copy(n2Edges.begin(), n2Edges.end(), std::inserter(g, g.end()));

                    // Add all connected edges to n2 to the queue
                    for (const Edge &e : n2Edges) {
                        queue.emplace(e.sourceId ? e.targetId : e.sourceId);
                    }
                }
            }
            std::vector<Edge> gVec(g.begin(), g.end());
            subGraphs.emplace_back(gVec);
        }

        // int i = 0;
        // for (const Edge &e : edges) {
        //     bool found = 0;
        //     for (const auto &subGraph : subGraphs) {
        //         if (std::find(subGraph.begin(), subGraph.end(), e) != subGraph.end()) {
        //             found = true;
        //             break;
        //         }
        //     }
        //     if (!found) {
        //         printf("Edge %i not found (%u -> %u)\n", i, e.sourceId, e.targetId);
        //     }
        //     i++;
        // }

        // unsigned int numEdgesInSubgraphs = 0;
        // for (const auto &subGraph : subGraphs) {
        //     numEdgesInSubgraphs += subGraph.size();
        // }
        // assert(numEdgesInSubgraphs == edges.size()); // ensure the amount of edges is the same

        return subGraphs;
    }

    static void TestSubgraph() {
        std::vector<Edge> edges {
            cluster::Edge{ 0, 1, 1.0f },
            cluster::Edge{ 1, 2, 1.0f },
            cluster::Edge{ 0, 2, 1.0f },
            cluster::Edge{ 2, 3, 1.0f }
        };

        const auto subGraphs = cluster::findSubGraphs(edges);
        printf("%zu\n", subGraphs.size());
        assert(subGraphs.size() == 1);
    }

     static void TestSubgraph2(std::string edgeFn) {
         std::vector<Edge> edges;

        std::ifstream edgeFile(edgeFn);
        std::string line;
        while(std::getline(edgeFile, line)) {
            std::stringstream linestream(line);
            std::string substr;
            const char delim = ' ';

            Edge e;
            std::getline(linestream, substr, delim);
            e.sourceId = std::stoul(substr);
            std::getline(linestream, substr, delim);
            e.targetId = std::stoul(substr);

            edges.emplace_back(e);
            
        }
        edgeFile.close();

        const auto subGraphs = findSubGraphs(edges);
        printf("Found %zu subgraphs from %zu edges\n", subGraphs.size(), edges.size());
     }
};


