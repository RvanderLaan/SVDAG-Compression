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

class cluster {

public:
    struct Edge {
        unsigned int sourceId;
        unsigned int targetId;
        float weight;
    };

    static inline std::vector<std::vector<unsigned int>> clusterSubgraphs(
            const std::vector<Edge> &edges, int id = 0) {
        const auto subGraphs = findSubGraphs(edges);
        printf("Found %zu subgraphs of %zu edges\n", subGraphs.size(), edges.size());
        std::vector<std::vector<unsigned int>> clusters;

        for (const auto& subGraph : subGraphs) {
            if (subGraph.size() <= 9999) {
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
        for (const Edge &e : edges) {
            if (edgeMap.count(e.sourceId) == 0) {
                edgeMap[e.sourceId] = std::vector<Edge>{e};
            } else {
                edgeMap[e.sourceId].emplace_back(e);
            }
        }

        // Flood fill through all edges to find subgraphs
        std::stack<unsigned int> queue; // Queue of nodes to visit
        while (!edgeMap.empty()) { // Until all edges have been visitied
            unsigned int n = edgeMap.begin()->first;
            const auto nEdges = edgeMap[n];
            edgeMap.erase(n);

            std::vector<Edge> g = std::vector<Edge>();
            g.insert(g.end(), nEdges.begin(), nEdges.end());
            subGraphs.emplace_back(g);

            for (const Edge &e : nEdges) {
                queue.emplace(e.targetId);
            }

            while (!queue.empty()) {
                unsigned int n2 = queue.top(); queue.pop();
                const auto &n2EdgesIt = edgeMap.find(n2);

                if (n2EdgesIt != edgeMap.end()) {
                    const std::vector<Edge> n2Edges = n2EdgesIt->second;
                    edgeMap.erase(n2);

                    // Add edges to the subgraph
                    g.insert(g.end(), n2Edges.begin(), n2Edges.end());

                    // Add all connected edges to n2 to the queue
                    for (const Edge &e : n2EdgesIt->second) {
                        queue.emplace(e.targetId);
                    }
                }
            }
        }
        return subGraphs;
    }
};


