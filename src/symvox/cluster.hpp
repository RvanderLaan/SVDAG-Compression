//
// Created by remi on 19-09-19.
//

#pragma once

#include <vector>
#include <map>
#include <set>
#include <stack>
#include <algorithm>
#include <math.h>
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

    // Find only those edges that are in this specific subCluster
    // NOTE: Cluster needs to be sorted
    static inline void getClusterEdges(
        const std::vector<Edge> &subGraph,
        const std::vector<unsigned int> &cluster,
        std::vector<Edge> &clusterEdges
    ) {
        for (const auto &e : subGraph) {
            // Todo: This could be quite slow for large clusters. Make an edgeMap again?
            if (std::binary_search(cluster.begin(), cluster.end(), e.sourceId)
                && std::binary_search(cluster.begin(), cluster.end(), e.targetId)) {
                clusterEdges.emplace_back(e);
            }
        }
    }

    static inline std::vector<std::vector<unsigned int>> clusterSubgraphs(
            const std::vector<Edge> &edges,
            float inflation,
            int id = 0) {
        printf("Splitting %zu edges into subgraphs... ", edges.size()); fflush(stdout);
        const auto subGraphs = findSubGraphs(edges);
        printf("Found %zu subgraphs. ", subGraphs.size()); fflush(stdout);
        std::vector<std::vector<unsigned int>> clusters;
        std::vector<Edge> clusterEdges;

        // Threshold of cluster size for when to run MCL or when to just create a cluster of the whole subgraph
        double graphSizeThresholdToCluster = std::max(sqrt(edges.size()) / 4.0, 2.0);

        int stepLogger = (int) round(subGraphs.size() / 10.f);
        int i = 0;
        for (const auto& subGraph : subGraphs) {
            if (stepLogger > 0 && (i++ % stepLogger) == 0) {
                printf("%.0f%%..", round(100.0f * (i / (float) subGraphs.size()))); fflush(stdout);
            }
            if (subGraph.size() <= graphSizeThresholdToCluster) {
                // For extremly small clusters, the subgraph is one cluster
                std::set<unsigned int> clusterSet;
                for (const auto& e : subGraph) {
                    clusterSet.emplace(e.sourceId);
                    clusterSet.emplace(e.targetId);
                }
                std::vector<unsigned int> clusterVec(clusterSet.begin(), clusterSet.end());
                
                // Set cluster center at index 0
                setClusterCenterToBeginning(subGraph, clusterVec);

                clusters.emplace_back(clusterVec);
            } else {
                // Else, perform MCL
                auto subClusters = MCL(subGraph, inflation, id);

                // Set cluster center at index 0 for all subgraphs
                for (auto &cluster : subClusters) {
                    // Find only those edges that are in this specific subCluster
                    clusterEdges.clear();
                    getClusterEdges(subGraph, cluster, clusterEdges);
                    setClusterCenterToBeginning(clusterEdges, cluster);
                }

                clusters.insert(clusters.end(), subClusters.begin(), subClusters.end());
            }
        }
        printf("Total of %zu clusters. ", clusters.size());
        return clusters;
    }

    /**
     * Finds the node with the highest sum of edge weights and puts it at index 0 of the cluster
     * @param edges
     * @param cluster
     */
    static void setClusterCenterToBeginning(
            const std::vector<Edge> &edges,
            std::vector<unsigned int> &cluster
    ) {
        if (cluster.size() == 1) return;
        // Find cluster representative, with highest total edge weight
        std::map<unsigned int, float> summedEdgeWeights;
        for (const auto &e : edges) {
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
        auto max = std::max_element(summedEdgeWeights.begin(), summedEdgeWeights.end(),
                [](const auto& a, const auto& b)->bool{ return a.second < b.second; });
        auto maxIt = std::find(cluster.begin(), cluster.end(), max->first);

#if 1
        // DEBUG: Check whether the node with highest edge weight is in the cluster
        if (maxIt == cluster.end()) {
//            printf("Cluster: ");
//            for (auto i = cluster.begin(); i != cluster.end(); ++i)
//                printf("%u ", *i);
//            printf("\nEdges: ");
//            for (auto i = edges.begin(); i != edges.end(); ++i)
//                printf("%u->%u ", i->sourceId, i->targetId);
//            printf("\n");

//            printf("Cluster size: %zu, num edges: %zu\n", cluster.size(), edges.size());
//            printf("Cluster center not found?! This is likely an unconnected cluster, using first node as fallback. %u, %f\n", max->first, max->second);
            //            exit(-1);

            // Todo: Make each node a seperate cluster, or use the --force-connected=y option for MCL
            // Currently unconnected clusters are still used as a cluster with a random node as representative

            if (cluster.size() > 2)
                printf("Found unconnected cluster of size %zu. ", cluster.size());
            maxIt = cluster.begin();
        }
#endif

        unsigned int maxIndex = std::distance(cluster.begin(), maxIt);
        unsigned int tmp = cluster[0];
        cluster[0] = max->first;
        cluster[maxIndex] = tmp;
    }

    static inline void writeEdgesToFile(const std::vector<Edge> &edges, std::string fname) {
        std::ofstream edgesFile(fname, std::ios::out | std::ios::trunc);

        for (auto const& edge : edges) {
            edgesFile << std::to_string(edge.sourceId) << " "
                   << std::to_string(edge.targetId) << " "
                   << edge.weight << "\n"; stdout;
        }
        edgesFile.close();
    }

    static inline std::vector<std::vector<unsigned int>> MCL(
            const std::vector<Edge> &edges,
            float inflation,
            int id = 0
    ) {
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
        std::string I = std::to_string(inflation);
        std::string cmd = "$HOME/local/bin/mcl " + fnIn + " --abc -I " + I + " -q x -V all -te " + nT + " -o " + fnOut; // -scheme 1
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

            std::vector<unsigned int> cluster;

            // mcl uses tab, hipmcl uses space
            const char delim = '\t'; // ' '

            while (std::getline(linestream, substr, delim)) {
                cluster.emplace_back(std::stoul(substr));
            }
            std::sort(cluster.begin(), cluster.end());
            clusters.emplace_back(cluster);
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
                        queue.emplace(n2 == e.sourceId ? e.targetId : e.sourceId);
                    }
                }
            }
            std::vector<Edge> gVec(g.begin(), g.end());

#if 0
//            printf("DEBUG: Check for duplicates in other graphs\n");
            unsigned int num_dupes = 0;
            for (const auto &og : subGraphs) {
                for (const auto &oe : og) {
                    if (g.find(oe) != g.end()) {
                        num_dupes++;
                    }
                }
            }
            if (num_dupes > 0) printf("DUPES: %u ", num_dupes); fflush(stdout);
#endif

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

        unsigned int numEdgesInSubgraphs = 0;
        std::vector<unsigned int> bins(16);
        unsigned int maxBin = 10;
        for (const auto &g : subGraphs) {
            numEdgesInSubgraphs += g.size();

            bins[std::min(g.size() / maxBin, bins.size() - 1)] += 1;
        }
        printf("Total of %u edges in all subgraphs\n", numEdgesInSubgraphs);
        printf("Subgraph sizes for %s \n", edgeFn.c_str());

        for (int i = 0; i < bins.size(); i++) {
            printf("%i, ", i + 1);
        }
        printf("\n");
        for (int i = 0; i < bins.size(); i++) {
            printf("%u, ", bins[i]);
        }
        printf("\n");
     }
};



     // Batch small subgraphs together to avoid IO overhead of writing/reading to/from disk
        // int maxBatchSize = 50000;

        // std::vector<Edge> batch;
        // batch.reserve(maxBatchSize);

        /**
         * Flow:
         * 1. Write to file 1
         * 2. For every graph
         * 2.1. Write to file 2 in thread
         * 2.2. Cluster file 1
         * 2.3. Wait for both...
         * 2.4. Write to file 1
         * 2.5. Cluster file 2
         * 2.6. Wait for both...
         */

        // int stepLogger = subGraphs.size() / 10;
        // int i = 0;
        // for (const auto& subGraph : subGraphs) {
        //     if ((i++ % stepLogger) == 0) {
        //         printf("%.0f%%..", round(100.0f * (i / (float) subGraphs.size()))); fflush(stdout);
        //     }
        //     // Todo: For extremly small clusters, the subgraph is one cluster? Needs experiment
        //     // if (subGraph.size() <= graphSizeThresholdToCluster) {
        //     if (subGraph.size() > maxBatchSize) {
        //         auto subClusters = MCL(subGraph, inflation, id);
        //         // Set cluster center at index 0 for all subgraphs
        //         for (auto &cluster : subClusters) {
        //             clusterEdges.clear();
        //             getClusterEdges(subGraph, cluster, clusterEdges);
        //             setClusterCenterToBeginning(clusterEdges, cluster);
        //         }
        //         clusters.insert(clusters.end(), subClusters.begin(), subClusters.end());

        //     } else if (batch.size() + subGraph.size() < maxBatchSize) {
        //         // Add subgraph to batch if it fits
        //         batch.insert(batch.end(), subGraph.begin(), subGraph.end());
        //     } else {
        //         // Cluster batch if it is full
        //         auto subClusters = MCL(batch, inflation, id);
        //         for (auto &cluster : subClusters) {
        //             clusterEdges.clear();
        //             getClusterEdges(batch, cluster, clusterEdges);
        //             setClusterCenterToBeginning(clusterEdges, cluster);
        //         }
        //         clusters.insert(clusters.end(), subClusters.begin(), subClusters.end());

        //         // Add current subGraph back into the batch
        //         batch.clear();
        //         batch.insert(batch.end(), subGraph.begin(), subGraph.end());
        //     }
        // }
        // // Cluster what remains in the batch
        // if (batch.size() > 0) {
        //     auto subClusters = MCL(batch, inflation, id);
        //     for (auto &cluster : subClusters) {
        //         clusterEdges.clear();
        //         getClusterEdges(batch, cluster, clusterEdges);
        //         setClusterCenterToBeginning(clusterEdges, cluster);
        //     }
        //     clusters.insert(clusters.end(), subClusters.begin(), subClusters.end());
        // }
