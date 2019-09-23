//
// Created by remi on 19-09-19.
//

#pragma once

#include <vector>
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

    static inline std::vector<std::vector<unsigned int>> MCL(const std::vector<Edge> &edges, int id = 0) {
        std::vector<std::vector<unsigned int>> clusters;
        if (edges.size() == 0) return clusters;

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
        std::string cmd = "$HOME/local/bin/mcl " + fnIn + " --abc -I 2.0 -scheme 1 -te " + nT + " -o " + fnOut; // -q x -V all
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
};


