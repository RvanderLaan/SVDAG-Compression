//
// Created by remi on 19-09-19.
//

#pragma once

#include <vector>


#include <Eigen/Core>
#include <Eigen/Sparse>
#include <fstream>

typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> Triplet;

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
        std::string cmd = "$HOME/local/bin/mcl " + fnIn + " --abc -I 2.0 -te " + nT + " -o " + fnOut; // -q x -V all
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

    cluster(unsigned int nNodes) : mat(nNodes, nNodes) { }

    inline std::vector<std::vector<unsigned int>> FindClusters(const std::vector<Edge> &edges, unsigned int nNodes) {

        int nEdges = edges.size();

        std::map<unsigned int, unsigned int> lookup; // map of node index to matrix index
        std::vector<Triplet> triplets(nEdges);

        unsigned int i = 0;
        for (auto const& edge : edges) {
            unsigned int matSourceId = lookup.size();
            if (!lookup.insert(std::make_pair(edge.sourceId, matSourceId)).second) { // exists
                matSourceId = lookup[edge.sourceId];
            }
            unsigned int matTargetId = lookup.size();
            if (!lookup.insert(std::make_pair(edge.targetId, matTargetId)).second) { // exists
                matSourceId = lookup[edge.targetId];
            }
            triplets.emplace_back(Triplet(matSourceId, matTargetId, edge.weight));
        }
        mat.setFromTriplets(triplets.begin(), triplets.end());

    }

private:
    SpMat mat;

    void normalize(SpMat& in) {
        for (unsigned int i = 0; i < in.cols(); ++i) {
            double sum = 0;
            for (unsigned int j = 0; j < in.rows(); ++j) {
                sum += in.coeff(i, j);
            }
            // Todo: Divide by sum
        }
    }

    void expand() {

    }

    void inflate(SpMat& in, double inflation) {
        auto lam = [inflation](double x) -> double { return std::pow(x, inflation); };
        in = in.unaryExpr(lam);

    }
};


