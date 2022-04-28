#include <iostream>
#include <fstream>
#include <queue>
#include <chrono>
#include "hnswlib/hnswlib.h"
#include <unordered_set>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "hnsw.h"

using namespace std;
using namespace hnswlib;






void hnsw_impl(string stage, string using_dataset, size_t data_size, size_t attr_size, size_t label_dim){
    string path_project = "..";
#if PLATG
    string label = "plat/";
#else
    string label = "base/";
#endif
    string path_graphindex = path_project + "/graphindex/" + label;

    string pre_index = path_graphindex + using_dataset;
    if (access(pre_index.c_str(), R_OK|W_OK)){
        if (mkdir(pre_index.c_str(), S_IRWXU) != 0) {
            printf("Error, dir %s create failed \n", pre_index.c_str());
            exit(1);
        }
    }

	size_t subset_size_milllions = data_size;
	size_t efConstruction = 200;
	size_t M = 20;
    size_t k = 10;

    size_t vecsize = subset_size_milllions * 1000000;
    size_t qsize, vecdim, gt_maxnum;
    string path_index, path_gt, path_q, path_data;

    std::map<string, size_t> mappmt;
    mappmt["subset_size_milllions"] = subset_size_milllions;
    mappmt["efConstruction"] = efConstruction;
    mappmt["M"] = M;
    mappmt["k"] = k;
    mappmt["vecsize"] = vecsize;
    mappmt["attr_size"] = attr_size;
    mappmt["label_dim"] = label_dim;

    std::map<string, string> mapstr;

    string hnsw_index = pre_index + "/" + using_dataset + to_string(subset_size_milllions) +
                        "m_ef" + to_string(efConstruction) + "m" + to_string(M) + "_attr_" + to_string(attr_size) + ".bin";
    mapstr["index"] = hnsw_index;
    CheckDataset(using_dataset, mappmt, mapstr);
    

    hnswMultiGraph<DTSET, DTRES> hnsw_multi_graph(mappmt, mapstr);

    if (stage == "build" || stage == "both")
        hnsw_multi_graph.build_index(attr_size);

    if (stage == "search" || stage == "both")
        hnsw_multi_graph.search_index(attr_size);

    return;
}

