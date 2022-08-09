
#include <iostream>
#include <fstream>
#include <queue>
#include <chrono>
#include "hnswlib/hnswlib.h"
#include <unordered_set>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <cmath>

using namespace std;
using namespace hnswlib;

class Timer {
    std::chrono::high_resolution_clock::time_point time_begin;
public:
    Timer() {
        time_begin = std::chrono::high_resolution_clock::now();
    }

    float getElapsedTimeus() {
        std::chrono::high_resolution_clock::time_point time_end = std::chrono::high_resolution_clock::now();
        return (std::chrono::duration_cast<std::chrono::microseconds>(time_end - time_begin).count());
    }

    float getElapsedTimes() {
        return 1e-6 * getElapsedTimeus();
    }

    void reset() {
        time_begin = std::chrono::high_resolution_clock::now();
    }

};



/*
* Author:  David Robert Nadeau
* Site:    http://NadeauSoftware.com/
* License: Creative Commons Attribution 3.0 Unported License
*          http://creativecommons.org/licenses/by/3.0/deed.en_US
*/

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))

#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif


/**
* Returns the peak (maximum so far) resident set size (physical
* memory use) measured in bytes, or zero if the value cannot be
* determined on this OS.
*/
static size_t getPeakRSS() {
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ((fd = open("/proc/self/psinfo", O_RDONLY)) == -1)
        return (size_t)0L;      /* Can't open? */
    if (read(fd, &psinfo, sizeof(psinfo)) != sizeof(psinfo))
    {
        close(fd);
        return (size_t)0L;      /* Can't read? */
    }
    close(fd);
    return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
    return (size_t) (rusage.ru_maxrss * 1024L);
#endif

#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L;          /* Unsupported. */
#endif
}


/**
* Returns the current resident set size (physical memory use) measured
* in bytes, or zero if the value cannot be determined on this OS.
*/
static size_t getCurrentRSS() {
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
    /* OSX ------------------------------------------------------ */
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
        (task_info_t)&info, &infoCount) != KERN_SUCCESS)
        return (size_t)0L;      /* Can't access? */
    return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    /* Linux ---------------------------------------------------- */
    long rss = 0L;
    FILE *fp = NULL;
    if ((fp = fopen("/proc/self/statm", "r")) == NULL)
        return (size_t) 0L;      /* Can't open? */
    if (fscanf(fp, "%*s%ld", &rss) != 1) {
        fclose(fp);
        return (size_t) 0L;      /* Can't read? */
    }
    fclose(fp);
    return (size_t) rss * (size_t) sysconf(_SC_PAGESIZE);

#else
    /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    return (size_t)0L;          /* Unsupported. */
#endif
}


template<typename DTres>
static void
get_gt(unsigned *massQA, size_t qsize, size_t &gt_maxnum, size_t vecdim,
        vector<std::priority_queue<std::pair<DTres, labeltype >>> &answers, size_t k) {
    (vector<std::priority_queue<std::pair<DTres, labeltype >>>(qsize)).swap(answers);
    cout << qsize << "\n";
    for (int i = 0; i < qsize; i++) {
        for (int j = 0; j < k; j++) {
            answers[i].emplace(0.0f, massQA[gt_maxnum * i + j]);
        }
    }
}

void SplitString(const string &s, vector<int> &v, const string &c)
{
    string::size_type pos1, pos2;
    pos2 = s.find(c);
    pos1 = 0;

    while (string::npos != pos2)
    {
        v.push_back(atoi(s.substr(pos1, pos2 - pos1).c_str()));

        pos1 = pos2 + c.size();
        pos2 = s.find(c, pos1);
    }
    if (pos1 != s.length())
        v.push_back(atoi(s.substr(pos1).c_str()));
}

void SplitString(const std::string &s, std::vector<std::string> &v, const std::string &c)
{
    std::string::size_type pos1, pos2;
    pos2 = s.find(c);
    pos1 = 0;

    while (std::string::npos != pos2)
    {
        v.push_back(s.substr(pos1, pos2 - pos1).c_str());

        pos1 = pos2 + c.size();
        pos2 = s.find(c, pos1);
    }
    if (pos1 != s.length())
        v.push_back(s.substr(pos1).c_str());
}

void load_data_txt(string filename, unsigned &num, unsigned &dim, std::vector<std::vector<std::string>> &data)
{
    std::string temp;
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cout << "open file error : " << filename << std::endl;
        exit(-1);
    }
    getline(file, temp);
    std::vector<int> tmp2;
    SplitString(temp, tmp2, " ");
    num = tmp2[0];
    dim = tmp2[1];
    data.resize(num);
    int groundtruth_count = 0;
    while (getline(file, temp))
    {
        SplitString(temp, data[groundtruth_count], " ");
        groundtruth_count++;
    }
    std::cout << "Load " << data.size() << " Data from " << filename << std::endl;
    file.close();
}



template<typename DTset, typename DTres>
class hnswMultiGraph {
public:
    size_t efConstruction;
    size_t M1;
    size_t M2;
    size_t vecsize;
    size_t vecdim;

    string path_data;
    string path_attr;
    string index;

    size_t k;
    size_t qsize;
    size_t gt_maxnum;

    string path_q;
    string path_q_attr;
    string path_gt;

    unsigned attr_num, attr_dim, bucket_num;

    HierarchicalNSW<DTres> *appr_alg;
    std::vector<HierarchicalNSW<DTres> *> appr_alg_bucket;
    HierarchicalNSW<DTres> *appr_alg_center;

    DTset *massB;
    DTset *massbucket;
    DTset *masscenter;
    DTset *massQ;

    char * data_memory;

    std::vector<std::vector<string>> attr_data;
    std::vector<std::vector<string>> attr_query;
    std::vector<unsigned> bucket_id;
    std::vector<std::vector<unsigned>> bucket_id_query;
    std::vector<unsigned> bucket_size;
    std::vector<std::vector<unsigned>> bucket;
    std::vector<unsigned> center_id_graph;
    std::vector<unsigned> center_id_data;

    size_t offsetLevel0;
    size_t max_elements;
    size_t cur_element_count;
    size_t size_data_per_element;
    size_t data_offset;
    size_t bucket_links_offset;
    size_t label_offset;
    size_t bucketid_offset;
    int maxlevel;
    hnswlib::tableint enterpoint_node;
    unsigned search_enterpoint_node_num;
    size_t M;
    size_t maxM;
    size_t maxM0;
    double mult;
    size_t ef_construction;

    unsigned enterpoint_num;

    // according to the index order in the global graph
    std::vector<unsigned> list_bucket_id;
    std::vector<hnswlib::linklistsizeint> list_nei_num;
    std::vector<std::vector<hnswlib::tableint>> list_nei;




    hnswMultiGraph(map<string, size_t> &mappmt, map<string, string> &mapstr) {
        efConstruction = mappmt["efConstruction"];
        M1 = mappmt["M1"];
        M2 = mappmt["M2"];
        vecsize = mappmt["vecsize"];
        vecdim = mappmt["vecdim"];

        path_data = mapstr["path_data"];
        path_attr = mapstr["path_attr"];
        index = mapstr["index"];

        k = mappmt["k"];
        qsize = mappmt["qsize"];
        gt_maxnum = mappmt["gt_maxnum"];

        path_q = mapstr["path_q"];
        path_q_attr = mapstr["path_q_attr"];
        path_gt = mapstr["path_gt"];
    }

    ~hnswMultiGraph() {}

    void testcode(std::vector<std::bitset<BIT_SET_SIZE>> bit_table) {
    }

    


    void test_approx(std::vector<std::priority_queue<std::pair<DTres, labeltype >>> &result, float &time_total, float &time_lookup_table) {
        std::vector<bool> bucket_id_query_tmp(bucket_num);
        enterpoint_num = 1;
        bool strategy = attr_dim < (log2(bucket_num) - 1) * 2 ? 0 : 1;

    //     omp_set_num_threads(3);
    // #pragma omp parallel for
        for (int i = 0; i < qsize; i++) {
            Timer stopw = Timer();
            bucket_id_query_tmp.assign(bucket_num, false);
            for (int j = 0; j < attr_dim; j++) {
                bucket_id_query_tmp[bucket_id_query[i][j]] = true;
            }
            time_lookup_table += (double)stopw.getElapsedTimeus();

            result.push_back(appr_alg->searchKnn(massQ + vecdim * i, k, bucket_id_query_tmp, enterpoint_num, strategy));
            time_total += (double)stopw.getElapsedTimeus();
        }
    }

    float get_recall(vector<std::priority_queue<std::pair<DTres, labeltype >>> &answers, std::vector<std::priority_queue<std::pair<DTres, labeltype >>> &result) {
        size_t correct = 0;
        size_t total = 0;

        for (int i = 0; i < qsize; i++) {
            std::priority_queue<std::pair<DTres, labeltype >> gt(answers[i]);
            unordered_set<labeltype> g;

            while (gt.size()) {
                g.insert(gt.top().second);
                gt.pop();
            }
            {
                total += g.size();
                while (result[i].size()) {
                    if (g.find(result[i].top().second) != g.end()) {
                        correct++;
                    }
                    result[i].pop();
                }
            }
        }
        return 1.0f * correct / total;
    }


    void test_vs_recall(vector<std::priority_queue<std::pair<DTres, labeltype >>> &answers) {
        vector<size_t> efs;
        for (int i = 10; i <= 200; i += 10)
            efs.push_back(i);

        cout    << "efs\t" 
                << "R@" << k << "\t" 
                << "QPS\t" 
                << "time_us\t" 
                // << "time_lt_us\t" 
                << "NDC_avg\t" 
                << endl;

        // for (int i = 0; i < 5;) {
        //     appr_alg->setEf(20 * (++i));
        //     appr_alg->metric_hops = 0;
        //     appr_alg->metric_distance_computations = 0;
        //     std::vector<std::priority_queue<std::pair<DTres, labeltype >>> result;

        //     float time_us_total = 0.0;
        //     float time_us_lookup_table = 0.0;
        //     test_approx(result, time_us_total, time_us_lookup_table);
        // }

        for (size_t ef : efs) {
            appr_alg->setEf(ef);
            appr_alg->metric_hops = 0;
            appr_alg->metric_distance_computations = 0;
            std::vector<std::priority_queue<std::pair<DTres, labeltype >>> result;

            float time_us_total = 0.0;
            float time_us_lookup_table = 0.0;
            test_approx(result, time_us_total, time_us_lookup_table);
            float time_us_per_query = time_us_total / qsize;
            float time_lt_us_per_query = time_us_lookup_table / qsize;

            float hop_avg = 1.0f * appr_alg->metric_hops / qsize;
            float NDC_avg = 1.0f * appr_alg->metric_distance_computations / qsize;

            float recall = get_recall(answers, result);

            cout    << ef << "\t" 
                    << recall << "\t" 
                    << 1 / time_us_per_query * 1e6 << "\t" 
                    << time_us_per_query << "\t" 
                    // << time_lt_us_per_query << "\t" 
                    << NDC_avg << "\t" 
                    << endl;

            if (recall > 1.0) {
                cout << recall << "\t" << time_us_per_query << " us\n";
                break;
            }
        }
    }

    inline bool exists_test(const std::string &name) {
        ifstream f(name.c_str());
        return f.good();
    }

    


    void build_index_1(HierarchicalNSW<DTres> *appr_alg, size_t vecsize, size_t vecdim){
#if PLATG
        unsigned center_id = compArrayCenter<DTset>(massB, vecsize, vecdim);
        appr_alg->addPoint((void *) (massB + center_id * vecdim), (size_t) center_id);
#else
        appr_alg->addPoint((void *) (massB), (size_t) 0);
#endif
        cout << "Building index:\n";
        int j1 = 0;
        Timer stopw = Timer();
        Timer stopw_full = Timer();
        size_t report_every = vecsize / 10;
#pragma omp parallel for
        for (size_t i = 1; i < vecsize; i++) {
#pragma omp critical
            {
                j1++;
                if (j1 % report_every == 0) {
                    cout << j1 / (0.01 * vecsize) << " %, "
                            << report_every / (1000.0 * stopw.getElapsedTimes()) << " kips " << " Mem: "
                            << getCurrentRSS() / 1000000 << " Mb \n";
                    stopw.reset();
                }
            }
#if PLATG
            size_t ic;
            if (i <= center_id)
                ic = i - 1;
            else
                ic = i;
            appr_alg->addPoint((void *) (massB + ic * vecdim), ic);
#else
            appr_alg->addPoint((void *) (massB + i * vecdim), i);
#endif
        }
        cout << "Build time:" << stopw_full.getElapsedTimes() << "  seconds\n";
    }


    unsigned build_index_2(HierarchicalNSW<DTres> *appr_alg, size_t vecsize, size_t vecdim, std::vector<unsigned> bucket){
        unsigned center_id_bucket = compArrayCenter<DTset>(massbucket, vecsize, vecdim);
        unsigned center_id = bucket[center_id_bucket];
        appr_alg->addPoint((void *) (massB + center_id * vecdim), (size_t) center_id);

        cout << "Building index:\n";
        int j1 = 0;
        Timer stopw = Timer();
        Timer stopw_full = Timer();
        size_t report_every = vecsize;
#pragma omp parallel for
        for (size_t j = 0; j < vecsize; j++) {
#pragma omp critical
            {
                j1++;
                if (j1 % report_every == 0) {
                    cout << j1 / (0.01 * vecsize) << " %, "
                        << report_every / (1000.0 * stopw.getElapsedTimes()) << " kips " << " Mem: "
                        << getCurrentRSS() / 1000000 << " Mb \n";
                    stopw.reset();
                }
            }

            size_t ic = bucket[j]; 
            if (j != center_id_bucket)
                appr_alg->addPoint((void *) (massB + ic * vecdim), ic);

        }
        cout << "Build time:" << stopw_full.getElapsedTimes() << "  seconds\n";

        return center_id;
    }

//     void build_index_3(HierarchicalNSW<DTres> *appr_alg, size_t vecsize, size_t vecdim, std::vector<unsigned> centers){
//         unsigned center_id_in_centers = compArrayCenter<DTset>(masscenter, vecsize, vecdim);
//         unsigned center_id = centers[center_id_in_centers];
//         appr_alg->addPoint((void *) (massB + center_id * vecdim), (size_t) center_id);

//         cout << "Building index:\n";
//         int j1 = 0;
//         Timer stopw = Timer();
//         Timer stopw_full = Timer();
//         size_t report_every = vecsize;
// #pragma omp parallel for
//         for (size_t j = 0; j < vecsize; j++) {
// #pragma omp critical
//             {
//                 j1++;
//                 if (j1 % report_every == 0) {
//                     cout << j1 / (0.01 * vecsize) << " %, "
//                         << report_every / (1000.0 * stopw.getElapsedTimes()) << " kips " << " Mem: "
//                         << getCurrentRSS() / 1000000 << " Mb \n";
//                     stopw.reset();
//                 }
//             }

//             size_t ic = centers[j]; 
//             if (j != center_id_in_centers)
//                 appr_alg->addPoint((void *) (massB + ic * vecdim), ic);

//         }
//         cout << "Build time:" << stopw_full.getElapsedTimes() << "  seconds\n";

//     }


    

    void config_index_new(){
        ///
        size_t size_links = (M1 + M2) * 2 * sizeof(hnswlib::tableint) + sizeof(hnswlib::linklistsizeint);

        offsetLevel0 = appr_alg->offsetLevel0_;
        max_elements = appr_alg->max_elements_;
        cur_element_count = appr_alg->cur_element_count;
        size_data_per_element = size_links + appr_alg->data_size_ + sizeof(hnswlib::labeltype) + sizeof(unsigned);
        data_offset = size_links;
        label_offset = size_links + appr_alg->data_size_;
        bucketid_offset = size_links + appr_alg->data_size_ + sizeof(hnswlib::labeltype);
        maxlevel = appr_alg->maxlevel_;
        enterpoint_node = appr_alg->enterpoint_node_;
        search_enterpoint_node_num = bucket_num;
        maxM = (M1 + M2);
        maxM0 = (M1 + M2) * 2;
        M = (M1 + M2);
        mult = appr_alg->mult_;
        ef_construction = appr_alg->ef_construction_;

        data_memory = (char *) malloc(vecsize * size_data_per_element);


#pragma omp parallel for
        for (size_t i = 0; i < vecsize; i++) {
            ///
            hnswlib::linklistsizeint linklistsize[2];
            hnswlib::tableint link_id;
            hnswlib::labeltype link_id_data;
            unsigned link_bucket_id;
            unsigned bucket_vecsize = bucket_size[i];
            
            std::vector<hnswlib::tableint> link_id_table;
            
            linklistsize[1] = *(hnswlib::linklistsizeint *)(appr_alg->data_level0_memory_ + i * appr_alg->size_data_per_element_);
            
            for (int j = 0; j < linklistsize[1]; j++) {
                memcpy(&link_id, appr_alg->data_level0_memory_ + i * appr_alg->size_data_per_element_ + sizeof(hnswlib::linklistsizeint) + j * sizeof(hnswlib::tableint), sizeof(hnswlib::tableint));
                if (find(link_id_table.begin(), link_id_table.end(), link_id) == link_id_table.end())
                    link_id_table.push_back(link_id);
            }
            
            for (int j = 0; j < list_nei_num[i]; j++) {
                link_id = list_nei[i][j];
                if (find(link_id_table.begin(), link_id_table.end(), link_id) == link_id_table.end())
                    link_id_table.push_back(link_id);
            }

            linklistsize[0] = std::min(link_id_table.size(), (M1 + M2) * 2);
            memcpy(data_memory + i * size_data_per_element, linklistsize, sizeof(hnswlib::linklistsizeint));
            for (int j = 0; j < linklistsize[0]; j++) {
                link_id = link_id_table[j];
                link_bucket_id = list_bucket_id[i];
                memcpy(data_memory + i * size_data_per_element + sizeof(hnswlib::linklistsizeint) + j * sizeof(hnswlib::tableint), &link_id, sizeof(hnswlib::tableint));
            }
            memcpy(data_memory + i * size_data_per_element + data_offset, appr_alg->data_level0_memory_ + i * appr_alg->size_data_per_element_ + appr_alg->offsetData_, appr_alg->data_size_ + sizeof(hnswlib::labeltype));
            memcpy(data_memory + i * size_data_per_element + bucketid_offset, &link_bucket_id, sizeof(unsigned));
        }

        
    }

    void save_index(const std::string location){
        std::ofstream output(location, std::ios::binary);
        writeBinaryPOD(output, offsetLevel0);
        writeBinaryPOD(output, max_elements);
        writeBinaryPOD(output, cur_element_count);
        writeBinaryPOD(output, size_data_per_element);
        writeBinaryPOD(output, label_offset);
        writeBinaryPOD(output, bucketid_offset);
        writeBinaryPOD(output, data_offset);
        writeBinaryPOD(output, bucket_links_offset);
        writeBinaryPOD(output, maxlevel);
        writeBinaryPOD(output, enterpoint_node);
        writeBinaryPOD(output, search_enterpoint_node_num);
        for (unsigned i = 0; i < search_enterpoint_node_num; i++) {
            writeBinaryPOD(output, center_id_graph[i]);
        }
        writeBinaryPOD(output, maxM);

        writeBinaryPOD(output, maxM0);
        writeBinaryPOD(output, M);
        writeBinaryPOD(output, mult);
        writeBinaryPOD(output, ef_construction);

        output.write(data_memory, cur_element_count * size_data_per_element);

        output.close();
        free(data_memory);
    }

    void build_index(size_t attr_size, bool isSave = true){

        if (exists_test(index)){
            printf("Index %s is existed \n", index.c_str());
            return;
        } else 
        {
            
            
            massB = new DTset[vecsize * vecdim]();
            cout << "Loading base data:\n";
            LoadBinToArray<DTset>(path_data, massB, vecsize, vecdim);

#if FMTINT
            L2SpaceI l2space(vecdim);
#else
            L2Space l2space(vecdim);
#endif

            auto s0 = std::chrono::high_resolution_clock::now();
            appr_alg = new HierarchicalNSW<DTres>(&l2space, vecsize, M1, efConstruction);
            auto s = std::chrono::high_resolution_clock::now();
            build_index_1(appr_alg, vecsize, vecdim);
            auto e = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_index = e - s;

            

            ///            
            cout << "Loading base data attribute:\n";
            load_data_txt(path_attr, attr_num, attr_dim, attr_data);
            
            for (size_t i = 0; i < vecsize; i++) {
                bucket_id.push_back(atoi(attr_data[i][0].c_str()));
            }

            bucket_num = attr_size;

            for (size_t i = 0; i < bucket_num; i++) {
                std::vector<unsigned> bucket_tmp;
                bucket.push_back(bucket_tmp);
            }
            for (size_t i = 0; i < vecsize; i++) {
                bucket[bucket_id[i]].push_back(i);
                unsigned l_bucket_id;
                list_bucket_id.push_back(l_bucket_id);
                hnswlib::linklistsizeint l_nei_num;
                list_nei_num.push_back(l_nei_num);
                std::vector<hnswlib::tableint> l_nei;
                list_nei.push_back(l_nei);
            }

            for (size_t i = 0; i < bucket_num; i++) {
                bucket_size.push_back(bucket[i].size());
                unsigned bucket_vecsize = bucket_size[i];
                HierarchicalNSW<DTres> *appr_alg_bucket_new = new HierarchicalNSW<DTres>(&l2space, bucket_vecsize, M2, efConstruction);
                massbucket = new DTset[bucket_vecsize * vecdim]();
                for (size_t j = 0; j < bucket_vecsize; j++) {
                    memcpy(massbucket + j * vecdim, massB + bucket[i][j] * vecdim, vecdim * sizeof(DTset));
                }
                s = std::chrono::high_resolution_clock::now();
                center_id_data.push_back(build_index_2(appr_alg_bucket_new, bucket_vecsize, vecdim, bucket[i]));
                e = std::chrono::high_resolution_clock::now();
                time_index += e - s;

                auto search = appr_alg->label_lookup_.find(center_id_data[i]);
                center_id_graph.push_back(search->second);

                hnswlib::labeltype local_id_data;
                hnswlib::tableint link_id;
                hnswlib::labeltype link_id_data;
                for (size_t j = 0; j < bucket_vecsize; j++) {
                    local_id_data = appr_alg_bucket_new->getExternalLabel(j);
                    auto global_search = appr_alg->label_lookup_.find(local_id_data);
                    hnswlib::tableint global_id = global_search->second;
                    list_bucket_id[global_id] = i;
                    list_nei_num[global_id] = *(hnswlib::linklistsizeint *)(appr_alg_bucket_new->data_level0_memory_ + j * appr_alg_bucket_new->size_data_per_element_);
                    for (int k = 0; k < list_nei_num[global_id]; k++) {
                        memcpy(&link_id, appr_alg_bucket_new->data_level0_memory_ + j * appr_alg_bucket_new->size_data_per_element_ + sizeof(hnswlib::linklistsizeint) + k * sizeof(hnswlib::tableint), sizeof(hnswlib::tableint));
                        link_id_data = appr_alg_bucket_new->getExternalLabel(link_id);
                        auto link_search = appr_alg->label_lookup_.find(link_id_data);
                        hnswlib::tableint link_id_graph = link_search->second;
                        list_nei[global_id].push_back(link_id_graph);
                    }
                }

                delete appr_alg_bucket_new;

                delete[] massbucket;
            }


            if (isSave) {
                s = std::chrono::high_resolution_clock::now();
                // config_index();
                config_index_new();
                e = std::chrono::high_resolution_clock::now();
                time_index += e - s;
                printf("Build index time : %lf\n", time_index.count());
                s = std::chrono::high_resolution_clock::now();
                save_index(index);
                e = std::chrono::high_resolution_clock::now();
                time_index = e - s;
                printf("Save index time : %lf\n", time_index.count());
            }
            else {
                printf("Build index time : %lf\n", time_index.count());
            }
            auto e0 = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_all = e0 - s0;
            printf("Build index total time : %lf\n", time_all.count());
            
            delete[] massB;
            delete appr_alg;
            
            printf("Build index %s is succeed \n", index.c_str());
            
        }
    }


    void search_index(size_t attr_size){

        if (!exists_test(index)){
            printf("Error, index %s is unexisted \n", index.c_str());
            exit(1);
        } else {

            unsigned *massQA = new unsigned[qsize * gt_maxnum];
            massQ = new DTset[qsize * vecdim];
            
            bucket_num = attr_size;

            cout << "Loading GT:\n";
            LoadBinToArray<unsigned>(path_gt, massQA, qsize, gt_maxnum);
            cout << "Loading queries:\n";
            LoadBinToArray<DTset>(path_q, massQ, qsize, vecdim);
            cout << "Loading queries attribute:\n";
            load_data_txt(path_q_attr, attr_num, attr_dim, attr_query);
            if (qsize != attr_num) {
                printf("Error, number of queries(%d) and attributes(%d) do not match \n", qsize, attr_num);
                exit(1);
            }

            

            for (size_t i = 0; i < attr_num; i++) {
                std::vector<unsigned> tmp;
                for (size_t j = 0; j < attr_dim; j++) {
                    tmp.push_back(atoi(attr_query[i][j].c_str()));
                }
                bucket_id_query.push_back(tmp);
            }

#if FMTINT
            L2SpaceI l2space(vecdim);
#else
            L2Space l2space(vecdim);
#endif
            appr_alg = new HierarchicalNSW<DTres>(&l2space, index, false);

            vector<std::priority_queue<std::pair<DTres, labeltype >>> answers;
            cout << "Parsing gt:\n";
            get_gt(massQA, qsize, gt_maxnum, vecdim, answers, k);

            cout << "Comput recall: \n";
            test_vs_recall(answers);

            printf("Search index %s is succeed \n", index.c_str());

            delete[] massQA;
            delete[] massQ;
            delete appr_alg;
        }
    }


};