#include <iostream>
#include <fstream>
#include <queue>
#include <chrono>
#include "hnswlib/hnswlib.h"
#include <unordered_set>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

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

template<typename DTset, typename DTres>
static float
test_approx(DTset *massQ, size_t qsize, HierarchicalNSW<DTres> &appr_alg, size_t vecdim,
            vector<std::priority_queue<std::pair<DTres, labeltype >>> &answers, size_t k) {
    size_t correct = 0;
    size_t total = 0;

//     omp_set_num_threads(3);
// #pragma omp parallel for
    for (int i = 0; i < qsize; i++) {
        std::priority_queue<std::pair<DTres, labeltype >> result = appr_alg.searchKnn(massQ + vecdim * i, k);

        std::priority_queue<std::pair<DTres, labeltype >> gt(answers[i]);
        unordered_set<labeltype> g;

        while (gt.size()) {
            g.insert(gt.top().second);
            gt.pop();
        }

// #pragma omp critical
        {
            total += g.size();
            while (result.size()) {
                if (g.find(result.top().second) != g.end()) {
                    correct++;
                }
                result.pop();
            }
        }
    }
    return 1.0f * correct / total;
}

template<typename DTset, typename DTres>
static void
test_vs_recall(DTset *massQ, size_t qsize, HierarchicalNSW<DTres> &appr_alg, size_t vecdim,
               vector<std::priority_queue<std::pair<DTres, labeltype >>> &answers, size_t k) {
    vector<size_t> efs;// = { 10,10,10,10,10 };
    for (int i = 10; i <= 100; i += 10)
        efs.push_back(i);
    cout << "efs\t" << "R@" << k << "\t" << "time_us\t" << "NDC_avg" << endl;

    for (size_t ef : efs) {
        appr_alg.setEf(ef);
        appr_alg.metric_hops = 0;
        appr_alg.metric_distance_computations = 0;

        Timer stopw = Timer();
        float recall = test_approx(massQ, qsize, appr_alg, vecdim, answers, k);
        float time_us_per_query = (double) stopw.getElapsedTimeus() / qsize;

        float hop_avg = 1.0f * appr_alg.metric_hops / qsize;
        float NDC_avg = 1.0f * appr_alg.metric_distance_computations / qsize;

        cout << ef << "\t" << recall << "\t" << time_us_per_query << "\t" << NDC_avg << endl;

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


template<typename DTset, typename DTres>
void build_index(map<string, size_t> &mappmt, map<string, string> &mapstr, bool isSave = true){
    //
    size_t efConstruction = mappmt["efConstruction"];
    size_t M = mappmt["M"];
    size_t vecsize = mappmt["vecsize"];
    size_t vecdim = mappmt["vecdim"];
    size_t qsize = mappmt["qsize"];

    string path_data = mapstr["path_data"];
    string index = mapstr["index"];

    if (exists_test(index)){
        printf("Index %s is existed \n", index.c_str());
        return;
    } else {

        DTset *massB = new DTset[vecsize * vecdim]();
        cout << "Loading base data:\n";
        LoadBinToArray<DTset>(path_data, massB, vecsize, vecdim);

#if FMTINT
        L2SpaceI l2space(vecdim);
#else
        L2Space l2space(vecdim);
#endif

        HierarchicalNSW<DTres> *appr_alg = new HierarchicalNSW<DTres>(&l2space, vecsize, M, efConstruction);

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
        delete[] massB;
        if (isSave)
            appr_alg->saveIndex(index);

        printf("Build index %s is succeed \n", index.c_str());
    }
}

template<typename DTset, typename DTres>
void search_index(map<string, size_t> &mappmt, map<string, string> &mapstr){
    //
    size_t k = mappmt["k"];
    size_t vecsize = mappmt["vecsize"];
    size_t qsize = mappmt["qsize"];
    size_t vecdim = mappmt["vecdim"];
    size_t gt_maxnum = mappmt["gt_maxnum"];

    string path_q = mapstr["path_q"];
    string index = mapstr["index"];
    string path_gt = mapstr["path_gt"];

    if (!exists_test(index)){
        printf("Error, index %s is unexisted \n", index.c_str());
        exit(1);
    } else {

        unsigned *massQA = new unsigned[qsize * gt_maxnum];
        DTset *massQ = new DTset[qsize * vecdim];

        cout << "Loading GT:\n";
        LoadBinToArray<unsigned>(path_gt, massQA, qsize, gt_maxnum);
        cout << "Loading queries:\n";
        LoadBinToArray<DTset>(path_q, massQ, qsize, vecdim);

#if FMTINT
        L2SpaceI l2space(vecdim);
#else
        L2Space l2space(vecdim);
#endif
        HierarchicalNSW<DTres> *appr_alg = new HierarchicalNSW<DTres>(&l2space, index, false);

        vector<std::priority_queue<std::pair<DTres, labeltype >>> answers;
        cout << "Parsing gt:\n";
        get_gt(massQA, qsize, gt_maxnum, vecdim, answers, k);

        cout << "Comput recall: \n";
        test_vs_recall(massQ, qsize, *appr_alg, vecdim, answers, k);

        printf("Search index %s is succeed \n", index.c_str());
    }
}


void hnsw_impl(string stage, string using_dataset, size_t data_size){
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

    std::map<string, string> mapstr;

    string hnsw_index = pre_index + "/" + using_dataset + to_string(subset_size_milllions) +
                        "m_ef" + to_string(efConstruction) + "m" + to_string(M) + ".bin";
    mapstr["index"] = hnsw_index;
    CheckDataset(using_dataset, mappmt, mapstr);

    if (stage == "build" || stage == "both")
        build_index<DTSET, DTRES>(mappmt, mapstr);

    if (stage == "search" || stage == "both")
        search_index<DTSET, DTRES>(mappmt, mapstr);

    return;
}
