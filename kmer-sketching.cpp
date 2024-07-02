#include "kmer.hpp"
#include "ani_estimator.hpp"
#include "fasta_processing.hpp"
#include <chrono>
#include <cilk/cilk.h>

// Helper function to print a list of strings
void print_strings(const std::vector<std::string> &string_list){
    for (std::string s : string_list){
        std::cout << s << std::endl;
    }
}

// Sketching function
frac_min_hash fmh;
inline bool sketching_condition(const kmer &test_kmer){
    const int c = 200;
    return (fmh(test_kmer) % c == 0);
}



int main(int argc, char *argv[]){
    auto t0 = std::chrono::high_resolution_clock::now();

    initialise_contiguous_kmer_array();
    initialise_reversing_kmer_array();

    auto t1 = std::chrono::high_resolution_clock::now();

    const int kmer_size = 20;
    kmer_bitset mask = contiguous_kmer(kmer_size);
    const int kmer_num_indices = (mask.count() / NUCLEOTIDE_BIT_SIZE);

    if (LOGGING) std::clog << INFO_LOG << " kmer size = " << kmer_size << std::endl;
    if (LOGGING) std::clog << INFO_LOG << " kmer mask = " << mask << std::endl;

    std::vector<kmer_set> kmer_set_data(argc-1);

    auto t_preprocess_string = std::chrono::high_resolution_clock::now();

    cilk_for (int i = 1; i < argc; ++i){
        kmer_set_data[i-1] = kmer_set_from_fasta_file(
            argv[i],
            mask,
            kmer_size,
            sketching_condition
        );
    }

    auto t_postprocess_kmers = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken for sketching = " << std::chrono::duration<double,std::milli>(t_postprocess_kmers-t_preprocess_string).count() << " ms" << std::endl;
    // auto t_preinsert_kmers = std::chrono::high_resolution_clock::now();


    // auto t_postinsert_kmers = std::chrono::high_resolution_clock::now();
    // std::cout << "Time taken for insertion = " << std::chrono::duration<double,std::milli>(t_postinsert_kmers-t_preinsert_kmers).count() << " ms" << std::endl;

    
    auto t2 = std::chrono::high_resolution_clock::now();

    int data_size = kmer_set_data.size();
    for (int i = 0; i < data_size; ++i){
        std::cout << "Comparing files " << i << " and " << i+1 << std::endl;
        int intersection = kmer_set_intersection(kmer_set_data[i],kmer_set_data[(i+1)%data_size]);
        double contain = containment(intersection,kmer_set_data[i].kmer_set_size());
        double ani_estimate = binomial_estimator(contain,kmer_num_indices);
        std::cout << "Intersection = " << intersection << "\nContainment = " << contain << "\nANI Estimate = " << ani_estimate << std::endl;
    }
    
    auto t3 = std::chrono::high_resolution_clock::now();

    std::cout << "Total time taken for initialisation = " << std::chrono::duration<double,std::milli>(t1-t0).count() << " ms" << std::endl;
    std::cout << "Total time taken for string processing and sketching = " << std::chrono::duration<double,std::milli>(t2-t1).count() << " ms" << std::endl;
    std::cout << "Total time taken for set comparison = " << std::chrono::duration<double,std::milli>(t3-t2).count() << " ms" << std::endl;


}