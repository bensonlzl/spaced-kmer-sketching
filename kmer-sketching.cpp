#include "kmer.hpp"
#include "ani_estimator.hpp"
#include "fasta_processing.hpp"
#include <chrono>

// Helper function to print a list of strings
void print_strings(const std::vector<std::string> &string_list){
    for (std::string s : string_list){
        std::cout << s << std::endl;
    }
}

// Sketching function, example
frac_min_hash fmh(0);
inline bool sketching_condition(const kmer &test_kmer){
    const int c = 200;
    return (fmh(test_kmer) % c == 0);
}

void test_compute_pairwise_ANI_estimation_contiguous_kmers(const int kmer_size, const int num_files, char *filenames[]){
    kmer_bitset mask = contiguous_kmer(kmer_size);
    const int kmer_num_indices = (mask.count() / NUCLEOTIDE_BIT_SIZE);

    if (LOGGING) std::clog << INFO_LOG << " kmer size = " << kmer_size << std::endl;
    if (LOGGING) std::clog << INFO_LOG << " kmer mask = " << mask << std::endl;
    
    auto t_preprocess_string = std::chrono::high_resolution_clock::now();

    std::vector<kmer_set> kmer_set_data = parallel_kmer_sets_from_fasta_files(
        num_files,
        filenames,
        mask,
        kmer_size,
        sketching_condition
    );
    auto t_postprocess_kmers = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken for sketching = " << std::chrono::duration<double,std::milli>(t_postprocess_kmers-t_preprocess_string).count() << " ms" << std::endl;
    int data_size = kmer_set_data.size();
    for (int i = 0; i < data_size; ++i){
        std::cout << "Comparing files " << i << " and " << i+1 << std::endl;
        int intersection = kmer_set_intersection(kmer_set_data[i],kmer_set_data[(i+1)%data_size]);
        double contain = containment(intersection,kmer_set_data[i].kmer_set_size());
        double ani_estimate = binomial_estimator(contain,kmer_num_indices);
        std::cout << "Intersection = " << intersection << "\nContainment = " << contain << "\nANI Estimate = " << ani_estimate << std::endl;
    }
}

int main(int argc, char *argv[]){
    initialise_contiguous_kmer_array();
    initialise_reversing_kmer_array();
    test_compute_pairwise_ANI_estimation_contiguous_kmers(20,argc-1,argv+1); // test on all files given in argv
}