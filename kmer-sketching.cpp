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

// Sketching function
frac_min_hash fmh;
inline bool sketching_condition(const kmer &test_kmer){
    const int c = 200;
    return (fmh(test_kmer) % c == 0);
}


// Helper function to update the kmer window
inline void update_kmer_window(kmer_bitset &current_kmer_window, const uint8_t &nucleotide_bits, const int &window_length){
    current_kmer_window <<= NUCLEOTIDE_BIT_SIZE;
    current_kmer_window[0] = (nucleotide_bits & 0x1);
    current_kmer_window[1] = ((nucleotide_bits & 0x2) >> 1);
}

// Function to convert a string it into a list of kmers
// Assumes that the nucleotide string only contains acgt/ACGT
void nucleotide_string_to_kmers(
    std::vector<kmer> &kmer_list,
    const std::vector<uint8_t> &nucleotide_string, 
    const kmer_bitset &mask, 
    const int &window_length,
    const std::function<bool(const kmer)> &sketching_cond
){
    // If the string length is too short, no kmers in this string
    int nucleotide_string_length = nucleotide_string.size();
    if (nucleotide_string_length < window_length){
        return;
    }

    // Initialise an empty kmer
    kmer_bitset current_kmer_window(KMER_BITSET_SIZE);

    // Current rightmost index of the kmer
    int idx;

    // Create the first kmer window
    for (idx = 0; idx + 1 < window_length; ++idx){
        update_kmer_window(current_kmer_window,nucleotide_string[idx],window_length);
    }

    // Shift the window by one each time and add each new kmer
    for (idx = 0; idx + window_length - 1 < nucleotide_string_length; ++idx){
        update_kmer_window(current_kmer_window,nucleotide_string[idx+window_length-1],window_length);
        kmer constructed_kmer = {
            window_length,
            current_kmer_window, 
            mask,
            current_kmer_window & mask,
        };
        kmer canon_kmer = canonical_kmer(constructed_kmer);
        if (sketching_cond(canon_kmer)) kmer_list.push_back(canon_kmer);
    }
}


std::vector<kmer> nucleotide_string_list_to_kmers(
    const std::vector<std::vector<uint8_t>> &nucleotide_strings, 
    const kmer_bitset &mask, 
    const int &window_length,
    const std::function<bool(const kmer)> &sketching_cond
){
    std::vector<kmer> return_kmers;
    for (std::vector<uint8_t> const &s : nucleotide_strings){
        nucleotide_string_to_kmers(return_kmers,s,mask,window_length,sketching_cond);
    }
    return return_kmers;
}

void nucleotide_string_list_to_kmers_by_reference(
    std::vector<kmer> &kmer_list,
    const std::vector<std::vector<uint8_t>> &nucleotide_strings, 
    const kmer_bitset &mask, 
    const int &window_length,
    const std::function<bool(const kmer)> &sketching_cond
){
    for (std::vector<uint8_t> const &s : nucleotide_strings){
        nucleotide_string_to_kmers(kmer_list,s,mask,window_length,sketching_cond);
    }
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

    std::vector<std::vector<std::vector<uint8_t>> > data_strings(argc-1);
    std::vector<std::vector<kmer> > kmer_lists(argc-1);
    std::vector<kmer_set> kmer_set_data(argc-1);

    auto t_preprocess_string = std::chrono::high_resolution_clock::now();

    for (int i = 1; i < argc; ++i){
        data_strings[i-1] = cut_nucleotide_strings(strings_from_fasta(argv[i]));
    }

    auto t_postprocess_string = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken for string processing = " << std::chrono::duration<double,std::milli>(t_postprocess_string-t_preprocess_string).count() << " ms" << std::endl;
    auto t_preprocess_kmers = std::chrono::high_resolution_clock::now();
    
    for (int i = 1; i < argc; ++i){
        nucleotide_string_list_to_kmers_by_reference(
            kmer_lists[i-1],
            data_strings[i-1],
            mask,
            kmer_size,
            sketching_condition
        );
    }

    auto t_postprocess_kmers = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken for sketching = " << std::chrono::duration<double,std::milli>(t_postprocess_kmers-t_preprocess_kmers).count() << " ms" << std::endl;
    auto t_preinsert_kmers = std::chrono::high_resolution_clock::now();


    for (int i = 1; i < argc; ++i){
        kmer_set_data[i-1].insert_kmers(
            kmer_lists[i-1]
        );
    }

    auto t_postinsert_kmers = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken for insertion = " << std::chrono::duration<double,std::milli>(t_postinsert_kmers-t_preinsert_kmers).count() << " ms" << std::endl;

    
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