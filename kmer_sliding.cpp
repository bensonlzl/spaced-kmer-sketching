#include "kmer.hpp"

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
    const int window_length,
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



void nucleotide_string_list_to_kmers_by_reference(
    std::vector<kmer> &kmer_list,
    const std::vector<std::vector<uint8_t>> &nucleotide_strings, 
    const kmer_bitset &mask, 
    const int window_length,
    const std::function<bool(const kmer)> &sketching_cond
){
    for (std::vector<uint8_t> const &s : nucleotide_strings){
        nucleotide_string_to_kmers(kmer_list,s,mask,window_length,sketching_cond);
    }
}

std::vector<kmer> nucleotide_string_list_to_kmers(
    const std::vector<std::vector<uint8_t>> &nucleotide_strings, 
    const kmer_bitset &mask, 
    const int window_length,
    const std::function<bool(const kmer)> &sketching_cond
){
    std::vector<kmer> return_kmers;
    nucleotide_string_list_to_kmers_by_reference(
        return_kmers,
        nucleotide_strings,
        mask,
        window_length,
        sketching_cond
    );
    return return_kmers;
}
