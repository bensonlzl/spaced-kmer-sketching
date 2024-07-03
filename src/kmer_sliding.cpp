#include "kmer.hpp"
#include "fasta_processing.hpp"

// Helper function to update the kmer window
inline void update_kmer_window(kmer_bitset &current_kmer_window, const uint8_t &nucleotide_bits, const int &window_length)
{
    current_kmer_window <<= NUCLEOTIDE_BIT_SIZE;
    current_kmer_window[0] = (nucleotide_bits & 0x1);
    current_kmer_window[1] = ((nucleotide_bits & 0x2) >> 1);
}

inline void update_reversed_kmer_window(kmer_bitset &current_kmer_window, const uint8_t &nucleotide_bits, const int &window_length)
{
    current_kmer_window >>= NUCLEOTIDE_BIT_SIZE;
    current_kmer_window[NUCLEOTIDE_BIT_SIZE * window_length - 2] = (nucleotide_bits & 0x1);
    current_kmer_window[NUCLEOTIDE_BIT_SIZE * window_length - 1] = ((nucleotide_bits & 0x2) >> 1);
}

// Function to convert a string it into a list of kmers
// Assumes that the nucleotide string only contains acgt/ACGT
void nucleotide_string_to_kmers(
    std::vector<kmer> &kmer_list,
    const acgt_string &nucleotide_string,
    const kmer_bitset &mask,
    const int window_length,
    const std::function<bool(const kmer)> &sketching_cond)
{
    // If the string length is too short, no kmers in this string
    int nucleotide_string_length = nucleotide_string.size();
    if (nucleotide_string_length < window_length)
    {
        return;
    }

    // Initialise an empty kmer
    kmer_bitset current_kmer_window(KMER_BITSET_SIZE);

    // Create the first kmer window
    for (int idx = 0; idx + 1 < window_length; ++idx)
    {
        update_kmer_window(current_kmer_window, nucleotide_string[idx], window_length);
    }

    // Shift the window by one each time and add each new kmer
    for (int idx = 0; idx + window_length - 1 < nucleotide_string_length; ++idx)
    {
        update_kmer_window(current_kmer_window, nucleotide_string[idx + window_length - 1], window_length);
        kmer constructed_kmer = {
            window_length,
            current_kmer_window,
            mask,
            current_kmer_window & mask,
        };
        kmer canon_kmer = canonical_kmer(constructed_kmer);
        if (sketching_cond(canon_kmer))
            kmer_list.push_back(canon_kmer);
    }
}

void nucleotide_string_to_kmers_inbuilt_reverse(
    std::vector<kmer> &kmer_list,
    const acgt_string &nucleotide_string,
    const kmer_bitset &mask,
    const int window_length,
    const std::function<bool(const kmer)> &sketching_cond)
{
    int nucleotide_string_length = nucleotide_string.size();
    if (nucleotide_string_length < window_length)
    {
        return;
    }

    acgt_string complement_nucleotide_string = nucleotide_string;

    // Initialise an empty kmer
    kmer_bitset current_kmer_window(KMER_BITSET_SIZE), reversed_current_kmer_window(KMER_BITSET_SIZE);
    const kmer_bitset reversed_mask = (reverse_kmer_bitset(mask) >> ((MAX_KMER_LENGTH - window_length) * NUCLEOTIDE_BIT_SIZE));

    // Create the first kmer window
    for (int idx = 0; idx + 1 < window_length; ++idx)
    {
        update_kmer_window(current_kmer_window, nucleotide_string[idx], window_length);
        update_reversed_kmer_window(reversed_current_kmer_window, nucleotide_string[idx] ^ 0x3, window_length);
    }

    // Shift the window by one each time and add each new kmer
    for (int idx = 0; idx + window_length - 1 < nucleotide_string_length; ++idx)
    {
        int next_idx = idx + window_length - 1;
        update_kmer_window(current_kmer_window, nucleotide_string[next_idx], window_length);
        update_reversed_kmer_window(reversed_current_kmer_window, nucleotide_string[next_idx] ^ 0x3, window_length);
        if (DEBUG)
            std::cout << "Current kmer            :" << current_kmer_window << std::endl;
        if (DEBUG)
            std::cout << "Current complement kmer :" << reversed_current_kmer_window << std::endl;
        kmer_bitset masked_main_strand = current_kmer_window & mask, masked_complement_strand = reversed_current_kmer_window & reversed_mask;
        kmer canon_kmer = ((masked_main_strand < masked_complement_strand) ? kmer(window_length, current_kmer_window, mask, current_kmer_window & mask) : kmer(window_length, reversed_current_kmer_window, reversed_mask, reversed_current_kmer_window & reversed_mask));
        if (sketching_cond(canon_kmer))
            kmer_list.push_back(canon_kmer);
    }
}

// Helper functions to compute the kmers in a list of nucleotide strings
// This version computes by reference to avoid copying, appends the kmers to a given list of kmers
void nucleotide_string_list_to_kmers_by_reference(
    std::vector<kmer> &kmer_list,
    const std::vector<acgt_string> &nucleotide_strings,
    const kmer_bitset &mask,
    const int window_length,
    const std::function<bool(const kmer)> &sketching_cond)
{
    for (acgt_string const &s : nucleotide_strings)
    {
        nucleotide_string_to_kmers_inbuilt_reverse(kmer_list, s, mask, window_length, sketching_cond);
    }
}

// Helper functions to compute the kmers in a list of nucleotide strings
// This version explicity returns a vector of kmers
std::vector<kmer> nucleotide_string_list_to_kmers(
    const std::vector<acgt_string> &nucleotide_strings,
    const kmer_bitset &mask,
    const int window_length,
    const std::function<bool(const kmer)> &sketching_cond)
{
    std::vector<kmer> return_kmers;
    nucleotide_string_list_to_kmers_by_reference(
        return_kmers,
        nucleotide_strings,
        mask,
        window_length,
        sketching_cond);
    return return_kmers;
}
