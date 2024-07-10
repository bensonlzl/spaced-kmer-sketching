/**
 * @file kmer_sliding.cpp
 * @author Benson Lin (bensonlinzl@gmail.com)
 * @brief
 * @date 2024-07-04
 *
 * @copyright Copyright (c) 2024
 * 
 * This file contains functions related to the sliding window algorithm that computes the kmers in a given string of nucleotides
 *
 */
#include "kmer.hpp"
#include "fasta_processing.hpp"

/**
 * @brief
 * Helper function to update the kmer window
 * Shifts the kmer window to the LEFT by the size of the nucleotide and updates the bits
 *
 * @param current_kmer_window A reference to the current kmer window
 * @param nucleotide_bits The nucleotide bits of the current nucleotide
 * @param window_length The length of the window
 */
inline void update_kmer_window(kmer_bitset &current_kmer_window, const uint8_t &nucleotide_bits, const int &window_length)
{
    current_kmer_window <<= NUCLEOTIDE_BIT_SIZE;
    current_kmer_window[0] = (nucleotide_bits & 0x1);
    current_kmer_window[1] = ((nucleotide_bits & 0x2) >> 1);
}

/**
 * @brief
 * Helper function to update the COMPLEMENTED kmer window
 * Shifts the kmer window to the RIGHT by the size of the nucleotide and updates the bits
 *
 * @param current_kmer_window A reference to the current complement kmer window
 * @param nucleotide_bits The nucleotide bits of the current nucleotide
 * @param window_length The length of the window
 */
inline void update_complement_kmer_window(kmer_bitset &current_kmer_window, const uint8_t &nucleotide_bits, const int &window_length)
{
    current_kmer_window >>= NUCLEOTIDE_BIT_SIZE;
    current_kmer_window[NUCLEOTIDE_BIT_SIZE * window_length - 2] = (nucleotide_bits & 0x1);
    current_kmer_window[NUCLEOTIDE_BIT_SIZE * window_length - 1] = ((nucleotide_bits & 0x2) >> 1);
}

/**
 * @brief
 * Function to convert an ACGT string into a list of CANONICAL kmers by reference
 * Uses a sliding window with a kmer_bitset representing the current window in order to efficiently compute the current canonical kmer
 * NOTE: Has been replaced by the updated function below since the latter is much faster
 *
 * @param kmer_list reference to a list of kmers for appending new kmers
 * @param nucleotide_string ACGT string representing the ACGT nucleotides
 * @param mask spaced seed mask used
 * @param window_length window size of the kmer
 * @param sketching_cond boolean function on kmers to decide which kmers are used
 */
void nucleotide_string_to_kmers_OLD_reverse(
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

/**
 * @brief
 * Function to convert an ACGT string into a list of kmers by reference
 * Uses a sliding window with a kmer_bitset representing the current window in order to efficiently compute the current canonical kmer
 * Implicitly constructs the complement strand of nucleotide_string in order to efficiently compute the reverse complement
 *
 * @param kmer_list reference to a list of kmers for appending new kmers
 * @param nucleotide_string ACGT string representing the ACGT nucleotides
 * @param mask spaced seed mask used
 * @param window_length window size of the kmer
 * @param sketching_cond boolean function on kmers to decide which kmers are used
 */
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

    // Initialise an empty kmer for both the main strand and the complement strand
    kmer_bitset current_kmer_window(KMER_BITSET_SIZE), reversed_current_kmer_window(KMER_BITSET_SIZE);

    // Initialise the reverse mask (is this necessary?)
    // const kmer_bitset reversed_mask = (reverse_kmer_bitset(mask) >> ((MAX_KMER_LENGTH - window_length) * NUCLEOTIDE_BIT_SIZE));

    // Create the first kmer window
    for (int idx = 0; idx + 1 < window_length; ++idx)
    {
        uint8_t cur_nucleotide = nucleotide_string[idx];
        update_kmer_window(current_kmer_window, cur_nucleotide, window_length);

        // nucleotide_string[idx] ^ 0x3 computes the complementary nucleotide
        update_complement_kmer_window(reversed_current_kmer_window, cur_nucleotide ^ 0x3, window_length);
    }

    // Shift the window by one each time and add each new kmer
    for (int idx = 0; idx + window_length - 1 < nucleotide_string_length; ++idx)
    {
        // Get the index of the next nucleotide character
        int next_idx = idx + window_length - 1;
        uint8_t cur_nucleotide = nucleotide_string[next_idx];

        update_kmer_window(current_kmer_window, cur_nucleotide, window_length);
        update_complement_kmer_window(reversed_current_kmer_window, cur_nucleotide ^ 0x3, window_length);

        if (DEBUG)
            std::cout << "Current kmer            :" << current_kmer_window << std::endl;
        if (DEBUG)
            std::cout << "Current complement kmer :" << reversed_current_kmer_window << std::endl;

        // Compute the masked kmers for both the main strand and the complement strand
        kmer_bitset masked_main_strand = current_kmer_window & mask;
        kmer_bitset masked_reverse_complement_strand = reversed_current_kmer_window & mask; // use the same mask on the reverse complement strand

        // Determine the canonical kmer by lexicographical comparison of the main strand and the reverse complement

        kmer_bitset *canonical_kmer_bits, *masked_canonical_kmer_bits;

        if (masked_main_strand < masked_reverse_complement_strand)
        {
            canonical_kmer_bits = &current_kmer_window;
            masked_canonical_kmer_bits = &masked_main_strand;
        }
        else
        {
            canonical_kmer_bits = &reversed_current_kmer_window;
            masked_canonical_kmer_bits = &masked_reverse_complement_strand;
        }

        // std::cout << "CANONICAL KMER " <<  *canonical_kmer_bits << '\n';
        // std::cout << "MASK " <<  mask << '\n';
        // std::cout << "MASKED CANONICAL KMER " <<  *masked_canonical_kmer_bits << '\n';
        // std::cout << std::endl;

        kmer canon_kmer(window_length, *canonical_kmer_bits, mask, *masked_canonical_kmer_bits);
        if (sketching_cond(canon_kmer))
            kmer_list.push_back(canon_kmer);
    }
}

/**
 * @brief
 * Helper function to compute the kmers in a list of nucleotide strings
 * This version computes by reference to avoid copying, appends the kmers to a given list of kmers
 *
 * @param kmer_list reference to a list of kmers for appending new kmers
 * @param nucleotide_strings list of ACGT strings representing the ACGT nucleotides
 * @param mask spaced seed mask used
 * @param window_length window size of the kmer
 * @param sketching_cond boolean function on kmers to decide which kmers are used
 */
void nucleotide_string_list_to_kmers_by_reference(
    std::vector<kmer> &kmer_list,
    const std::vector<acgt_string> &nucleotide_strings,
    const kmer_bitset &mask,
    const int window_length,
    const std::function<bool(const kmer)> &sketching_cond)
{
    for (acgt_string const &s : nucleotide_strings)
    {
        nucleotide_string_to_kmers(kmer_list, s, mask, window_length, sketching_cond);
    }
}

/**
 * @brief
 * Helper function to compute the kmers in a list of nucleotide strings
 * This version explicity returns a vector of kmers
 *
 * @param kmer_list reference to a list of kmers for appending new kmers
 * @param nucleotide_strings list of ACGT strings representing the ACGT nucleotides
 * @param mask spaced seed mask used
 * @param window_length window size of the kmer
 * @param sketching_cond boolean function on kmers to decide which kmers are used
 * @return returns a list of kmers in the nucleotide strings
 */
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
