/**
 * @file kmer_bitset.cpp
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2024-07-04
 *
 * @copyright Copyright (c) 2024
 *
 */
#include "kmer.hpp"

/**
 * @brief
 * To save time in initialising bitsets for common operations involving a prefix of bits,
 * contiguous_kmer_array is initialised to store kmer_bitset with exactly 2l 1s
 * for each length 0 <= l <= MAX_KMER_LENGTH/
 */
static kmer_bitset contiguous_kmer_array[MAX_KMER_LENGTH + 1];

/**
 * @brief
 * Initialisation function to create the prefixes in contiguous_kmer_array
 * TO DO : Refactor this to be performed at compile time? Probably not possible since boost::dynamic_bitset is not a literal type
 */
void initialise_contiguous_kmer_array()
{
    for (int i = 0; i <= MAX_KMER_LENGTH; ++i)
    {
        contiguous_kmer_array[i] = kmer_bitset(KMER_BITSET_SIZE); // Initialise the bitset with the right size
    }
    for (int i = 1; i <= MAX_KMER_LENGTH; ++i)
    {
        contiguous_kmer_array[i] |= contiguous_kmer_array[i - 1]; // Take the previous prefix
        for (int j = 0; j < NUCLEOTIDE_BIT_SIZE; ++j)
        {
            contiguous_kmer_array[i][j + (i - 1) * NUCLEOTIDE_BIT_SIZE] = 1; // Add more 1s
        }
    }
}

/**
 * @brief
 * Helper function to obtain the correct contiguous kmer given the length
 *
 * @param kmer_length length of the required kmer
 * @return kmer_bitset of the correct length
 */
kmer_bitset contiguous_kmer(const int kmer_length)
{
    if (kmer_length > MAX_KMER_LENGTH)
        throw std::runtime_error("Given k-mer length exceeds maximum k-mer length");
    return contiguous_kmer_array[kmer_length];
}

/**
 * @brief
 * To reverse a kmer in LOG_KMER_BITSET_SIZE operations, we will use some bit hacks
 *
 * For each power of two from 2 to KMER_BITSET_SIZE / 2,
 * compute a bitset consisting of alternating 1s and 0s with run length equal to that power of two.
 */
kmer_bitset reversing_kmer_array[LOG_KMER_BITSET_SIZE];
kmer_bitset invert_reversing_kmer_array[LOG_KMER_BITSET_SIZE];

/**
 * @brief
 * Initialisation function to create the bitsets needed for reversing kmers
 * TO DO : Refactor this to be performed at compile time? Probably not possible since boost::dynamic_bitset is not a literal type
 */
void initialise_reversing_kmer_array()
{
    // First initiaise the bitsets to be of the right size
    for (int i = 0; i < LOG_KMER_BITSET_SIZE; ++i)
    {
        reversing_kmer_array[i] = kmer_bitset(KMER_BITSET_SIZE);
        invert_reversing_kmer_array[i] = kmer_bitset(KMER_BITSET_SIZE);
    }

    // Now we set the 1 bits in each bitset
    int array_idx = 0;
    for (int gap_size = NUCLEOTIDE_BIT_SIZE; gap_size < KMER_BITSET_SIZE; gap_size *= 2, ++array_idx)
    { // Double the gap size on each iteration
        for (uint32_t gap_idx = 0; gap_idx < KMER_BITSET_SIZE / gap_size; ++gap_idx)
        { // Iterate over each block of size (gap_size)
            for (int pos_idx = gap_idx * gap_size; pos_idx < (gap_idx + 1) * gap_size; ++pos_idx)
            {
                reversing_kmer_array[array_idx][pos_idx] = (gap_idx & 0x1);                                  // Alternate between blocks of 0s and 1s
                invert_reversing_kmer_array[array_idx][pos_idx] = ~reversing_kmer_array[array_idx][pos_idx]; // Alternate between blocks of 1s and 0s
            }
        }
    }
}

/**
 * @brief
 * Function to reverse a kmer_bitset (used for reverse complementation)
 * Creates a new kmer_bitset with the reversed bits
 *
 * @param kbs kmer_bitset to be reversed
 * @return reversed kmer_bitset
 */
kmer_bitset reverse_kmer_bitset(const kmer_bitset &kbs)
{
    int array_idx = 0;
    kmer_bitset cur = kbs;

    // log KMER_BITSET_SIZE iterations of reversing in blocks of powers of 2
    for (int gap_size = NUCLEOTIDE_BIT_SIZE; gap_size < KMER_BITSET_SIZE; gap_size *= 2, ++array_idx)
    {
        // Get one half of the bits by AND-ing with reversing_kmer_array and
        // the other half of the bis by AND-ing with invert_reversing_kmer_array,
        // then swap them by performing appropriate bitshifts
        cur = (((cur & reversing_kmer_array[array_idx]) >> gap_size) | ((cur & (invert_reversing_kmer_array[array_idx])) << gap_size));
    }
    return cur;
}

/**
 * @brief
 * Helper function that generates a random spaced seed mask across a window of size window_size
 * with exactly kmer_size characters that are used (rest are ignored)
 * Also accepts an optional random seed to randomise the seed
 *
 * @param window_size total span of the spaced seed
 * @param kmer_size number of characters to be used in the spaced seed
 * @param random_seed random seed to generate a different spaced seed, default is 0 as defined in kmer.hpp
 * @return randomly generated kmer_bitset as a mask
 */
kmer_bitset generate_random_spaced_seed_mask(
    const int window_size,
    const int kmer_size,
    size_t random_seed)
{
    kmer_bitset mask(KMER_BITSET_SIZE);

    std::vector<int> bit_indices(window_size);
    std::iota(bit_indices.begin(), bit_indices.end(), 0);
    std::shuffle(bit_indices.begin(), bit_indices.end(), std::mt19937(random_seed));

    for (int i = 0; i < kmer_size; ++i)
    {
        for (int j = 0; j < NUCLEOTIDE_BIT_SIZE; ++j)
        {
            mask[NUCLEOTIDE_BIT_SIZE * bit_indices[i] + j] = 1;
        }
    }

    return mask;
}
