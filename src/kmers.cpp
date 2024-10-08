/**
 * @file kmers.cpp
 * @author your name (you@domain.com)
 * @brief
 * @date 2024-07-04
 *
 * @copyright Copyright (c) 2024
 *
 */
#include "kmer.hpp"

constexpr int KMERS_DEBUG = DEBUG | 0;

// TO DO: Support spaced seeds
// Currently only supports palindromic masks
kmer reverse_complement(kmer k)
{
    kmer_bitset rc_bits = (reverse_kmer_bitset(k.kmer_bits).flip()) >> ((MAX_KMER_LENGTH - k.window_length) * NUCLEOTIDE_BIT_SIZE);
    if (KMERS_DEBUG)
        std::cout << k.kmer_bits << " reverse complemented to " << rc_bits << std::endl;
    if (KMERS_DEBUG)
        std::cout << "Masked bits " << (k.kmer_bits & k.mask) << " reverse complemented to " << (rc_bits & k.mask) << std::endl;
    return {
        k.window_length,
        rc_bits,
        k.mask,
        rc_bits & k.mask};
}

// Computes the canonical kmer for a kmer k
kmer canonical_kmer(kmer k)
{
    kmer rc = reverse_complement(k);
    return ((k.masked_bits < rc.masked_bits) ? k : rc);
}