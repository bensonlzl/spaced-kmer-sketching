// STL includes
#include <bitset>
#include <cstring>
#include <stdexcept>
#include <cstdint>
#include <cstdlib>
#include <algorithm>
#include <unordered_map>
#include <iterator>
#include <vector>
#include <iostream>
#include <fstream>
#include <ranges>
#include <functional>
#include <random>

// Boost includes
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>

// OpenCilk needed for parallel processing of files
#include <cilk/cilk.h>

#include "logging.hpp"

// We will be using the boost dynamic bitsets for representing the kmers
typedef boost::dynamic_bitset<> kmer_bitset;

// ONLY MODIFY THIS
// const int LOG_KMER_UINT32_SIZE = 2;
constexpr int LOG_KMER_BITSET_SIZE = 6; // Make this as small as is necessary, increasing it makes the program slower
// Setting it to 6 gives 64-bit bitsets -> 32-mers, 7 -> 64-mers, 8 -> 128-mers, 9 -> 256-mers, 10 -> 512-mers

// DO NOT MODIFY THESE
constexpr int NUCLEOTIDE_BIT_SIZE = 2;
constexpr int KMER_BITSET_SIZE = (1 << LOG_KMER_BITSET_SIZE);
constexpr int MAX_KMER_LENGTH = (KMER_BITSET_SIZE / NUCLEOTIDE_BIT_SIZE);

// Functions for initialising contiguous kmers and the kmer reversing functions
void initialise_contiguous_kmer_array();
kmer_bitset contiguous_kmer(const int kmer_length);
void initialise_reversing_kmer_array();
kmer_bitset reverse_kmer_bitset(const kmer_bitset &kbs);
kmer_bitset generate_random_spaced_seed_mask(
    const int window_size,
    const int kmer_size,
    size_t random_seed = 0
);

// Struct to store information about the kmer
struct kmer
{
    int window_length;       // Length of the whole kmer_window
    kmer_bitset kmer_bits;   // Raw bits in the kmer
    kmer_bitset mask;        // Mask used for the kmer
    kmer_bitset masked_bits; // Masked bits of the kmer

    bool operator==(const kmer &other) const
    {
        return (masked_bits == other.masked_bits) && (mask == other.mask);
    }
};

// Functions for computing canonical kmers
kmer reverse_complement(kmer k);
kmer canonical_kmer(kmer k);

// Helper functions to compute a list of kmers from a list of nucleotide strings
std::vector<kmer> nucleotide_string_list_to_kmers(
    const std::vector<std::vector<uint8_t>> &nucleotide_strings,
    const kmer_bitset &mask,
    const int window_length,
    const std::function<bool(const kmer)> &sketching_cond);
void nucleotide_string_list_to_kmers_by_reference(
    std::vector<kmer> &kmer_list,
    const std::vector<std::vector<uint8_t>> &nucleotide_strings,
    const kmer_bitset &mask,
    const int window_length,
    const std::function<bool(const kmer)> &sketching_cond);

// Struct for computing kmer hashes using std::hash
// This is primary used in the hash map to check if an identical kmer exists
struct kmer_hash
{
    std::hash<kmer_bitset> kmer_bitset_std_hash;
    std::hash<int> int_std_hash;
    inline size_t operator()(const kmer &k) const
    {
        size_t k_hash = (kmer_bitset_std_hash(k.masked_bits) ^
                         kmer_bitset_std_hash(k.mask) ^
                         int_std_hash(k.window_length));
        return k_hash;
    }
};

// Struct for FracMinHash, initialised with a nonce to create a different hash
// This is primarily used in FracMinHash to determine which kmers are kept in the sketching process
struct frac_min_hash
{
    boost::hash<kmer_bitset> kmer_bitset_boost_hash;
    boost::hash<int> int_boost_hash;
    int nonce;

    frac_min_hash(int n) : nonce(int_boost_hash(n)) {}

    // Hash function
    inline size_t operator()(const kmer &k) const
    {
        size_t k_hash = (kmer_bitset_boost_hash(k.masked_bits) ^ kmer_bitset_boost_hash(k.mask) ^ int_boost_hash(k.window_length) ^ nonce);
        return k_hash;
    }
};

// hash table for kmers
typedef std::unordered_map<kmer, int, kmer_hash> kmer_hash_table;

// Custom struct to store the kmers in a set
struct kmer_set
{
    kmer_hash_table kmer_hashes;

    // Helper function for inserting kmers
    void insert_kmers(const std::vector<kmer> &kmers)
    {
        for (kmer k : kmers)
        {
            kmer_hashes[k] = 1;
        }
    }

    // Helper function for computing the size of the kmer set
    inline int kmer_set_size() const
    {
        return kmer_hashes.size();
    }
};
// Helper function for computing kmer set intersection
int kmer_set_intersection(const kmer_set &ks1, const kmer_set &ks2);

// Helper functions to compute kmer sets from fasta files
kmer_set kmer_set_from_fasta_file(
    const char fasta_filename[],
    const kmer_bitset &mask,
    const int window_length,
    const std::function<bool(const kmer)> &sketching_cond);
std::vector<kmer_set> kmer_sets_from_fasta_files(
    int num_files,
    char *fasta_filenames[],
    const kmer_bitset &mask,
    const int window_length,
    const std::function<bool(const kmer)> &sketching_cond);
// This version processes the fasta files in parallel using OpenCilk
std::vector<kmer_set> parallel_kmer_sets_from_fasta_files(
    int num_files,
    char *fasta_filenames[],
    const kmer_bitset &mask,
    const int window_length,
    const std::function<bool(const kmer)> &sketching_cond);
std::vector<int> compute_pairwise_kmer_set_intersections(
    const std::vector<kmer_set*> &kmer_sets_1,
    const std::vector<kmer_set*> &kmer_sets_2
);
std::vector<int> parallel_compute_pairwise_kmer_set_intersections(
    const std::vector<kmer_set*> &kmer_sets_1,
    const std::vector<kmer_set*> &kmer_sets_2
);