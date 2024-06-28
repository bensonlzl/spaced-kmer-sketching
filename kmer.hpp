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
#include <boost/dynamic_bitset.hpp>

// #include "xxhash.hpp"
// #include <ext/pb_ds/assoc_container.h>

typedef boost::dynamic_bitset<> kmer_bitset;

// ONLY MODIFY THIS
// const int LOG_KMER_UINT32_SIZE = 2;
const int LOG_KMER_BITSET_SIZE = 8;

// DO NOT MODIFY THESE
const int NUCLEOTIDE_BIT_SIZE = 2;
// const int KMER_UINT32_SIZE = (1 << LOG_KMER_UINT32_SIZE);
// const int KMER_BITSET_SIZE = (KMER_UINT32_SIZE * 32) ;
// const int MAX_KMER_LENGTH = (KMER_BITSET_SIZE / NUCLEOTIDE_BIT_SIZE);
// const int LOG_KMER_BITSET_SIZE = (LOG_KMER_UINT32_SIZE + 5);

const int KMER_BITSET_SIZE = (1 << LOG_KMER_BITSET_SIZE) ;
const int MAX_KMER_LENGTH = (KMER_BITSET_SIZE / NUCLEOTIDE_BIT_SIZE);

const std::string INFO_LOG = "[INFO] ";
const int LOGGING = 1;
const int DEBUG = 0;

void initialise_contiguous_kmer_array();
kmer_bitset contiguous_kmer(const int kmer_length);
void initialise_reversing_kmer_array();
kmer_bitset reverse_kmer_bitset(kmer_bitset kbs);
// std::basic_string<char32_t> kmer_bitset_to_basic_string(const kmer_bitset &kbs);
// bool kmer_bitset_lexi_smaller(const kmer_bitset &kb1, const kmer_bitset &kb2);
// inline void mask_kmer_bitset(kmer_bitset &kbs, const kmer_bitset &mask);

// Struct to store information about the kmer 
struct kmer{
    int window_length;          // Length of the whole kmer_window
    kmer_bitset kmer_bits;          // Raw bits in the kmer
    kmer_bitset mask;               // Mask used for the kmer
    kmer_bitset masked_bits;        // Masked bits of the kmer

    bool operator==(const kmer &other)const{
        return (masked_bits == other.masked_bits) && (mask == other.mask);
    }
};
kmer reverse_complement(kmer k);
kmer canonical_kmer(kmer k);

struct kmer_hash {
    std::hash<kmer_bitset> kmer_bitset_std_hash;
    std::hash<int> int_std_hash;
    inline uint64_t operator()(const kmer &k) const {
        uint64_t k_hash = (
            kmer_bitset_std_hash(k.masked_bits) ^
            kmer_bitset_std_hash(k.mask) ^
            int_std_hash(k.window_length)
        );
        return k_hash;
    }
};

typedef std::unordered_map<kmer,int,kmer_hash> kmer_hash_table;

// Custom struct to store the kmers
struct kmer_set{
    kmer_hash_table kmer_hashes;

    kmer_set(const std::vector<kmer> kmers) {
        if (DEBUG) std::cout << "INSERTION BEGIN"  << '\n';
        for (kmer k : kmers){
            if (DEBUG) std::cout << "INSERTING KMER " << k.masked_bits << '\n';
            kmer_hashes[k] = 1;
        }
        if (DEBUG) std::cout << "INSERTION COMPLETED" << '\n';
    }

    inline int kmer_set_size() const {
        return kmer_hashes.size();
    }
};
int kmer_set_intersection(const kmer_set &ks1, const kmer_set &ks2);
double kmer_set_containment(const kmer_set &ks1, const kmer_set &ks2);