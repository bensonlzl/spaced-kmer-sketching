#include "kmer.hpp"

// Hash table for the kmers

// Compute the number of elements in the intersection
int kmer_set_intersection(const kmer_set &ks1, const kmer_set &ks2) {
    int inters = 0;
    for (auto it : ks2.kmer_hashes){
        if (ks1.kmer_hashes.find(it.first) != ks1.kmer_hashes.end()){
            inters++;
        }
    }
    return inters;
}

