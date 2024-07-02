#include "kmer.hpp"
#include "fasta_processing.hpp"

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

kmer_set kmer_set_from_fasta_file(
    const char fasta_filename[],
    const kmer_bitset &mask,
    const int kmer_size,
    const std::function<bool(const kmer)> &sketching_cond
){
    kmer_set ks;
    ks.insert_kmers(
        nucleotide_string_list_to_kmers(
            nucleotide_strings_from_fasta_file(fasta_filename),
            mask,
            kmer_size,
            sketching_cond
        )
    );
    return ks;
}

std::vector<kmer_set> kmer_sets_from_fasta_files(
    int num_files,
    char *fasta_filenames[],
    const kmer_bitset &mask,
    const int kmer_size,
    const std::function<bool(const kmer)> &sketching_cond
){
    std::vector<kmer_set> kmer_sets(num_files);
    for (int i = 0; i < num_files; ++i){
        kmer_sets[i] = kmer_set_from_fasta_file(
            fasta_filenames[i],
            mask,
            kmer_size,
            sketching_cond
        );
    }
    return kmer_sets;
}

std::vector<kmer_set> parallel_kmer_sets_from_fasta_files(
    int num_files,
    char *fasta_filenames[],
    const kmer_bitset &mask,
    const int kmer_size,
    const std::function<bool(const kmer)> &sketching_cond
){
    std::vector<kmer_set> kmer_sets(num_files);
    cilk_for (int i = 0; i < num_files; ++i){
        kmer_sets[i] = kmer_set_from_fasta_file(
            fasta_filenames[i],
            mask,
            kmer_size,
            sketching_cond
        );
    }
    return kmer_sets;
}