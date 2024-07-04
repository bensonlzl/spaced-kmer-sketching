#include "kmer.hpp"
#include "fasta_processing.hpp"


// Compute the number of elements in the intersection
/***
 * Helper function to compute the number of kmers in the intersection of two kmer sets
 * 
 * @param ks1 kmer_set containing the first set of kmers
 * @param ks2 kmer_set cotnaining the second set of kmers
 * @return number of kmers in the intersection of the two kmer sets
 */
int kmer_set_intersection(const kmer_set &ks1, const kmer_set &ks2)
{
    // Minor optimization: Since we are iterating over ks2, we swap the kmer_sets to ensure that ks2 is the smaller of the two sets
    if (ks1.kmer_set_size() < ks2.kmer_set_size()) return kmer_set_intersection(ks2,ks1);

    // Counter of the number of intersections
    int inters = 0;

    // Iterate over ks2 and check if the kmer is in ks1
    for (auto it : ks2.kmer_hashes) 
    {
        if (ks1.kmer_hashes.find(it.first) != ks1.kmer_hashes.end())
        {
            inters++;
        }
    }
    return inters;
}

/***
 * Helper function that generates a kmer_set from a .fasta file
 * Composes a number of other helper functions together
 * 
 * @param fasta_filename path to the .fasta file to be read
 * @param mask spaced seed mask used 
 * @param window_length window size of the kmer
 * @param sketching_cond boolean function on kmers to decide which kmers are used
 * @return kmer_set containing the kmers sketched from that file 
 */
kmer_set kmer_set_from_fasta_file(
    const char fasta_filename[],
    const kmer_bitset &mask,
    const int window_length,
    const std::function<bool(const kmer)> &sketching_cond)
{
    kmer_set ks;
    ks.insert_kmers(
        nucleotide_string_list_to_kmers(
            nucleotide_strings_from_fasta_file(fasta_filename),
            mask,
            window_length,
            sketching_cond));
    return ks;
}

/***
 * Iterates over a list of filenames and creates a kmer_set for each file
 * 
 * @param num_files number of files to be processed
 * @param fasta_filenames pointer to the list of file names to be read
 * @param mask spaced seed mask used 
 * @param window_length window size of the kmer
 * @param sketching_cond boolean function on kmers to decide which kmers are used 
 * @return a list of kmer_sets corresponding to the file names given
 */
std::vector<kmer_set> kmer_sets_from_fasta_files(
    int num_files,
    char *fasta_filenames[],
    const kmer_bitset &mask,
    const int window_length,
    const std::function<bool(const kmer)> &sketching_cond)
{
    std::vector<kmer_set> kmer_sets(num_files);
    for (int i = 0; i < num_files; ++i)
    {
        kmer_sets[i] = kmer_set_from_fasta_file(
            fasta_filenames[i],
            mask,
            window_length,
            sketching_cond);
    }
    return kmer_sets;
}

/***
 * Parallel version of kmer_sets_from_fasta_files
 * Uses a cilk_for to parallelize the for loop over the fasta files
 * 
 * @param num_files number of files to be processed
 * @param fasta_filenames pointer to the list of filenames to be read
 * @param mask spaced seed mask used 
 * @param kmer_size window size of the kmer
 * @param sketching_cond boolean function on kmers to decide which kmers are used
 * @return a list of kmer_sets corresponding to the file names given
 */
std::vector<kmer_set> parallel_kmer_sets_from_fasta_files(
    int num_files,
    char *fasta_filenames[],
    const kmer_bitset &mask,
    const int window_length,
    const std::function<bool(const kmer)> &sketching_cond)
{
    // Debug: If the PARALLEL_FILES flag is set to 0, use the serial version
    if (!PARALLEL_FILES) return kmer_sets_from_fasta_files(num_files,fasta_filenames,mask,window_length,sketching_cond);

    std::vector<kmer_set> kmer_sets(num_files);
    cilk_for(int i = 0; i < num_files; ++i)
    {
        kmer_sets[i] = kmer_set_from_fasta_file(
            fasta_filenames[i],
            mask,
            window_length,
            sketching_cond);
    }
    return kmer_sets;
}