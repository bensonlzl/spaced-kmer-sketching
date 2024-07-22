/**
 * @file kmer-sketching.cpp
 * @author Benson Lin (bensonlinzl@gmail.com)
 * @brief
 * @date 2024-07-04
 *
 * @copyright Copyright (c) 2024
 */

#include "kmer.hpp"
#include "ani_estimator.hpp"
#include "fasta_processing.hpp"
#include "generators.hpp"

/**
 * @brief
 *
 * @param string_list
 */
void print_strings(const std::vector<std::string> &string_list)
{
    for (std::string s : string_list)
    {
        std::cout << s << std::endl;
    }
}

// Sketching function, example
frac_min_hash fmh(1);
inline bool sketching_condition(const kmer &test_kmer)
{
    const int c = 200;
    return (fmh(test_kmer) % c == 0);
}

/**
 * @brief 
 * Helper function that writes a list of estimated ANI values with the corresponding filenames and masks to a csv
 * 
 * @param filenames1 first list of filenames
 * @param filenames2 second list of filenames
 * @param estimated_values list of estimated ANI values
 * @param mask mask used 
 * @param output_filename name of the output csv 
 */
void write_to_csv(
    const std::vector<std::string> &filenames1,
    const std::vector<std::string> &filenames2,
    const std::vector<double> &estimated_values,
    const int window_size,
    const kmer_bitset &mask,
    const std::string &output_filename,
    bool is_append = false)
{
    // Open the output file
    std::ofstream output_file;
    if (is_append)
        output_file.open(output_filename, std::ios_base::app);
    else
        output_file.open(output_filename);

    if (!output_file.is_open())
    {
        std::cerr << "Error: Unable to open file " << output_filename << " for writing." << std::endl;
        return;
    }

    // Write header
    if (!is_append)
        output_file << "File 1,File 2,Estimated Value,Window Size,Mask" << std::endl;

    // Write data
    size_t numEntries = std::min(std::min(filenames1.size(), filenames2.size()), estimated_values.size());
    for (size_t i = 0; i < numEntries; ++i)
    {
        output_file << filenames1[i] << "," << filenames2[i] << "," << estimated_values[i] << "," << window_size << "," << mask << std::endl;
    }

    // Close the file
    output_file.close();
}

/**
 * @brief 
 * Wrapper for computing adjacent pairwise pairs of kmer_set pointers
 * 
 * @param kmer_set_pointers 
 * @return std::pair<std::vector<kmer_set *>,std::vector<kmer_set *>> 
 */
std::pair<std::vector<kmer_set *>,std::vector<kmer_set *>> compute_kmer_set_pointer_pairwise(
    std::vector<kmer_set *> kmer_set_pointers
){
    return generate_pairwise_from_vector<kmer_set *>(kmer_set_pointers);
}

/**
 * @brief 
 * Wrapper for computing all pairs of kmer_set pointers
 * 
 * @param kmer_set_pointers 
 * @return std::pair<std::vector<kmer_set *>,std::vector<kmer_set *>> 
 */
std::pair<std::vector<kmer_set *>,std::vector<kmer_set *>>  compute_kmer_set_pointer_all_pairs(
    std::vector<kmer_set *> kmer_set_pointers
){
    return generate_all_pairs_from_vector<kmer_set *>(kmer_set_pointers);
}

/**
 * @brief 
 * Wrapper for computing adjacent pairwise pairs of strings
 * 
 * @param kmer_set_pointers 
 * @return std::pair<std::vector<std::string>,std::vector<std::string>> 
 */
std::pair<std::vector<std::string>,std::vector<std::string>> compute_strings_pairwise(
    std::vector<std::string> kmer_set_pointers
){
    return generate_pairwise_from_vector<std::string>(kmer_set_pointers);
}

/**
 * @brief 
 * Wrapper for computing all pairs of strings
 * 
 * @param kmer_set_pointers 
 * @return std::pair<std::vector<std::string>,std::vector<std::string>> 
 */
std::pair<std::vector<std::string>,std::vector<std::string>> compute_strings_all_pairs(
    std::vector<std::string> kmer_set_pointers
){
    return generate_all_pairs_from_vector<std::string>(kmer_set_pointers);
}

/**
 * @brief 
 * Given a number of FASTA files, this function computes an ANI estimate and writes it to a .csv file
 * 
 * 
 * @tparam kmer_set_callable 
 * @tparam string_callable 
 * @param compute_kmer_set_pairs function to compute the kmer_set* pairs (must be compatible with compute_string_pairs)
 * @param compute_string_pairs function to compute the string pairs (must be compatible with compute_kmer_set_pairs)
 * @param window_size size of the kmer window
 * @param kmer_size number of characters to be used in the kmer
 * @param num_files number of FASTA files to be procesed
 * @param filenames array of char* representing the FASTA filenames
 * @param output_filename .csv file name
 * @param is_append bool to determine whether to
 */
template <typename kmer_set_callable, typename string_callable>
void test_compute_ANI_estimation_random_spaced_kmers(
    kmer_set_callable compute_kmer_set_pairs,
    string_callable compute_string_pairs,
    const int window_size,
    const int kmer_size,
    const int num_files,
    char *filenames[],
    const std::string &output_filename,
    bool is_append
){
    // kmer_bitset mask = contiguous_kmer(kmer_size);
    kmer_bitset mask = generate_random_spaced_seed_mask(window_size, kmer_size);
    const int kmer_num_indices = (mask.count() / NUCLEOTIDE_BIT_SIZE); // How many nucleotides are in the kmer

    auto t_preprocess_string = std::chrono::high_resolution_clock::now();

    std::vector<kmer_set> kmer_set_data = parallel_kmer_sets_from_fasta_files(
        num_files,
        filenames,
        mask,
        window_size,
        sketching_condition);
    auto t_postprocess_kmers = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken for sketching = " << std::chrono::duration<double, std::milli>(t_postprocess_kmers - t_preprocess_string).count() << " ms" << std::endl;
    std::vector<kmer_set *> kmer_sets_init;
    std::vector<std::string> kmer_filenames_init;

    for (int i = 0; i < kmer_set_data.size(); ++i)
    {
        kmer_sets_init.push_back(&kmer_set_data[i]);
        kmer_filenames_init.push_back(std::string(filenames[i]));
    }

    auto kmer_sets_pairwise = compute_kmer_set_pairs(kmer_sets_init);
    auto kmer_filenames_pairwise = compute_string_pairs(kmer_filenames_init);

    std::vector<int> intersection_vals = parallel_compute_pairwise_kmer_set_intersections(
        kmer_sets_pairwise.first, 
        kmer_sets_pairwise.second
    );

    int data_size = intersection_vals.size();

    std::vector<double> containment_vals(data_size), ani_estimate_vals(data_size);
    for (int i = 0; i < data_size; ++i)
    {
        containment_vals[i] = containment(intersection_vals[i], kmer_sets_pairwise.first[i]->kmer_set_size());
        ani_estimate_vals[i] = binomial_estimator(containment_vals[i], kmer_num_indices);
    }

    auto t_postcomparison = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken for comparison = " << std::chrono::duration<double, std::milli>(t_postcomparison - t_postprocess_kmers).count() << " ms" << std::endl;
    write_to_csv(
        kmer_filenames_pairwise.first,
        kmer_filenames_pairwise.second,
        ani_estimate_vals,
        window_size,
        mask,
        output_filename,
        is_append);
}

int main(int argc, char *argv[])
{
    initialise_contiguous_kmer_array();
    initialise_reversing_kmer_array();
    const std::string filename = std::string(argv[1]);//"../../data_temp/Single-Family-Cross-Genus.csv";
    test_compute_ANI_estimation_random_spaced_kmers(
        compute_kmer_set_pointer_all_pairs,
        compute_strings_all_pairs,
        10, 
        10, 
        argc - 2, 
        argv + 2,
        filename,false
    ); // test on all files given in argv
    for (int k = 11; k <= 40; ++k){
        test_compute_ANI_estimation_random_spaced_kmers(
            compute_kmer_set_pointer_all_pairs,
            compute_strings_all_pairs,
            k, k, argc - 2, argv + 2,filename,true); // test on all files given in argv
    }
    for (int k = 10; k <= 40; ++k){
        test_compute_ANI_estimation_random_spaced_kmers(
            compute_kmer_set_pointer_all_pairs,
            compute_strings_all_pairs,
            k+10, k, argc - 2, argv + 2,filename,true); // test on all files given in argv
    }
}