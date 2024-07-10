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
#include <chrono>

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
frac_min_hash fmh(0);
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
 * Helper function to compute all pairwise combinations in a list (vector)
 * 
 * @tparam T 
 * @param v list of elements (vector)
 * @return a pair of vectors that whose corresponding pairs of elements represent every pair of elements in v
 */
template<typename T>
std::pair<std::vector<T>, std::vector<T>> generate_pairs_from_vector(
    const std::vector<T> &v
){
    std::vector<T> v1, v2;

    for (size_t i = 0; i < v.size(); ++i){
        for (size_t j = 0; j < v.size(); ++j){
            v1.push_back(v[i]);
            v2.push_back(v[j]);
        }
    }

    return make_pair(v1,v2);
}

/**
 * @brief 
 * Helper function to compute pairwise ANI estimation (adjacent pairs) with a random space seed
 * 
 * @param window_size total window size of the spaced seed
 * @param kmer_size number of characters used in the kmer
 * @param num_files number of files
 * @param filenames list of filenames (as a char* array)
 * @param output_filename name of the output file
 */
void test_compute_adjacent_pairwise_ANI_estimation_random_spaced_kmers(
    const int window_size,
    const int kmer_size,
    const int num_files,
    char *filenames[],
    const std::string &output_filename)
{
    // kmer_bitset mask = contiguous_kmer(kmer_size);
    kmer_bitset mask = generate_random_spaced_seed_mask(window_size, kmer_size);
    const int kmer_num_indices = (mask.count() / NUCLEOTIDE_BIT_SIZE); // How many nucleotides are in the kmer

    auto t_preprocess_string = std::chrono::high_resolution_clock::now();

    std::vector<kmer_set> kmer_set_data = parallel_kmer_sets_from_fasta_files(
        num_files,
        filenames,
        mask,
        kmer_size,
        sketching_condition);
    auto t_postprocess_kmers = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken for sketching = " << std::chrono::duration<double, std::milli>(t_postprocess_kmers - t_preprocess_string).count() << " ms" << std::endl;
    std::vector<kmer_set *> kmer_sets_1, kmer_sets_2;
    std::vector<std::string> kmer_filenames_1, kmer_filenames_2;
    int data_size = kmer_set_data.size();
    for (int i = 0; i < data_size; ++i)
    {
        kmer_sets_1.push_back(&kmer_set_data[i]);
        kmer_sets_2.push_back(&kmer_set_data[(i + 1) % data_size]);

        kmer_filenames_1.push_back(std::string(filenames[i]));
        kmer_filenames_2.push_back(std::string(filenames[(i + 1) % data_size]));
    }
    std::vector<int> intersection_vals = parallel_compute_pairwise_kmer_set_intersections(kmer_sets_1, kmer_sets_2);
    std::vector<double> containment_vals(data_size), ani_estimate_vals(data_size);
    for (int i = 0; i < data_size; ++i)
    {
        containment_vals[i] = containment(intersection_vals[i], kmer_sets_1[i]->kmer_set_size());
        ani_estimate_vals[i] = binomial_estimator(containment_vals[i], kmer_num_indices);
    }
    auto t_postcomparison = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken for comparison = " << std::chrono::duration<double, std::milli>(t_postcomparison - t_postprocess_kmers).count() << " ms" << std::endl;
    write_to_csv(
        kmer_filenames_1,
        kmer_filenames_2,
        ani_estimate_vals,
        window_size,
        mask,
        output_filename);
}



/**
 * @brief 
 * Helper function to compute all pairwise ANI estimation with a random space seed
 * 
 * @param window_size total window size of the spaced seed
 * @param kmer_size number of characters used in the kmer
 * @param num_files number of files
 * @param filenames list of filenames (as a char* array)
 * @param output_filename name of the output file
 */
void test_compute_all_pairwise_ANI_estimation_random_spaced_kmers(
    const int window_size,
    const int kmer_size,
    const int num_files,
    char *filenames[],
    const std::string &output_filename,
    bool is_append = false)
{
    // kmer_bitset mask = contiguous_kmer(kmer_size);
    kmer_bitset mask = generate_random_spaced_seed_mask(window_size, kmer_size);
    const int kmer_num_indices = (mask.count() / NUCLEOTIDE_BIT_SIZE); // How many nucleotides are in the kmer

    auto t_preprocess_string = std::chrono::high_resolution_clock::now();

    std::vector<kmer_set> kmer_set_data = parallel_kmer_sets_from_fasta_files(
        num_files,
        filenames,
        mask,
        kmer_size,
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

    auto kmer_sets_pairwise = generate_pairs_from_vector(kmer_sets_init);
    auto kmer_filenames_pairwise = generate_pairs_from_vector(kmer_filenames_init);

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

    if (DEBUG)
    {
        for (int i = 0; i < data_size; ++i)
        {
            std::cout << "INTERSECTION = " << intersection_vals[i] << " | CONTAINMENT  = " << containment_vals[i] << " | ANI ESTIMATE = " << ani_estimate_vals[i] << std::endl;
        }
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
    test_compute_all_pairwise_ANI_estimation_random_spaced_kmers(32, 32, argc - 2, argv + 2,filename,false); // test on all files given in argv
    for (int k = 10; k <= 10; ++k){
        test_compute_all_pairwise_ANI_estimation_random_spaced_kmers(k, k, argc - 2, argv + 2,filename,true); // test on all files given in argv
    }
    for (int k = 10; k <= 10; ++k){
        test_compute_all_pairwise_ANI_estimation_random_spaced_kmers(k+10, k, argc - 2, argv + 2,filename,true); // test on all files given in argv
    }
}