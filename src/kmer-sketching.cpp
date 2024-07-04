#include "kmer.hpp"
#include "ani_estimator.hpp"
#include "fasta_processing.hpp"
#include <chrono>

// Helper function to print a list of strings
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

void write_to_csv(
    const std::vector<std::string> &filenames1, 
    const std::vector<std::string> &filenames2, 
    const std::vector<double> &estimated_values, 
    const kmer_bitset &mask,
    const std::string &output_filename
) {
    // Open the output file
    std::ofstream output_file(output_filename);
    if (!output_file.is_open()) {
        std::cerr << "Error: Unable to open file " << output_filename << " for writing." << std::endl;
        return;
    }
    
    // Write header
    output_file << "File 1,File 2,Estimated Value,Mask" << std::endl;
    
    // Write data
    size_t numEntries = std::min(std::min(filenames1.size(), filenames2.size()), estimated_values.size());
    for (size_t i = 0; i < numEntries; ++i) {
        output_file << filenames1[i] << "," << filenames2[i] << "," << estimated_values[i] << "," << mask << std::endl;
    }
    
    // Close the file
    output_file.close();
}

void test_compute_pairwise_ANI_estimation_contiguous_kmers(const int window_size, const int kmer_size, const int num_files, char *filenames[])
{
    // kmer_bitset mask = contiguous_kmer(kmer_size);
    kmer_bitset mask = generate_random_spaced_seed_mask(window_size,kmer_size);
    const int kmer_num_indices = (mask.count() / NUCLEOTIDE_BIT_SIZE); // How many nucleotides are in the kmer

    if (LOGGING)
        std::clog << INFO_LOG << " kmer size = " << kmer_size << std::endl;
    if (LOGGING)
        std::clog << INFO_LOG << " kmer mask = " << mask << std::endl;

    auto t_preprocess_string = std::chrono::high_resolution_clock::now();

    std::vector<kmer_set> kmer_set_data = parallel_kmer_sets_from_fasta_files(
        num_files,
        filenames,
        mask,
        kmer_size,
        sketching_condition);
    auto t_postprocess_kmers = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken for sketching = " << std::chrono::duration<double, std::milli>(t_postprocess_kmers - t_preprocess_string).count() << " ms" << std::endl;

    std::vector<kmer_set*> kmer_sets_1, kmer_sets_2;
    std::vector<std::string> kmer_filenames_1, kmer_filenames_2;
    

    int data_size = kmer_set_data.size();

    for (int i = 0; i < data_size; ++i){
        kmer_sets_1.push_back(&kmer_set_data[i]);
        kmer_sets_2.push_back(&kmer_set_data[(i+1)%data_size]);

        kmer_filenames_1.push_back(std::string(filenames[i]));
        kmer_filenames_2.push_back(std::string(filenames[(i+1)%data_size]));
    }

    std::vector<int> intersection_vals = parallel_compute_pairwise_kmer_set_intersections(kmer_sets_1,kmer_sets_2);
    std::vector<double> containment_vals(data_size), ani_estimate_vals(data_size);

    for (int i = 0; i < data_size; ++i)
    {
        // std::cout << "Comparing files " << kmer_filenames_1[i] << " and " << kmer_filenames_2[i] << std::endl;
        containment_vals[i] = containment(intersection_vals[i], kmer_sets_1[i]->kmer_set_size());
        ani_estimate_vals[i] = binomial_estimator(containment_vals[i], kmer_num_indices);
        // std::cout << "Intersection = " << intersection_vals[i] << "\nContainment = " << containment_vals[i] << "\nANI Estimate = " << ani_estimate_vals[i] << std::endl;
    }

    auto t_postcomparison = std::chrono::high_resolution_clock::now();
    std::cout << "Time taken for comparison = " << std::chrono::duration<double, std::milli>(t_postcomparison - t_postprocess_kmers).count() << " ms" << std::endl;

    write_to_csv(
        kmer_filenames_1,
        kmer_filenames_2,
        ani_estimate_vals,
        mask,
        "test_spaced.csv"
    );
}



int main(int argc, char *argv[])
{
    initialise_contiguous_kmer_array();
    initialise_reversing_kmer_array();
    test_compute_pairwise_ANI_estimation_contiguous_kmers(30,20, argc - 1, argv + 1); // test on all files given in argv
}