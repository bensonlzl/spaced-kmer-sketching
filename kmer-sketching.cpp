#include "kmer.hpp"
#include <chrono>

/*
This implementation supports kmers of length up to 128 as the bitsets used will have length 256

Support for spaced seeds is done using a mask of the same size. 
- Default contiguous kmers are implemented by setting the first 2*k bits of the mask to 1 and the rest to 0
*/


// Helper function to print a list of strings
void print_strings(const std::vector<std::string> &string_list){
    for (std::string s : string_list){
        std::cout << s << std::endl;
    }
}

const kmer_bitset nucleotide_to_bits_table[4] = {
    kmer_bitset(KMER_BITSET_SIZE,0),
    kmer_bitset(KMER_BITSET_SIZE,1),
    kmer_bitset(KMER_BITSET_SIZE,2),
    kmer_bitset(KMER_BITSET_SIZE,3),
};

// Inlined function to convert nucleotide letters to 2-bit words
inline uint8_t nucleotide_to_bits(const char nucleotide){
    uint8_t ret_val = 0;
    switch (nucleotide){
        case 'a':
            ret_val = 0;
            break;
        case 'A':
            ret_val = 0;
            break;
        case 'c':
            ret_val = 1;
            break;
        case 'C':
            ret_val = 1;
            break;
        case 'g':
            ret_val = 2;
            break;
        case 'G':
            ret_val = 2;
            break;
        case 't':
            ret_val = 3;
            break;
        case 'T':
            ret_val = 3;
            break;
        default:
            ret_val = 4;
            break;
    }
    return ret_val;
}

// Sketching function
frac_min_hash fmh;

inline bool sketching_condition(const kmer test_kmer){
    const int c = 200;
    return (fmh(test_kmer) % c == 0);
}

// Inlined function to return the associated kmer_bitset
inline kmer_bitset nucleotide_to_btiset(const char nucleotide){
    return nucleotide_to_bits_table[nucleotide_to_bits(nucleotide)];
}

// Function to read the input file and convert it into a list of strings
// Referenced from https://rosettacode.org/wiki/FASTA_format#C++
std::vector<std::string> strings_from_fasta(const char fasta_filename[]){
    std::ifstream fasta_file(fasta_filename);
    if (!fasta_file.good()){
        std::cerr << "Unable to open " << fasta_filename << ". \n Exiting..." << std::endl;
        exit(1);
    }

    std::vector<std::string> return_strings;

    std::string line, name, content;
    for (std::string line; std::getline(fasta_file, line);){
        if (line.empty() || line[0] == '>'){
            if (!name.empty()){
                if (LOGGING) std::clog << INFO_LOG << "Read " << name << " from file " << fasta_filename << std::endl;
                return_strings.push_back(content);
            }
            if (!line.empty()){
                name = line.substr(1);
            }
            content.clear();
        }
        else if (!name.empty()){
            if (line.find(' ') != std::string::npos){
                name.clear();
                content.clear();
            }
            else{
                content += line;
            }
        }
    }
    if (!name.empty()){
        if (LOGGING) std::clog << INFO_LOG << "Read " << name << " from file " << fasta_filename << std::endl;
        return_strings.push_back(content);
    }    

    return return_strings;
}


// Helper function to add nucleotide strings
void add_nucleotide_strings(std::vector<std::vector<uint8_t>> &return_strings, const std::string &raw_string){
    std::vector<uint8_t> cur_nucleotides;
    int raw_string_length = raw_string.length();
    for (int idx = 0; idx < raw_string_length; ++idx){
        uint8_t nucleotide_bits = nucleotide_to_bits(raw_string[idx]);
        if (nucleotide_bits & 0x4){
            if (!cur_nucleotides.empty()){
                return_strings.push_back(cur_nucleotides);
            }
            cur_nucleotides.clear();
        }
        else{
            cur_nucleotides.push_back(nucleotide_bits);
        }
    }
    if (!cur_nucleotides.empty()){
        return_strings.push_back(cur_nucleotides);
    }
}

// Function to split the strings at non-nucleotide characters
std::vector<std::vector<uint8_t>> cut_nucleotide_strings(const std::vector<std::string> &raw_strings){
    std::vector<std::vector<uint8_t>> return_strings;
    for (std::string const &raw_string : raw_strings){
        add_nucleotide_strings(return_strings,raw_string);
    }
    return return_strings;
}

// Helper function to update the kmer window
inline void update_kmer_window(kmer_bitset &current_kmer_window, const uint8_t &nucleotide_bits, const int &window_length){
    current_kmer_window <<= NUCLEOTIDE_BIT_SIZE;
    // current_kmer_window |= nucleotide_to_btiset(nucleotide);
    current_kmer_window[0] = (nucleotide_bits & 0x1);
    current_kmer_window[1] = ((nucleotide_bits & 0x2) >> 1);
    // current_kmer_window &= contiguous_kmer(window_length);
}

// Function to convert a string it into a list of kmers
// Assumes that the nucleotide string only contains acgt/ACGT
std::vector<kmer> nucleotide_string_to_kmers(
    const std::vector<uint8_t> &nucleotide_string, 
    const kmer_bitset mask, 
    const int window_length,
    std::function<bool(const kmer)> sketching_cond
){
    std::vector<kmer> kmer_list;

    // If the string length is too short, no kmers in this string
    int nucleotide_string_length = nucleotide_string.size();
    if (nucleotide_string_length < window_length){
        return kmer_list;
    }

    // Initialise an empty kmer
    kmer_bitset current_kmer_window(KMER_BITSET_SIZE);

    // Current rightmost index of the kmer
    int idx;

    // Create the first kmer window
    for (idx = 0; idx + 1 < window_length; ++idx){
        update_kmer_window(current_kmer_window,nucleotide_string[idx],window_length);
    }

    // Shift the window by one each time and add each new kmer
    for (idx = 0; idx + window_length - 1 < nucleotide_string_length; ++idx){
        update_kmer_window(current_kmer_window,nucleotide_string[idx+window_length-1],window_length);
        kmer constructed_kmer = {
            window_length,
            current_kmer_window, 
            mask,
            current_kmer_window & mask,
        };
        kmer canon_kmer = canonical_kmer(constructed_kmer);
        if (sketching_cond(canon_kmer)) kmer_list.push_back(canon_kmer);
    }
    
    return kmer_list;
}


std::vector<kmer> nucleotide_string_list_to_kmers(
    const std::vector<std::vector<uint8_t>> nucleotide_strings, 
    const kmer_bitset mask, 
    const int window_length,
    std::function<bool(const kmer)> sketching_cond
){
    std::vector<kmer> return_kmers;
    for (std::vector<uint8_t> const &s : nucleotide_strings){
        std::vector<kmer> constructed_kmers = nucleotide_string_to_kmers(s,mask,window_length,sketching_cond);
        return_kmers.insert(return_kmers.end(),constructed_kmers.begin(),constructed_kmers.end());
    }
    return return_kmers;
}


kmer_set kmer_set_from_file(const char filename[], const kmer_bitset mask, const int window_size, std::function<bool(const kmer)> sketching_cond){
    return kmer_set(
        nucleotide_string_list_to_kmers(
            cut_nucleotide_strings(strings_from_fasta(filename)),
            mask,
            window_size,
            sketching_cond
        )
    );
}


int main(int argc, char *argv[]){
    auto t0 = std::chrono::high_resolution_clock::now();

    initialise_contiguous_kmer_array();
    initialise_reversing_kmer_array();

    auto t1 = std::chrono::high_resolution_clock::now();

    const int kmer_size = 20;
    kmer_bitset mask = contiguous_kmer(kmer_size);

    if (LOGGING) std::clog << INFO_LOG << " kmer size = " << kmer_size << std::endl;
    if (LOGGING) std::clog << INFO_LOG << " kmer mask = " << mask << std::endl;

    std::vector<kmer_set> kmer_set_data;
    for (int i = 1; i < argc; ++i){
        kmer_set_data.push_back(
            kmer_set_from_file(argv[i],mask,kmer_size,sketching_condition)
        );
    }
    
    auto t2 = std::chrono::high_resolution_clock::now();

    int data_size = kmer_set_data.size();
    for (int i = 0; i < data_size; ++i){
        std::cout << "Comparing files " << i << " and " << i+1 << std::endl;
        std::cout << kmer_set_intersection(kmer_set_data[i],kmer_set_data[(i+1)%data_size]) << std::endl;
    }
    
    auto t3 = std::chrono::high_resolution_clock::now();

    std::cout << "Time taken for initialisation = " << std::chrono::duration<double,std::milli>(t1-t0).count() << " ms" << std::endl;
    std::cout << "Time taken for sketching = " << std::chrono::duration<double,std::milli>(t2-t1).count() << " ms" << std::endl;
    std::cout << "Time taken for set comparison = " << std::chrono::duration<double,std::milli>(t3-t2).count() << " ms" << std::endl;


}