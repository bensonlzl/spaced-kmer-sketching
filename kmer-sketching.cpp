#include "kmer.hpp"


/*
This implementation supports kmers of length up to 128 as the bitsets used will have length 256

Support for spaced seeds is done using a mask of the same size. 
- Default contiguous kmers are implemented by setting the first 2*k bits of the mask to 1 and the rest to 0
*/



// inline std::basic_string<char32_t> kmer_bitset_to_basic_string(const kmer_bitset &kb){
//     std::basic_string<char32_t> return_ints;
//     for (int num_shifts = 0; num_shifts < KMER_UINT32_SIZE; ++num_shifts){
//         return_ints.push_back(
//             (
//                 contiguous_kmer(32 / NUCLEOTIDE_BIT_SIZE) 
//                 & ((kb & contiguous_kmer((num_shifts+1) * 32 / NUCLEOTIDE_BIT_SIZE))  >> (num_shifts * 32))
//             ).to_ulong()
//         );
//     }
//     return return_ints;
// }

// bool kmer_bits_lexi_smaller(const kmer_bitset &kb1, const kmer_bitset &kb2){
//     std::basic_string<char32_t> kb_string_1 = kmer_bitset_to_basic_string(kb1);
//     std::basic_string<char32_t> kb_string_2 = kmer_bitset_to_basic_string(kb2);
//     for (int num_shifts = 0; num_shifts < KMER_UINT32_SIZE; ++num_shifts){
//         if (kb_string_1[KMER_UINT32_SIZE - num_shifts - 1] < kb_string_2[KMER_UINT32_SIZE - num_shifts - 1]) return true;
//         if (kb_string_1[KMER_UINT32_SIZE - num_shifts - 1] > kb_string_2[KMER_UINT32_SIZE - num_shifts - 1]) return false;
//     }
//     return false;
// }


// Helper function to print a list of strings
void print_strings(const std::vector<std::string> string_list){
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
void add_nucleotide_strings(std::vector<std::string> &return_strings, const std::string &raw_string){
    std::string cur_string;
    int raw_string_length = raw_string.length();
    for (int idx = 0; idx < raw_string_length; ++idx){
        if (nucleotide_to_bits(raw_string[idx]) & 0x4){
            if (!cur_string.empty()){
                return_strings.push_back(cur_string);
            }
            cur_string.clear();
        }
        else{
            cur_string.push_back(raw_string[idx]);
        }
    }
    if (!cur_string.empty()){
        return_strings.push_back(cur_string);
    }
}

// Function to split the strings at non-nucleotide characters
std::vector<std::string> cut_nucleotide_strings(const std::vector<std::string> &raw_strings){
    std::vector<std::string> return_strings;
    for (std::string raw_string : raw_strings){
        add_nucleotide_strings(return_strings,raw_string);
    }
    return return_strings;
}

// Helper function to update the kmer window
inline void update_kmer_window(kmer_bitset &current_kmer_window, const char &nucleotide, const int &window_length){
    current_kmer_window <<= NUCLEOTIDE_BIT_SIZE;
    uint8_t nucleotide_bits = nucleotide_to_bits(nucleotide);
    // current_kmer_window |= nucleotide_to_btiset(nucleotide);
    current_kmer_window[0] = (nucleotide_bits & 0x1);
    current_kmer_window[1] = ((nucleotide_bits & 0x2) >> 1);
    // current_kmer_window &= contiguous_kmer(window_length);
}

// Function to convert a string it into a list of kmers
// Assumes that the nucleotide string only contains acgt/ACGT
std::vector<kmer> nucleotide_string_to_kmers(const std::string &nucleotide_string, const kmer_bitset mask, const int window_length){
    std::vector<kmer> kmer_list;

    // If the string length is too short, no kmers in this string
    int nucleotide_string_length = nucleotide_string.length();
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
        kmer_list.push_back(canonical_kmer(constructed_kmer));
    }
    
    return kmer_list;
}


std::vector<kmer> nucleotide_string_list_to_kmers(const std::vector<std::string> nucleotide_strings, const kmer_bitset mask, const int window_length){
    std::vector<kmer> return_kmers;
    for (std::string s : nucleotide_strings){
        std::vector<kmer> constructed_kmers = nucleotide_string_to_kmers(s,mask,window_length);
        return_kmers.insert(return_kmers.end(),constructed_kmers.begin(),constructed_kmers.end());
    }
    return return_kmers;
}

// Sketching function
bool sketching_condition(const kmer test_kmer){
    return true; // TO DO
}

// Given a sketching function and list of kmers, 
// reduce the list of kmers down to a smaller list using the sketching condition
std::vector<kmer> sketch_filter(const std::vector<kmer> kmers, std::function<bool(const kmer)> sketching_cond){
    auto filtered_kmers = std::ranges::views::filter(kmers, sketching_cond);
    std::vector<kmer> return_kmers;
    for (kmer filtered_kmer : filtered_kmers){
        return_kmers.push_back(filtered_kmer);
    }
    return return_kmers;
}


int main(int argc, char *argv[]){

    if (LOGGING) std::clog << INFO_LOG << " kmer size = " << sizeof(kmer) << std::endl;
    initialise_contiguous_kmer_array();
    initialise_reversing_kmer_array();

    

    const int kmer_size = 60;
    kmer_bitset mask = contiguous_kmer(kmer_size);

    kmer_set kmer_set_1 = kmer_set(
        sketch_filter(
            nucleotide_string_list_to_kmers(
                cut_nucleotide_strings(strings_from_fasta(argv[1])),
                mask,
                kmer_size
            ),
            sketching_condition
        )
    );
    kmer_set kmer_set_2 = kmer_set(
        sketch_filter(
            nucleotide_string_list_to_kmers(
                cut_nucleotide_strings(strings_from_fasta(argv[2])),
                mask,
                kmer_size
            ),
            sketching_condition
        )
    );
    if (LOGGING) std::clog << "Comparing ks1 and ks2" << std::endl;
    std::cout << kmer_set_intersection(kmer_set_1,kmer_set_2) << std::endl;
    std::cout << kmer_set_1.kmer_set_size() << std::endl;
}