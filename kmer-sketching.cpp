#include "xxhash.hpp"
#include <bits/extc++.h>

#include <functional>
#include <cstdint>
#include <vector>
#include <bitset>
#include <cstring>
#include <fstream>
#include <iostream>
#include <ranges>




/*
This implementation supports kmers of length up to 128 as the bitsets used will have length 256

Support for spaced seeds is done using a mask of the same size. 
- Default contiguous kmers are implemented by setting the first 2*k bits of the mask to 1 and the rest to 0
*/

#define NUCLEOTIDE_SIZE 2
#define MAX_KMER_LENGTH 128 
#define KMER_BITSET_SIZE NUCLEOTIDE_SIZE * MAX_KMER_LENGTH // Must be divisible by 64!
#define KMER_CHAR32_SIZE KMER_BITSET_SIZE / 32

#define INFO_LOG "[INFO] "
#define DEBUG 0

typedef std::bitset<KMER_BITSET_SIZE> kmer_bitset;



kmer_bitset contiguous_kmer_array[MAX_KMER_LENGTH+1];

// TO DO : Replace this initialisation with a compile time initialisation
// kmer reverse_complement currently relies on calling this function
void initialise_contiguous_kmer_bitsets(){
    std::clog << INFO_LOG << "Initialising contiguous_kmer_array" << std::endl;
    for (int i = 1; i <= MAX_KMER_LENGTH; ++i){
        contiguous_kmer_array[i] |= contiguous_kmer_array[i-1];
        for (int j = 0; j < NUCLEOTIDE_SIZE; ++j){
            contiguous_kmer_array[i][j + (i-1) * NUCLEOTIDE_SIZE] = 1;
        }
    }
}

inline kmer_bitset contiguous_kmer(const int kmer_length){
    if (kmer_length > MAX_KMER_LENGTH) throw std::runtime_error("Given k-mer length exceeds maximum k-mer length");
    return contiguous_kmer_array[kmer_length];
}

inline std::basic_string<char32_t> kmer_bitset_to_basic_string(const kmer_bitset &kb){
    std::basic_string<char32_t> return_ints;
    for (int num_shifts = 0; num_shifts < KMER_CHAR32_SIZE; ++num_shifts){
        return_ints.push_back(
            (
                contiguous_kmer(32 / NUCLEOTIDE_SIZE) 
                & ((kb & contiguous_kmer((num_shifts+1) * 32 / NUCLEOTIDE_SIZE))  >> (num_shifts * 32))
            ).to_ulong()
        );
    }
    return return_ints;
}

bool kmer_bits_lexi_smaller(const kmer_bitset &kb1, const kmer_bitset &kb2){
    std::basic_string<char32_t> kb_string_1 = kmer_bitset_to_basic_string(kb1);
    std::basic_string<char32_t> kb_string_2 = kmer_bitset_to_basic_string(kb2);
    for (int num_shifts = 0; num_shifts < KMER_CHAR32_SIZE; ++num_shifts){
        if (kb_string_1[KMER_CHAR32_SIZE - num_shifts - 1] < kb_string_2[KMER_CHAR32_SIZE - num_shifts - 1]) return true;
        if (kb_string_1[KMER_CHAR32_SIZE - num_shifts - 1] > kb_string_2[KMER_CHAR32_SIZE - num_shifts - 1]) return false;
    }
    return false;
}

// Custom struct to represent a kmer
struct kmer{
    int window_length;          // Length of the whole kmer_window
    kmer_bitset kmer_bits;          // Raw bits in the kmer
    kmer_bitset mask;               // Mask used for the kmer
    kmer_bitset masked_bits;        // Masked bits of the kmer

    bool operator==(const kmer &other)const{
        return (masked_bits == other.masked_bits) && (mask == other.mask);
    }

    kmer reverse_complement(){
        kmer_bitset rc_bits = kmer_bits ^ contiguous_kmer(window_length);
        return {
            window_length,
            rc_bits,
            mask,
            rc_bits & mask
        };
    }

    kmer canonical_kmer(){
        kmer rc = this->reverse_complement();
        return (kmer_bits_lexi_smaller(masked_bits,rc.masked_bits) ? *this : rc);
    }
};



// Custom hash to hash the kmer, based on xxHash
struct kmer_hash {
    uint64_t operator()(const kmer &k) const {
        uint64_t k_hash = (
            xxh::xxhash3<64,char32_t>(kmer_bitset_to_basic_string(k.kmer_bits)) ^
            xxh::xxhash3<64,char32_t>(kmer_bitset_to_basic_string(k.mask)) ^
            xxh::xxhash3<64>((void*) &k.window_length, sizeof(int), (uint64_t) 0)
        );
        if (DEBUG) std::cout << "HASHED VALUE: " << k_hash << '\n';
        return k_hash;
    }
};

// Hash table for the kmers
typedef __gnu_pbds::gp_hash_table<kmer,int,kmer_hash> kmer_hash_table;

// Custom struct to store the kmers
struct kmer_set{
    kmer_hash_table kmer_hashes;

    kmer_set(const std::vector<kmer> kmers) {
        if (DEBUG) std::cout << "INSERTION BEGIN"  << '\n';

        for (kmer k : kmers){
            if (DEBUG) std::cout << "INSERTING KMER " << k.masked_bits << '\n';
            if (DEBUG) std::cout << "LENGTH " << KMER_BITSET_SIZE * sizeof(char) << " " << k.kmer_bits.to_string().data()  << '\n';
            kmer_hashes[k] = 1;

        }
        if (DEBUG) std::cout << "INSERTATION COMPLETED" << '\n';
    
    }

    inline int size(){
        return kmer_hashes.size();
    }

    // Compute the number of elements in the intersection
    int intersection(const kmer_set &other) const {
        int inters = 0;
        for (auto it : other.kmer_hashes){
            if (kmer_hashes.find(it.first) != kmer_hashes.end()){
                inters++;
            }
        }
        return inters;
    }
};





// Helper function to print a list of strings
void print_strings(const std::vector<std::string> string_list){
    for (std::string s : string_list){
        std::cout << s << std::endl;
    }
}


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
                std::clog << INFO_LOG << "Read " << name << " from file " << fasta_filename << std::endl;
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
        std::clog << INFO_LOG << "Read " << name << " from file " << fasta_filename << std::endl;
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
    current_kmer_window <<= NUCLEOTIDE_SIZE;
    current_kmer_window |= nucleotide_to_bits(nucleotide);
    current_kmer_window &= contiguous_kmer(window_length);
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
    kmer_bitset current_kmer_window;
    current_kmer_window.reset();

    // Current rightmost index of the kmer
    int idx;

    // Create the first kmer window
    for (idx = 0; idx + 1 < window_length; ++idx){
        update_kmer_window(current_kmer_window,nucleotide_string[idx],window_length);
    }

    // Shift the window by one each time and add each new kmer
    for (idx = 0; idx + window_length < nucleotide_string_length; ++idx){
        update_kmer_window(current_kmer_window,nucleotide_string[idx+window_length-1],window_length);
        kmer constructed_kmer = {
            window_length,
            current_kmer_window, 
            mask,
            current_kmer_window & mask,
        };
        kmer_list.push_back(constructed_kmer.canonical_kmer());
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
    auto filtered_kmers = std::ranges::filter_view(kmers, sketching_cond);
    std::vector<kmer> return_kmers;
    for (kmer filtered_kmer : filtered_kmers){
        return_kmers.push_back(filtered_kmer);
    }
    return return_kmers;
}

// Compute the intersection of two kmer sets
int intersection(const kmer_set ks1, const kmer_set ks2){
    return ks1.intersection(ks2);
}

int main(int argc, char *argv[]){
    initialise_contiguous_kmer_bitsets();

    std::clog << INFO_LOG << " kmer size = " << sizeof(kmer) << std::endl;
    std::clog << INFO_LOG << "KMER_INT_SIZE = " << KMER_CHAR32_SIZE << std::endl;

    const int kmer_size = 6;
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
    kmer_set kmer_set_3 = kmer_set(
        sketch_filter(
            nucleotide_string_list_to_kmers(
                cut_nucleotide_strings(strings_from_fasta(argv[1])),
                mask,
                kmer_size
            ),
            sketching_condition
        )
    );
    std::cout << intersection(kmer_set_1,kmer_set_2) << std::endl;
    std::cout << intersection(kmer_set_1,kmer_set_3) << std::endl;
    std::cout << kmer_set_1.size() << std::endl;
}