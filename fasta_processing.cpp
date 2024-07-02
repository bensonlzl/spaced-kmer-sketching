#include "fasta_processing.hpp"
#include "logging.hpp"

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