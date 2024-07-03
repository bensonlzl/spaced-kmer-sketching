/*
This file contains a number of functions related to processing the  .fasta files used to store the genomic data
*/

#include "fasta_processing.hpp"

#include "logging.hpp"

/***
 * Inlined function to convert nucleotide letters to 2-bit words
 * - A -> 0
 * - C -> 1
 * - G -> 2
 * - T -> 3
 * - Anything else -> 4
 *
 * The above values are chosen so that
 * 1) Complementation can be done with a single flip operation
 * 2) Lexicographical ordering is the same as if we used the original ACGT
 * 3) We can detect non-ACGT by checking the 3rd lowest bit
 *
 * @param nucleotide The nucleotide character to be read
 * @return 0,1,2,3,4 depending on the nucleotide character
 */
inline uint8_t nucleotide_to_bits(const char nucleotide)
{
    uint8_t ret_val = 0;
    switch (nucleotide)
    {
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

/***
 * Helper function to read a fasta file and convert it into a list of strings
 * Referenced from https://rosettacode.org/wiki/FASTA_format#C++
 *
 * @param fasta_filename Filename of the fasta file
 * @return A list of strings representing the nucleotide data in the file
 */
std::vector<std::string> strings_from_fasta(const char fasta_filename[])
{

    // Create the input file stream
    std::ifstream fasta_file(fasta_filename);

    // If unable to open the file, print and error and exit
    if (!fasta_file.good())
    {
        std::cerr << "Unable to open " << fasta_filename << ". \n Exiting..." << std::endl;
        exit(1);
    }

    // Vector of strings to be returned
    std::vector<std::string> return_strings;

    std::string line, name, content;
    for (std::string line; std::getline(fasta_file, line);)
    {
        if (line.empty() || line[0] == '>')
        {
            if (!name.empty())
            {
                if (LOGGING)
                    std::clog << INFO_LOG << "Read " << name << " from file " << fasta_filename << std::endl;
                return_strings.push_back(content);
            }
            if (!line.empty())
            {
                name = line.substr(1);
            }
            content.clear();
        }
        else if (!name.empty())
        {
            if (line.find(' ') != std::string::npos)
            {
                name.clear();
                content.clear();
            }
            else
            {
                content += line;
            }
        }
    }
    if (!name.empty())
    {
        if (LOGGING)
            std::clog << INFO_LOG << "Read " << name << " from file " << fasta_filename << std::endl;
        return_strings.push_back(content);
    }

    return return_strings;
}

/***
 * Helper function to add nucleotide strings
 * Takes in a reference to a vector of ACGT strings and the current string to be processed
 * Mutates the vector of ACGT strings by appending new ACGT strings
 *
 * @param return_strings List of ACGT strings to be mutated
 * @param raw_string String of nucleotide characters to be processed
 */
void add_nucleotide_strings(
    std::vector<acgt_string> &return_strings,
    const std::string &raw_string)
{
    // ACGT string storing the current nucleotides
    acgt_string cur_nucleotides;

    // Get the raw string length
    int raw_string_length = raw_string.length();

    // For each index of the raw string, we convert that
    for (int idx = 0; idx < raw_string_length; ++idx)
    {
        uint8_t nucleotide_bits = nucleotide_to_bits(raw_string[idx]);

        // If the nucleotide character is not ACGT, cut the ACGT string here
        if (nucleotide_bits & 0x4)
        {
            if (!cur_nucleotides.empty())
            {
                return_strings.push_back(cur_nucleotides);
            }
            cur_nucleotides.clear();
        }
        // Otherwise append the new nucleotide to the current vector
        else
        {
            cur_nucleotides.push_back(nucleotide_bits);
        }
    }
    // If there is still an ACGT string left, add this string as well
    if (!cur_nucleotides.empty())
    {
        return_strings.push_back(cur_nucleotides);
    }
}

/*



*/
/***
 * Helper function to split a list of strings at non-nucleotide characters\
 * Directly calls add_nucleotide_strings over every string in the given list of strings
 * and returns the list of ACGT strings generated
 *
 * @param raw_strings List of strings of nucleotide characters to be processed
 * @return List of ACGT strings in the strings
 */
std::vector<acgt_string> cut_nucleotide_strings(const std::vector<std::string> &raw_strings)
{
    std::vector<acgt_string> return_strings;
    for (std::string const &raw_string : raw_strings)
    {
        add_nucleotide_strings(return_strings, raw_string);
    }
    return return_strings;
}

/***
 * Helper function to compute the nucleotide strings from a fasta file given the name of the file
 * Composes the above helper functions together for convenience
 *
 * @param fasta_filename Filename of the fasta file
 * @return List of ACGT strings in file
 */
std::vector<acgt_string> nucleotide_strings_from_fasta_file(const char fasta_filename[])
{
    return cut_nucleotide_strings(strings_from_fasta(fasta_filename));
}