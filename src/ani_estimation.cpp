/*
This file contains functions involving ANI estimation from set intersections
*/

#include "ani_estimator.hpp"

// Helper function to compute the containment based on the intersection size and the set size
double containment(int intersection, int set_size)
{
    return ((double)intersection) / ((double)set_size);
}

// Helper function to compute (containment)^(1/k) as an estimate of the ANI
double binomial_estimator(double containment, int kmer_num_ones)
{
    return std::pow(containment, (((double)1.0) / ((double)kmer_num_ones)));
}