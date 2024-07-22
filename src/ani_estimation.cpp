/**
 * @file ani_estimation.cpp
 * @author your name (you@domain.com)
 * @brief
 * @date 2024-07-04
 *
 * This file contains functions involving ANI estimation from set intersections
 *
 * @copyright Copyright (c) 2024
 */

#include "ani_estimator.hpp"

constexpr int ANI_DEBUG = DEBUG | 0;

/**
 * @brief
 * Helper function to compute the containment based on the intersection size and the set size
 *
 * @param intersection number of elements in the intersection
 * @param set_size size of the set
 * @return double
 */
double containment(int intersection, int set_size)
{
    if (intersection == 0) return 0;
    else return ((double)intersection) / ((double)set_size);
}

/**
 * @brief
 * Helper function to compute (containment)^(1/k) as an estimate of the ANI
 *
 * @param containment containment value
 * @param kmer_num_ones number of positions used in the kmer
 * @return double
 */
double binomial_estimator(double containment, int kmer_num_ones)
{
    if (containment <= 0) return 0;
    else return std::pow(containment, (((double)1.0) / ((double)kmer_num_ones)));
}