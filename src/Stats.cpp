/**
 * Metagenomics Canopy Clustering Implementation
 *
 * Copyright (C) 2013, 2014 Piotr Dworzynski (piotr@cbs.dtu.dk), Technical University of Denmark
 *
 * This file is part of Metagenomics Canopy Clustering Implementation.
 *
 * Metagenomics Canopy Clustering Implementation is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Metagenomics Canopy Clustering Implementation is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <iostream>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <limits>
#include <utility>
#include <vector>
#include <algorithm>

#include "Stats.hpp"

using namespace std;

void precompute_pearson_data(int sample_data_length, const PRECISIONT* __restrict__ sample_data, PRECISIONT* __restrict__ precomputed_pearson_data){

    //Calculate sum and average of data samples
    PRECISIONT sum = 0, avg = 0;
    for(int i = 0; i < sample_data_length; i++)
        sum += sample_data[i];

    avg =  sum / sample_data_length;

    //Calculate standard deviation of data samples
    PRECISIONT factor_sum = 0;
    for(int i = 0; i < sample_data_length; i++)
        factor_sum += pow((sample_data[i] - avg),2);

    PRECISIONT stddev = 0;
    stddev = sqrt(factor_sum/sample_data_length);

    //Precompute pearson data
    for(int i = 0; i < sample_data_length; i++){
        if(fabs(stddev) < 2* std::numeric_limits<PRECISIONT>::min())
            precomputed_pearson_data[i] = 0;
        else
            precomputed_pearson_data[i] = (sample_data[i] - avg)/(stddev * sample_data_length);
    }
}

//Spearman is basically pearson on ranked values
void precompute_spearman_data(int sample_data_length, const PRECISIONT* __restrict__ sample_data, PRECISIONT* __restrict__ precomputed_pearson_data){ 

    //First rank sample data
    vector<pair<PRECISIONT, int> > samples_and_indexes_vector;
    for(int i = 0; i < sample_data_length; i++)
        samples_and_indexes_vector.push_back(make_pair(sample_data[i], i));

    sort(samples_and_indexes_vector.begin(), samples_and_indexes_vector.end());

    PRECISIONT* ranked_sample_data = new PRECISIONT[sample_data_length];

    PRECISIONT rank = 1.0;
    int i = 0;
    while(i < sample_data_length){
        int j = i;

        while((j < sample_data_length-1) & (fabs(samples_and_indexes_vector[j].first - samples_and_indexes_vector[j+1].first) < 2* std::numeric_limits<PRECISIONT>::min()))
            j+=1;

        int num_equal = j-i+1;

        for(int j=0; j<num_equal; j++)
            ranked_sample_data[samples_and_indexes_vector[i+j].second] = rank + (num_equal-1)/2.0;

        rank += num_equal;
        i += num_equal;
    }

    //Then precompute pearson correlation values on top of the rank data vector
    precompute_pearson_data(sample_data_length, ranked_sample_data, precomputed_pearson_data);

    delete[] ranked_sample_data;
}

PRECISIONT corr_from_precomputed(int n, const PRECISIONT* __restrict__  v1, const PRECISIONT* __restrict__  v2){
    PRECISIONT sum = 0;
    for(int i = 0; i < n; i++){
        sum += v1[i] * v2[i];
    }
    return sum*n;
}


