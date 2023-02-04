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
#include <fstream>
#include <ctime>
#include <random>

#include <algorithm>

#include <omp.h>

#include "Log.hpp"

#include "signal_handlers.hpp"
#include "prog_bar_misc.hpp"
#include "Canopy.hpp"
#include "Point.hpp"
#include "CanopyClustering.hpp"

Canopy* CanopyClusteringAlg::create_canopy(Point* origin, vector<Point*>& points, vector<Point*>& close_points, PRECISIONT max_neighbour_dist, PRECISIONT max_close_dist, bool set_close_points){

    std::vector<Point*> neighbours;

    if(set_close_points){
        Point* potential_neighbour; 

        //Go through all points and set the close points to contain the ones that are "close"
        close_points.clear();//Will not reallocate
        for(int i=0; i<points.size(); i++){

            potential_neighbour = points[i];

            PRECISIONT dist = get_distance_between_points(origin, potential_neighbour);

            if(dist < max_close_dist){

                close_points.push_back(potential_neighbour);

                if(dist < max_neighbour_dist){

                    neighbours.push_back(potential_neighbour);

                }
            } 
        }

    } else {
            
        Point* potential_neighbour;

        for(int i=0; i<close_points.size(); i++){

            potential_neighbour = close_points[i];

            PRECISIONT dist = get_distance_between_points(origin, potential_neighbour);

            if(dist < max_neighbour_dist){
                neighbours.push_back(potential_neighbour);
            }
        }

    }

    if(neighbours.size()){
        return new Canopy(neighbours);
    } else {
        return new Canopy(origin);
    }
}

Canopy* CanopyClusteringAlg::canopy_walk(Point* origin, vector<Point*>& points, vector<Point*>& close_points, PRECISIONT max_canopy_dist, PRECISIONT max_close_dist, PRECISIONT min_step_dist, PRECISIONT max_num_canopy_walks, int& num_canopy_jumps){

    Canopy *c1;
    Canopy *c2;


    c1 = create_canopy(origin, points, close_points, max_canopy_dist, max_close_dist, true);

    //special case for which there is no walking, return the canopy immediatelly
    if(max_num_canopy_walks == 0){
        return c1;                                                               
    }

    c2 = create_canopy(c1->center, points, close_points, max_canopy_dist, max_close_dist, false);
    
    PRECISIONT dist = get_distance_between_points(c1->center, c2->center);

    {
        _log(logDEBUG2) << "Canopy walking, first step";
        _log(logDEBUG2) << *c1;
        _log(logDEBUG2) << *c2;
        _log(logDEBUG2) << "First potential jump correlation dist: " << dist;
    }

    int num_canopy_jumps_local = 0;
    while((dist > min_step_dist) && (num_canopy_jumps_local <= max_num_canopy_walks )){
        delete c1;
        c1=c2;

        num_canopy_jumps_local++; 

#pragma omp atomic
        num_canopy_jumps++;

        c2=create_canopy(c1->center, points, close_points, max_canopy_dist, max_close_dist, false);
        dist = get_distance_between_points(c1->center, c2->center); 
        {
            _log(logDEBUG2) << "Canopy walking, step: " << num_canopy_jumps_local;
            _log(logDEBUG2) << *c1;
            _log(logDEBUG2) << *c2;
            _log(logDEBUG2) << "distance: " << dist;
        }
    }

    //Now we know that c1 and c2 are close enough and we should choose the one that has more neighbours
    Canopy* final_canopy; 
    if(c1->neighbours.size() > c2->neighbours.size()){
        final_canopy = c1;
        delete c2;
    } else {
        final_canopy = c2;
        delete c1;
    }
    return final_canopy;
}

void CanopyClusteringAlg::filter_clusters_by_size(std::vector<Canopy*>& canopies_to_filter){

    vector<int> canopy_indexes_to_remove;

    for(int i=0; i < canopies_to_filter.size(); i++){
        Canopy* canopy = canopies_to_filter[i];
        if(canopy->neighbours.size() < 2){
            canopy_indexes_to_remove.push_back(i);
        }
    }

    std::sort(canopy_indexes_to_remove.begin(), canopy_indexes_to_remove.end());
    std::reverse(canopy_indexes_to_remove.begin(), canopy_indexes_to_remove.end());

    for(int i=0; i < canopy_indexes_to_remove.size(); i++)
        canopies_to_filter.erase(canopies_to_filter.begin() + canopy_indexes_to_remove[i]);

}

void CanopyClusteringAlg::cag_filter_max_top3_sample_contribution(PRECISIONT max_single_data_point_proportion, std::vector<Canopy*>& canopies_to_filter){

    vector<int> canopy_indexes_to_remove;

    for(int i=0; i < canopies_to_filter.size(); i++){
        Point* ccenter = canopies_to_filter[i]->center;
        //if(! ccenter->check_if_single_point_proportion_is_smaller_than(max_single_data_point_proportion) ){
        if(! ccenter->check_if_top_three_point_proportion_is_smaller_than(max_single_data_point_proportion) ){
            delete ccenter;
            canopy_indexes_to_remove.push_back(i);
        }
    }

    std::sort(canopy_indexes_to_remove.begin(), canopy_indexes_to_remove.end());
    std::reverse(canopy_indexes_to_remove.begin(), canopy_indexes_to_remove.end());

    for(int i=0; i < canopy_indexes_to_remove.size(); i++)
        canopies_to_filter.erase(canopies_to_filter.begin() + canopy_indexes_to_remove[i]);

}

void CanopyClusteringAlg::filter_clusters_by_zero_medians(int min_num_non_zero_medians, std::vector<Canopy*>& canopies_to_filter){

    vector<int> canopy_indexes_to_remove;

    for(int i=0; i < canopies_to_filter.size(); i++){
        Point* ccenter = canopies_to_filter[i]->center;
        if(! ccenter->check_if_num_non_zero_samples_is_greater_than_x(min_num_non_zero_medians) ){
            delete ccenter;
            canopy_indexes_to_remove.push_back(i);
        }
    }

    std::sort(canopy_indexes_to_remove.begin(), canopy_indexes_to_remove.end());
    std::reverse(canopy_indexes_to_remove.begin(), canopy_indexes_to_remove.end());

    for(int i=0; i < canopy_indexes_to_remove.size(); i++)
        canopies_to_filter.erase(canopies_to_filter.begin() + canopy_indexes_to_remove[i]);

}



void CanopyClusteringAlg::shuffle_points(vector<Point*>& points, vector<string>& priority_read_names){

    //Copy the read names from the vector into a map, we will be checking if the reads from the input file are in it and which position they should have
    std::map<string, int> priority_read_name__to_position;
    for(int i=0; i<priority_read_names.size(); i++)
        priority_read_name__to_position[priority_read_names[i]] = i;

    //Sort the points vector so that those reads that are in priority_read_names come first and in the order of the priority_read_names
    sort(points.begin(), points.end(), [&priority_read_name__to_position](const Point* p1, const Point* p2) -> bool {
        //integers are compared through "<" to get ascending sort
        auto p1_map_it = priority_read_name__to_position.find(p1->id);
        auto p2_map_it = priority_read_name__to_position.find(p2->id);

        //If both are in the priority read map - then compare integers directly
        if((p1_map_it != priority_read_name__to_position.end()) && (p2_map_it != priority_read_name__to_position.end()))
            return p1_map_it->second < p2_map_it->second;
        //If first is in priority read map (we know the second must not be there) then the first one goes before the other
        else if(p1_map_it != priority_read_name__to_position.end())
            return true;
        else
            return false; //That is the second read was in the priority read map and first wasn't, then second one goes before the first one
    });

    //Go through all points and find the first point that is NOT in the priority read names (using the map)
    //This is a somewhat clever approach: find_if returns "An iterator to the first element in the range for which pred does not return false."
    auto first_non_priority_pint_it = find_if(points.begin(), points.end(), [&priority_read_name__to_position](const Point* p) -> bool {
        return priority_read_name__to_position.find(p->id) == priority_read_name__to_position.end();
    });

    //Now shuffle the non prioritized pionts
    std::random_device rd;
    std::mt19937 random_generator(rd());
    std::shuffle(first_non_priority_pint_it, points.end(), random_generator);

}

std::vector<Canopy*> CanopyClusteringAlg::multi_core_run_clustering_on(vector<Point*>& points, vector<string>& priority_read_names, int num_threads, PRECISIONT max_canopy_dist, PRECISIONT max_close_dist, PRECISIONT max_merge_dist, PRECISIONT min_step_dist, int max_num_canopy_walks, PRECISIONT stop_after_num_seeds_processed, bool create_canopy_size_stats, string canopy_size_stats_fp, string not_processed_points_fp, bool show_progress_bar, TimeProfile& time_profile){

    _log(logINFO) << "";
    _log(logINFO) << "Algorithm Parameters:";
    _log(logINFO) << "max_canopy_dist:\t " << max_canopy_dist;
    _log(logINFO) << "max_close_dist:\t " << max_close_dist;
    _log(logINFO) << "max_merge_dist:\t " << max_merge_dist;
    _log(logINFO) << "min_step_dist:\t " << min_step_dist;
    _log(logINFO) << "max_num_canopy_walks:\t " << max_num_canopy_walks;
    _log(logINFO) << "";
    _log(logINFO) << "Early stopping:";
    _log(logINFO) << "stop_after_num_seeds_processed:\t " << stop_after_num_seeds_processed;
    _log(logINFO) << "";
    _log(logINFO) << "Priority reads";
    _log(logINFO) << "Number of reads with first priority:\t" << priority_read_names.size();

    _log(logPROGRESS) << "";
    _log(logPROGRESS) << "############ Shuffling ############";
    time_profile.start_timer("Shuffling");
    shuffle_points(points, priority_read_names);
    time_profile.stop_timer("Shuffling");

    _log(logPROGRESS) << "";
    _log(logPROGRESS) << "############ Creating Canopies ############";
    _log(logPROGRESS) << "To make the program stop and generate output send an interrupt signal by either:" ;
    _log(logPROGRESS) << "\t * running command \"kill -INT [ canopy PID ]\"";
    _log(logPROGRESS) << "\t * pressing \"CTRL + C\" in this terminal";
    _log(logPROGRESS) << "";

#pragma clang diagnostic push
#pragma ide diagnostic ignored "TemplateArgumentsIssues" //Clion is messing up, the set declaration is fine
    boost::unordered_set<Point*> marked_points;//Points that should not be investigated as origins
#pragma clang diagnostic pop
    vector<unsigned int> canopy_size_per_origin_num;//Contains size of the canopy created from origin by it's number, so first origin gave canopy of size 5, second origin gave canopy of size 8 and so on
    int last_progress_displayed_at_num_points = 0;

    std::vector<Canopy*> canopy_vector;

    vector<Point*> close_points;
    close_points.reserve(points.size());

    //
    //Create canopies
    //
    time_profile.start_timer("Clustering");
        
    ofstream canopy_size_stats_file;

    if(create_canopy_size_stats)
        canopy_size_stats_file.open(canopy_size_stats_fp.c_str(), ios::out | ios::trunc);

    int canopy_stats_row_num = 0;

    int num_canopy_jumps = 0;
    int num_collisions = 0;
    int num_seeds_processed = 0;

    int first_non_processed_origin_due_interruption = points.size();

    //Disable stop criterion if set to zero
    if(stop_after_num_seeds_processed == 0){
        stop_after_num_seeds_processed = points.size();
    }

#pragma omp parallel for shared(marked_points, canopy_vector, num_canopy_jumps, canopy_size_per_origin_num, num_collisions, num_seeds_processed, canopy_stats_row_num, terminate_called, first_non_processed_origin_due_interruption) firstprivate(close_points, max_canopy_dist, max_close_dist, max_merge_dist, min_step_dist, last_progress_displayed_at_num_points) schedule(dynamic,100)
    for(int origin_i = 0; origin_i < points.size(); origin_i++){

        //Early stopping after num of points
        if(num_seeds_processed > stop_after_num_seeds_processed){
            continue;
        }

        //Stop if exit signal received
        if(terminate_called){
            if(first_non_processed_origin_due_interruption > origin_i){
#pragma omp critical
                {
                    first_non_processed_origin_due_interruption = origin_i;
                }
            }
            continue;
        }

        Point* origin = points[origin_i]; 

        if(marked_points.find(origin) != marked_points.end())
            continue;

        //Show progress bar
        {
            //Only master thread executes this
            if(omp_get_thread_num() == 0){
                if(log_level >= logPROGRESS && show_progress_bar){
                    if(marked_points.size() > last_progress_displayed_at_num_points + stop_after_num_seeds_processed/100){
                        printProgBar(marked_points.size(),stop_after_num_seeds_processed * points.size());
                        last_progress_displayed_at_num_points = stop_after_num_seeds_processed;
                    }
                }
            }
        }




        {
            _log(logDEBUG) << "Unmarked points count: " << points.size() - marked_points.size() << " Marked points count: " << marked_points.size();
            _log(logDEBUG) << "points.size: " << points.size() << " origin_i: " << origin_i << " origin->id: " << origin->id ;

            _log(logDEBUG1) << "Current canopy origin: " << origin->id;
        }

        Canopy* final_canopy = canopy_walk(origin, points, close_points, max_canopy_dist, max_close_dist, min_step_dist, max_num_canopy_walks, num_canopy_jumps);

#pragma omp critical
        {
            //Do not commit anything if by chance another thread marked the current origin
            if(marked_points.find(origin) == marked_points.end()){

                //Add canopy
                marked_points.insert(origin);

                canopy_vector.push_back(final_canopy);

                for(Point* n : final_canopy->neighbours){
                    marked_points.insert(n);
                }

                //Statistics showing size of canopies per analyzed origin
                if(canopy_size_stats_file.is_open()){
                    canopy_size_stats_file << canopy_stats_row_num++  << "\t" << points.size() - marked_points.size() << "\t" << final_canopy->neighbours.size() << "\t" << num_collisions << endl;
                }

            } else {
                num_collisions++;
                delete final_canopy;
            }

            num_seeds_processed += 1;

        }

    }
    if(canopy_size_stats_file.is_open())
        canopy_size_stats_file.close();

    time_profile.stop_timer("Clustering");
    
    if(terminate_called && (not_processed_points_fp != "")){
        time_profile.start_timer("Saving unprocessed points file");
        _log(logERR) << "Received signal, clustering was stopped early, saving non processed points in file: " << not_processed_points_fp; 
        cout << "first_non_processed_origin_due_interruption:" << first_non_processed_origin_due_interruption << endl; 

        ofstream not_processed_points_file;
        not_processed_points_file.open(not_processed_points_fp.c_str(), ios::out | ios::trunc);
        for(int i = first_non_processed_origin_due_interruption; i < points.size(); i++){
            Point* point = points[i];
            if(marked_points.find(point) == marked_points.end()){
                not_processed_points_file << point->id;
                for(int j = 0; j < point->num_data_samples; j++){
                    not_processed_points_file << "\t" << point->sample_data[j];
                }
                not_processed_points_file << "\n";
            }
        }
        not_processed_points_file.close();
        _log(logERR) << "Unprocessed points saved"; 
        time_profile.stop_timer("Saving unprocessed points file");
    }


    _log(logINFO) << "";
    _log(logINFO) << "Avg. number of canopy walks: " << num_canopy_jumps/((PRECISIONT)canopy_vector.size());
    _log(logINFO) << "Number of all canopies before merging: " << canopy_vector.size();

    _log(logPROGRESS) << "";
    _log(logPROGRESS) << "############Removing canopies of size 1 to speed-up merging#############";
    int num_single_sample_canopies = 0;
    for(int i=0; i < canopy_vector.size(); ){
        if(canopy_vector[i]->neighbours.size() == 1){
            delete canopy_vector[i];
            canopy_vector.erase(canopy_vector.begin() + i);
            num_single_sample_canopies++;
        } else {
            i++;
        }
    }

    _log(logINFO) << "";
    _log(logINFO) << "Number of canopies which are removed due to having only 1 sample: " << num_single_sample_canopies;
    _log(logINFO) << "Number of canopies left after removal of single sample canpies: " << canopy_vector.size();

    int original_number_of_canopies = canopy_vector.size();

    //
    // Merge Canopies
    //
    std::vector<Canopy*> merged_canopy_vector;

    time_profile.start_timer("Merging");
    _log(logPROGRESS) << "";
    _log(logPROGRESS) << "############Merging canopies#############";


    //Actual merge 
    while(canopy_vector.size()){

        std::vector<Canopy*> canopies_to_merge;

        //This is the canopy we will look for partners for
        Canopy* c = *canopy_vector.rbegin();
        canopy_vector.pop_back();

        canopies_to_merge.push_back(c);

        //Get indexes of those canopies that are nearby
//#pragma omp parallel for shared(canopies_to_merge) 
        for(int i = 0; i < canopy_vector.size(); i++){

            Canopy* c2 = canopy_vector[i]; 

            PRECISIONT dist = get_distance_between_points(c->center, c2->center);

            if(dist < max_merge_dist){
//#pragma omp critical
                {
                    canopies_to_merge.push_back(c2);
                }
            }

        }

        //There are several canopies to merge, do it
        if( canopies_to_merge.size() > 1 ){

            vector<Point*> all_points_from_merged_canopies;
            
            for(Canopy* canopy : canopies_to_merge){
                for(Point* n : canopy->neighbours){
                    if(std::find(all_points_from_merged_canopies.begin(), all_points_from_merged_canopies.end(), n) == all_points_from_merged_canopies.end()){ //If the element hasn't been added already
                        all_points_from_merged_canopies.push_back(n);
                    }
                }
            }

            //Create new canopy
            Point* temp_merged_canopy_origin = get_centroid_of_points(all_points_from_merged_canopies);
            Canopy* merged_canopy = canopy_walk(temp_merged_canopy_origin, all_points_from_merged_canopies, close_points, max_canopy_dist, max_close_dist, min_step_dist, max_num_canopy_walks, num_canopy_jumps);
            delete temp_merged_canopy_origin;


            canopy_vector.push_back(merged_canopy);

            
            //Removed merged canopies //TODO might be slow
            for(Canopy* canopy : canopies_to_merge){
                canopy_vector.erase(remove(canopy_vector.begin(), canopy_vector.end(), canopy), canopy_vector.end());
                delete canopy;
            }


        //If no canopies were merged remove the canopy we compared against the others
        } else {

            merged_canopy_vector.push_back(c);

            //Show progress bar
            {
                if(log_level >= logPROGRESS && show_progress_bar){
                    if(original_number_of_canopies - canopy_vector.size() % 1000)
                        printProgBar(original_number_of_canopies - canopy_vector.size(), original_number_of_canopies );
                }
            }
        }
    }
    time_profile.stop_timer("Merging");

    _log(logINFO) << "";
    _log(logINFO) << "Number of canopies after merging: " << merged_canopy_vector.size();

    return merged_canopy_vector;

}

