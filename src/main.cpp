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
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h> /* mmap() is defined in this header */
#include <fcntl.h>

#include <boost/program_options.hpp>
#include <boost/type_index.hpp>

#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/assign/std/vector.hpp>

#include <omp.h>

#include "precision_type.hpp"
#include "Point.hpp"
#include "CanopyClustering.hpp"
#include "Log.hpp"
#include "program_options_misc.hpp"
#include "signal_handlers.hpp"


using namespace std;
using namespace boost::program_options;
using namespace boost::assign;

ProfileMeasureType profile_measure;
DistanceMeasureType distance_measure;

int main(int argc, char* argv[])
{
    //
    //Initialization
    //
    
    //Set initial logging level
    log_level = logINFO;

    std::srand ( unsigned ( std::time(NULL) ) );

    //Preapre Time Profile
    TimeProfile time_profile;
    time_profile.start_timer("Total");

    //Prepare variables for command line input
    string input_file_path;
    string input_filter_file;
    string priority_reads_file_path; //This optional file contains names of the reads that will be used as canopy centers first, then clustering proceeds as before (points chosen by random)
    string output_clusters_file_path;
    string excut_out_file_path; //Optional file to which profiles from extended neighbourhoods of final profiles will be written
    string output_cluster_profiles_file;
    string output_cluster_prefix;
    string profile_measure_str;
    string distance_measure_str;
    int num_threads;
    double max_canopy_dist;
    const double max_close_dist = 0.6; //The value is hardcoded and the option to change it removed from CLI to not confuse users
    double max_merge_dist;
    const double min_step_dist = 0.001; //The value is hardcoded and the option to change it removed from CLI to not confuse users 
    string verbosity_option;
    int filter_min_obs;
    double filter_max_top3_sample_contribution;
    int cag_filter_min_sample_obs;
    double cag_filter_max_top3_sample_contribution;
    double stop_after_num_seeds_processed;
    double excut_threshold;
    bool dont_create_progress_stat_file;
    string progress_stat_file;
    string not_processed_profiles_file;
    bool show_progress_bar;
    bool print_time_statistics;
    bool die_on_kill;
    bool dont_use_mmap;
    int max_num_canopy_walks;
    int excut_min_size;
    int excut_max_size;
    vector<string> valid_profile_measure_values;
    valid_profile_measure_values += "median", "mean", "75Q", "80Q", "85Q", "90Q", "95Q";

    vector<string> valid_distance_measures;
    valid_distance_measures += "pearson", "spearman";

    //Define and read command line options
    options_description all_options_desc("Allowed options");
    options_description options_shown_in_help_desc("Allowed options");
    options_description options_not_shown_in_help_desc("Options not shown in help");
    options_description general_options_desc("General");
    options_description algorithm_param_options_desc("Algorithm Parameters");
    options_description filter_in_options_desc("Input filter parameters");
    options_description extended_cut_options_desc("Extended cut parameters");
    options_description filter_out_options_desc("Output filter parameters");
    options_description early_stop_options_desc("Early stopping");
    options_description misc_options_desc("Miscellaneous");


    general_options_desc.add_options()
        ("input_file_path,i", value<string>(&input_file_path), "Path to the input file")
        ("priority_reads_file_path", value<string>(&priority_reads_file_path)->default_value(""), "Path to (optional) file containing an ordered list (line by line) of read names according to which they should considered as cluster seeds or members. Reads from the input file not present in this file will be shuffled and considered second.")
        ("output_clusters_file_path,o", value<string>(&output_clusters_file_path)->default_value("clusters_out"), "Path to file to which clusters will be written")
        ("output_cluster_profiles_file,c", value<string>(&output_cluster_profiles_file)->default_value("profiles_out"), "Path to file to which cluster profiles will be written")
        ("cluster_name_prefix,p", value<string>(&output_cluster_prefix)->default_value("CAG"), "Prefix prepended to output cluster names")
        ("num_threads,n", value<int>(&num_threads)->default_value(4), "Number of cpu threads to use.")
        ("verbosity,v", value<string>(&verbosity_option)->default_value("info"), "Control how much information should be printed to the screen. Available levels according to their verbosity: error, progress, warn, info, debug, debug1.");

    algorithm_param_options_desc.add_options()
        ("distance_measure", value<string>(&distance_measure_str)->default_value("pearson"), "Distance measure used to cluster samples by, valid options are \"pearson\" and \"spearman\" for the respective correlation measures. \"pearson\" is used by default.")
        ("max_canopy_dist", value<double>(&max_canopy_dist)->default_value(0.1), "Max distance (f.e. pearson correlation difference) between a canopy center and a point included to the canopy")
        //This option is removed from CLI to avoid user confusion. The default value is hardcoded above
        //("max_close_dist", value<double>(&max_close_dist)->default_value(0.6), "Max pearson correlation difference between a canopy center and a point in which the point will be considered close to the canopy. As a heuristc, only points within this distance will be considered as potential neighbours during the canopy walk.")
        ("max_merge_dist", value<double>(&max_merge_dist)->default_value(0.1), "Max distance (f.e. pearson correlation difference) between two canopy centers in which the canopies should be merged. Please note, that the final canopy profiles are calculated after the merge step and consequently some final canopies might have profiles that are closer then max_merge_dist specifies.")
        //This option is removed from CLI to avoid user confusion. The default value is hardcoded above
        //("min_step_dist", value<double>(&min_step_dist)->default_value(0.001), "Min pearson correlation difference between canopy center and canopy centroid in which the centroid will be used as an origin for a new canpy (canopy walk). This is a stop criterion for canopy walk.")
        ("profile_measure", value<string>(&profile_measure_str)->default_value("75Q"), "Speicfies gene abundance measure should the algorithm use. Valid options are: \"median\", \"mean\", \"75Q\", \"80Q\", \"85Q\", \"90Q\", \"95Q\" where \"XXQ\" stands for XXth quantile measure.");

    filter_in_options_desc.add_options()
        ("filter_min_obs", value<int>(&filter_min_obs)->default_value(3), "Discard those profiles which have fewer than N non-zero samples. Setting it to 0 will disable the filter.")
        ("filter_max_top3_sample_contribution", value<double>(&filter_max_top3_sample_contribution)->default_value(0.9), "Discard those profiles for which top 3 samples constitute more than X fraction of the total signal. Setting it to 1 will disable the filter")
        ("input_filter_file", value<string>(&input_filter_file)->default_value(""), "The file to which profiles filtered out by either of the input filters will be written");

    filter_out_options_desc.add_options()
        ("cag_filter_min_sample_obs", value<int>(&cag_filter_min_sample_obs)->default_value(3), "Return only those canopies that have at least N non-zero cluster profile observations. Setting it to 0 will disable the filter.")
        ("cag_filter_max_top3_sample_contribution", value<double>(&cag_filter_max_top3_sample_contribution)->default_value(0.9), "Don't return canopies where top three(or less) samples constitute more than X fraction of the total profile signal. Setting it to 1 disables the filter.");

    extended_cut_options_desc.add_options()
        ("excut_out_file_path", value<string>(&excut_out_file_path)->default_value(""), "Path to output file with additional/extended neighbourhood of final/output profiles")
        ("excut_min_size", value<int>(&excut_min_size)->default_value(0), "Minimum size of a final/output profile for which extended cut neighbourhood will be returned")
        ("excut_max_size", value<int>(&excut_max_size)->default_value(1000000), "Maximum size of a final/output profile for which extended cut neighbourhood will be returned")
        ("excut_threshold", value<double>(&excut_threshold)->default_value(0.4), "Distance threshold below which samples will be considered to be in an extended neighbourhood of profile");

    early_stop_options_desc.add_options()
        ("stop_criteria", value<double>(&stop_after_num_seeds_processed)->default_value(50000), "Stop clustering after X number of seeds have been processed. Setting it to 0 will disable this stop criterion.");

    misc_options_desc.add_options()
        ("die_on_kill", bool_switch(&die_on_kill), "If set, after receiving a KILL signal, the program will die and no results will be produced. By default clustering will stop but clusters will be merged and partial results will be printed as usual.")
        ("dont_use_mmap", bool_switch(&dont_use_mmap), "If set, the program will not attempt to read in the entire file into memory but read it line by line. It will be slower but will potentially save lot of RAM.")
        ("not_processed_profiles_file", value<string>(&not_processed_profiles_file)->default_value(""), "Path to file to which unprocessed profiles will be dumped at KILL signal")
        ("print_time_statistics,t", bool_switch(&print_time_statistics), "Print wall clock time profiles of various analysis parts. This is not aggressive and won't increase compuatation time.")
        ("show_progress_bar,b", bool_switch(&show_progress_bar), "Show progress bar, nice if output is printed to console, don't use if you are redirecting to a file. Verbosity must be set to at least PROGRESS for it to have an effect.") 
        ("dont_create_progress_stat_file", bool_switch(&dont_create_progress_stat_file), "If set, the canopy progress file will not be created.")
        ("progress_stat_file", value<string>(&progress_stat_file)->default_value("canopy_progress.out"), "Name of the canopy size statistics file. To this file current progress after each processed seed profile will be dumped in format <index> <num_profiles_left> <this_canopy_size> <total_num_thread_collisions>")
        ("help", "write help message");

    options_not_shown_in_help_desc.add_options()
        //This option used to be removed from CLI to avoid user confusion and hardcoded above
        ("max_num_canopy_walks", value<int>(&max_num_canopy_walks)->default_value(6), "Max number of times the canopy will walk. This is a stop criterion for canopy walk.");

    all_options_desc.add(general_options_desc).add(algorithm_param_options_desc).add(filter_in_options_desc).add(filter_out_options_desc).add(early_stop_options_desc).add(misc_options_desc).add(options_not_shown_in_help_desc).add(extended_cut_options_desc);
    options_shown_in_help_desc.add(general_options_desc).add(algorithm_param_options_desc).add(filter_in_options_desc).add(filter_out_options_desc).add(early_stop_options_desc).add(misc_options_desc).add(extended_cut_options_desc);

    positional_options_description command_line_positional_desc;
    command_line_positional_desc.add("input_file_path",1);
    command_line_positional_desc.add("output_clusters_file_path",1);
    command_line_positional_desc.add("output_cluster_profiles_file",1);

    variables_map command_line_variable_map;
    store(command_line_parser(argc,argv).options(all_options_desc).positional(command_line_positional_desc).run(), command_line_variable_map);
    notify(command_line_variable_map);

    //
    //Verify command line input parameters
    //
    //verify_input_correctness(all_options_desc, command_line_variable_map);
    if (command_line_variable_map.count("help") || argc < 3) {
        cout << "Usage: canopy [options] PROFILES_INPUT_FILE CLUSTERS_OUTPUT_FILE" << endl << endl;;
        cout << options_shown_in_help_desc << "\n";
        exit(0);
    }

    check_if_file_is_readable("input_file_path", input_file_path);
    if(priority_reads_file_path != "")
        check_if_file_is_readable("priority_reads_file_path", priority_reads_file_path);
    check_if_file_is_writable("output_clusters_file_path", output_clusters_file_path);
    check_if_file_is_writable("output_cluster_profiles_file", output_cluster_profiles_file);
    check_if_file_is_writable("input_filter_file", input_filter_file);
    if(excut_out_file_path != "")
        check_if_file_is_writable("excut_out_file_path", excut_out_file_path);
    vector<string> valid_verbosities;
    valid_verbosities += "error", "progress", "warn", "info", "debug", "debug1", "debug2", "debug3";
    check_if_one_of("verbosity_option",verbosity_option, valid_verbosities);
    check_if_within_bounds("num_threads",num_threads,1,999);//Not exactly future proof, but let's put foolproofness first
    check_if_within_bounds("max_canopy_dist",max_canopy_dist,0.0,1.0);
    check_if_within_bounds("max_close_dist",max_close_dist,0.0,1.0);
    check_if_within_bounds("max_merge_dist",max_merge_dist,0.0,1.0);
    check_if_within_bounds("min_step_dist",min_step_dist,0.0,1.0);
    check_if_within_bounds("max_num_canopy_walks",max_num_canopy_walks,0,100);
    check_if_one_of("distance_measure", distance_measure_str, valid_distance_measures);
    check_if_one_of("profile_measure", profile_measure_str, valid_profile_measure_values);

    check_if_within_bounds("filter_min_obs",filter_min_obs,0,10000);
    check_if_within_bounds("filter_max_top3_sample_contribution",filter_max_top3_sample_contribution,0.0,1.0);
    check_if_within_bounds("cag_filter_min_sample_obs",cag_filter_min_sample_obs,0,10000);
    check_if_within_bounds("cag_filter_max_top3_sample_contribution",cag_filter_max_top3_sample_contribution,0.0,1.0);
    bool create_progress_stat_file = ! dont_create_progress_stat_file;
    if(create_progress_stat_file)
        check_if_file_is_writable("progress_stat_file",progress_stat_file);
    if(not_processed_profiles_file!= "")
        check_if_file_is_writable("not_processed_profiles_file",not_processed_profiles_file);
    
    //
    //Set appropriate distance measure method
    //
    if(distance_measure_str == "pearson"){
        distance_measure = PEARSON;
    } else if(distance_measure_str == "spearman"){
        distance_measure = SPEARMAN;
    } else {
        cout << "Unknown distance measure: \"" << distance_measure_str << "\"" << endl;
        cout << "This is most likely a programming error, please report this bug" << endl;
        exit(1);
    }


    //
    //Set appropriate profile measure method to the global var (ugh..)
    //
    if(profile_measure_str == "median"){
        profile_measure = MEDIAN;
    } else if(profile_measure_str == "mean"){
        profile_measure = MEAN;
    } else if(profile_measure_str == "75Q"){
        profile_measure = PERCENTILE_75;
    } else if(profile_measure_str == "80Q"){
        profile_measure = PERCENTILE_80;
    } else if(profile_measure_str == "85Q"){
        profile_measure = PERCENTILE_85;
    } else if(profile_measure_str == "90Q"){
        profile_measure = PERCENTILE_90;
    } else if(profile_measure_str == "95Q"){
        profile_measure = PERCENTILE_95;
    } else {
        cout << "Unknown type of profile measure method: \"" << profile_measure_str << "\"" << endl;
        cout << "This is most likely a programming error, please report this bug" << endl;
        exit(1);
    }

    //
    //Set user chosen logging level
    //
    if(verbosity_option == "error"){
        log_level = logERR;
    }else if(verbosity_option == "progress"){
        log_level = logPROGRESS;
    }else if(verbosity_option == "warn"){
        log_level = logWARN;
    }else if(verbosity_option == "info"){
        log_level = logINFO;
    }else if(verbosity_option == "debug"){
        log_level = logDEBUG;
    }else if(verbosity_option == "debug1"){
        log_level = logDEBUG1;
    }else if(verbosity_option == "debug2"){
        log_level = logDEBUG2;
    }else if(verbosity_option == "debug3"){
        log_level = logDEBUG3;
    }

    _log(logINFO) << "";
    _log(logINFO) << "Files:";
    _log(logINFO) << "input_file_path:\t " << input_file_path;
    _log(logINFO) << "priority_reads_file_path:\t" << priority_reads_file_path;
    _log(logINFO) << "output_cluster_profiles_file:\t " << output_cluster_profiles_file;
    _log(logINFO) << "progress_stat_file:\t " << progress_stat_file;
    _log(logINFO) << "not_processed_profiles_file:\t " << not_processed_profiles_file;
    _log(logINFO) << "input_filter_file:\t " << input_filter_file;
    _log(logINFO) << "";
    

    //Set signal handler
    if(die_on_kill) 
        signal(SIGINT, signal_callback_die_handler);
    else    
        signal(SIGINT, signal_callback_gentle_handler);

    //Set number of threads
    _log(logINFO) << "";
    _log(logINFO) << "General:";
    _log(logINFO) << "num_threads:\t " << num_threads;
    _log(logINFO) << "precision_type:\t " << boost::typeindex::type_id<PRECISIONT>().pretty_name();
    _log(logINFO) << "";

    omp_set_num_threads(num_threads);

    //
    //Parse priority point name file
    //
    time_profile.start_timer("Loading priority reads");
    vector<string> priority_read_names;
    if(priority_reads_file_path != ""){
        ifstream priority_reads_file (priority_reads_file_path);
        if(priority_reads_file.is_open())
        {
            string line;
            while(getline(priority_reads_file, line))
            {
                if(line.length()) {
                    priority_read_names.push_back(line);
                }
            }
            priority_reads_file.close();
        }
        _log(logINFO) << "";
        _log(logINFO) << "numbers of read priority reads:\t " << priority_read_names.size();
        _log(logINFO) << "";
        time_profile.stop_timer("Loading priority reads");
    }


    //
    //Parse point description file
    //
    
    //check if it is gz file
    bool is_gzipped_file = boost::algorithm::ends_with(input_file_path, "gz"); //In C++20 there is "ends_with" for string
    

    vector<Point*> points;
    vector<Point*> filtered_points;

    //TODO: this part could be improved by just fully commiting to boost iostreams instead of mixing and repeating code
    if((!dont_use_mmap) && (!is_gzipped_file)){//Read input file from mmap
        int point_file;
        char* point_file_mmap;
        struct stat statbuf;

        /* open the input file */
        point_file = open(input_file_path.c_str(), O_RDONLY);

        /* find size of input file */
        fstat(point_file,&statbuf);

        _log(logINFO) << "Reading file through mmap";
        time_profile.start_timer("File loading");

        point_file_mmap = (char*)mmap (0, statbuf.st_size, PROT_WRITE, MAP_PRIVATE, point_file, 0);

        time_profile.stop_timer("File loading");
        time_profile.start_timer("Reading profiles");

        _log(logINFO) << "File loaded into memory, generating profiles";

        char* line_start_ptr = point_file_mmap;
        char* line_end_ptr = point_file_mmap;
        char* mmap_end_ptr = point_file_mmap + statbuf.st_size;
        while(line_start_ptr < mmap_end_ptr){
            line_end_ptr = line_start_ptr;
            while(*line_end_ptr != '\n' && *line_end_ptr != '\r' && line_end_ptr < mmap_end_ptr){
                line_end_ptr++;
            }
            if(line_end_ptr != line_start_ptr){//Check if the line is not empty
                *line_end_ptr = '\0';
                //cout << line_start_ptr << endl;
                points.push_back(new Point(line_start_ptr));
            }
            line_start_ptr = ++line_end_ptr;
            die_if_true(terminate_called);
        }

        time_profile.stop_timer("Reading profiles");

        _log(logINFO) << "Profiles read, dropping file from memory";
        _log(logINFO) << "";

        /* drop the file from memory*/
        munmap(point_file_mmap, statbuf.st_size);
    } else if(is_gzipped_file){
        _log(logINFO) << "Reading file line by line from gzip";
        time_profile.start_timer("Loading file and reading profiles");

        std::ifstream point_file(input_file_path, std::ios_base::in | std::ios_base::binary);

        boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
        inbuf.push(boost::iostreams::gzip_decompressor());
        inbuf.push(point_file);
        //Convert streambuf to istream
        std::istream instream(&inbuf);

        std::string line;

        while (std::getline(instream, line)){
            if(line.length() < 2)
                break;

            points.push_back(new Point(line.c_str()));
            die_if_true(terminate_called);
        }

        time_profile.stop_timer("Loading file and reading profiles");
        _log(logINFO) << "";
    } else {

        _log(logINFO) << "Reading file line by line";
        time_profile.start_timer("Loading file and reading profiles");

        std::ifstream point_file(input_file_path);
        std::string line;

        while (std::getline(point_file, line)){
            if(line.length() < 2)
                break;

            points.push_back(new Point(line.c_str()));
            die_if_true(terminate_called);
        }

        time_profile.stop_timer("Loading file and reading profiles");
        _log(logINFO) << "";
    }
    
    _log(logINFO) << "Running basic validation of profiles";
    _log(logINFO) << "filter_max_top3_sample_contribution:\t " << filter_max_top3_sample_contribution;
    _log(logINFO) << "filter_min_obs:\t " << filter_min_obs;
    _log(logINFO) << "";

    time_profile.start_timer("Profiles validation");
    verify_proper_point_input_or_die(points);
    time_profile.stop_timer("Profiles validation");

    vector<Point*> points_filtered_out_due_to_three_point_proportion_filter;
    vector<Point*> points_filtered_out_due_to_num_non_zero_samples_filter;
    set<Point*> filtered_out_points;

    time_profile.start_timer("Input profiles filtering");

#pragma omp parallel for shared(points_filtered_out_due_to_three_point_proportion_filter, points_filtered_out_due_to_num_non_zero_samples_filter, filtered_points, filtered_out_points)
    for(int i = 0; i < points.size(); i++){
        //Both filters are set
        if((filter_min_obs > 0) && (filter_max_top3_sample_contribution< 0.9999)){
            bool point_is_valid = true;
                
            if( ! points[i]->check_if_num_non_zero_samples_is_greater_than_x(filter_min_obs))
            {
#pragma omp critical
                {
                    points_filtered_out_due_to_num_non_zero_samples_filter.push_back(points[i]);
                    filtered_out_points.insert(points[i]);
                    point_is_valid = false;
                }
            }

            if( ! points[i]->check_if_top_three_point_proportion_is_smaller_than(filter_max_top3_sample_contribution))
            {
#pragma omp critical
                {
                    points_filtered_out_due_to_three_point_proportion_filter.push_back(points[i]);
                    filtered_out_points.insert(points[i]);
                    point_is_valid = false;
                }
            }

            if(point_is_valid){
#pragma omp critical
                filtered_points.push_back(points[i]);
            }
        } else if (filter_min_obs > 0){ 
            if(points[i]->check_if_num_non_zero_samples_is_greater_than_x(filter_min_obs)){
#pragma omp critical
                filtered_points.push_back(points[i]);
            }
            else 
            {
#pragma omp critical
                {
                    points_filtered_out_due_to_num_non_zero_samples_filter.push_back(points[i]);
                    filtered_out_points.insert(points[i]);
                }
            }
        } else if (filter_max_top3_sample_contribution< 0.9999){ 
            if(points[i]->check_if_top_three_point_proportion_is_smaller_than(filter_max_top3_sample_contribution)){ 
#pragma omp critical
                filtered_points.push_back(points[i]);
            }
            else 
            {
#pragma omp critical
                {
                    points_filtered_out_due_to_three_point_proportion_filter.push_back(points[i]);
                    filtered_out_points.insert(points[i]);
                }
            }
        }
    }

    if(input_filter_file != ""){
        ofstream filtered_point_file;
        filtered_point_file.open(input_filter_file.c_str(), ios::out | ios::trunc);
        filtered_point_file << "#filtered_profile_id\tinput_filter_name\n";
        for(int i = 0; i < points_filtered_out_due_to_num_non_zero_samples_filter.size(); i++){
            filtered_point_file << points[i]->id << "\t" << "min_observations_filter" << "\n";
        }
        for(int i = 0; i < points_filtered_out_due_to_three_point_proportion_filter.size(); i++){
            filtered_point_file << points[i]->id << "\t" << "max_top3_sample_contribution_filter" << "\n";
        }
        filtered_point_file.close();
    }

    time_profile.stop_timer("Input profiles filtering");
    _log(logINFO) << "Number of profiles filtered out due to three sample signal contribution filter: " << points_filtered_out_due_to_three_point_proportion_filter.size();
    _log(logINFO) << "Number of profiles filtered out due to non zero samples number filter: " << points_filtered_out_due_to_num_non_zero_samples_filter.size();
    _log(logINFO) << "Number of profiles filtered out: " << points.size() - filtered_points.size(); 

    _log(logINFO) << "Finished input profiles processing";
    
    _log(logINFO) << "Number of profiles after filtering: " << filtered_points.size();

    //Sometimes these filters will remove over 50% of the dataset, let's free filtered out points now
    _log(logINFO) << "Relseasing filtered out points";
    for(Point* point: filtered_out_points)
        delete point;
    _log(logINFO) << "Relseasing filtered out points: Done";

    //Do not use "points" or points_filtered_out_due_to_three_point_proportion_filter or points_filtered_out_due_to_num_non_zero_samples_filter at this point
    
    
    die_if_true(terminate_called);
    die_if_true(filtered_points.size() < 1);

    //
    //This will precompute values for quicker pearson correlation calculation
    //This is a bit clumsy but helps prevent huge memory spikes
    //
    _log(logINFO) << "Precomputing pearson correlation data to speed up distance calculations";
    time_profile.start_timer("Precomputing pearson correlation data");
    for(int i = 0; i < filtered_points.size(); i++){
        filtered_points[i]->allocate_and_precompute_corr_data(distance_measure);
    }
    time_profile.stop_timer("Precomputing pearson correlation data");
    
    die_if_true(terminate_called);


    //
    //Run Canopy Clustering
    //
    std::vector<Canopy*> canopies;

    canopies = CanopyClusteringAlg::multi_core_run_clustering_on(filtered_points, priority_read_names, num_threads, max_canopy_dist, max_close_dist, max_merge_dist, min_step_dist, max_num_canopy_walks, stop_after_num_seeds_processed, create_progress_stat_file, progress_stat_file, not_processed_profiles_file, show_progress_bar, time_profile);

    _log(logINFO) << "Finished clustering";

    //
    //Filter out canopies
    //

    if(cag_filter_min_sample_obs){
        time_profile.start_timer("Filtering canopies by minimum number of sample detections");
        CanopyClusteringAlg::filter_clusters_by_zero_medians(cag_filter_min_sample_obs, canopies);
        _log(logINFO) << "Finished filtering for minimum number of sample detections, number of canopies:" << canopies.size();
        time_profile.stop_timer("Filtering canopies by minimum number of sample detections");
    }


    if(cag_filter_max_top3_sample_contribution< 0.99999){ //It's due to a double comparison
        time_profile.start_timer("Filtering canopies by three sample signal contribution proportion");
        CanopyClusteringAlg::cag_filter_max_top3_sample_contribution(cag_filter_max_top3_sample_contribution, canopies);
        _log(logINFO) << "Finished filtering by three sample signal contribution proportion, number of canopies:" << canopies.size();
        time_profile.stop_timer("Filtering canopies by three sample signal contribution proportion");
    }

    {
        time_profile.start_timer("Filtering canopies by size");
        CanopyClusteringAlg::filter_clusters_by_size(canopies);
        _log(logINFO) << "Finished filtering by size(number of neighbours must be bigger than 1), number of canopies:" << canopies.size();
        time_profile.stop_timer("Filtering canopies by size");
    }

    _log(logPROGRESS) << "";
    _log(logPROGRESS) << "####################Writing Results####################" ;
    ofstream out_file;

    sort(canopies.begin(), canopies.end(), compare_canopy_ptrs_by_canopy_size);

    int num_digits = ceil(log10(canopies.size()));
    cout << std::setfill('0');

    //
    //Number canopies that survived filtering so that the naming is consistent
    //
    vector<std::tuple<string, Canopy*> > named_canopies; 

    int i = 1;
    for(Canopy* c : canopies){
        stringstream canopy_name_buffer; 
        canopy_name_buffer << output_cluster_prefix << std::setw(num_digits) << std::setfill('0') << i;
        named_canopies.push_back(make_tuple(canopy_name_buffer.str(), c));
        i++;
    }

    out_file.open(output_clusters_file_path.c_str(), ios::out | ios::trunc);
    for(auto [canopy_name, canopy] : named_canopies){

        for(Point* p : canopy->neighbours){
            out_file << canopy_name << "\t" << p->id << "\n";
        }
    }
    out_file.close();

    out_file.open(output_cluster_profiles_file.c_str(), ios::out | ios::trunc);
    for(auto [canopy_name, canopy] : named_canopies){

        out_file << canopy_name << "\t";
        
        for(int j=0; j < canopy->center->num_data_samples; j++){
            out_file << canopy->center->sample_data[j] << "\t" ;
        }

        i++;
        out_file << "\n";
    }
    out_file.close();

    //
    // Extended cut step
    //
    if(excut_out_file_path != ""){

        _log(logPROGRESS) << "Producing extended cut profiles";
        _log(logDEBUG) << "excut_min_size: " << excut_min_size << " excut_max_size: " << excut_max_size;
        _log(logDEBUG) << "excut_threshold: " << excut_threshold;
        out_file.open(excut_out_file_path.c_str(), ios::out | ios::trunc);

        for(auto [canopy_name, canopy] : named_canopies){

            int num_samples = canopy->neighbours.size();
            if((num_samples < excut_max_size) && (num_samples > excut_min_size)){
                for(auto profile : filtered_points){
                    double distance = get_distance_between_points(canopy->center, profile);
                    if(distance < excut_threshold){

                        out_file << canopy_name << "\t" << distance << "\t" << profile->id << "\t";
                        
                        for(int j=0; j < canopy->center->num_data_samples; j++){
                            out_file << canopy->center->sample_data[j] << "\t" ;
                        }

                        i++;
                        out_file << "\n";

                    }
                }
            }
        }

        out_file.close();

    } else {
        _log(logPROGRESS) << "Not producing extended cut profiles";
    }


    //
    //Clean up
    //
    for(Canopy* c : canopies)
        delete c;


    time_profile.stop_timer("Total");
    //Write output statistics
    if(print_time_statistics){
        cout << time_profile << endl;
    }


    return 0;
}
