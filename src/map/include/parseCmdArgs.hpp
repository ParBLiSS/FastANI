/**
 * @file    parseCmdArgs.hpp
 * @brief   Functionality related to command line parsing for indexing and mapping
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef PARSE_CMD_HPP 
#define PARSE_CMD_HPP

#include <iostream>
#include <string>
#include <fstream>
#include <cassert>

//Own includes
#include "map/include/map_parameters.hpp"
#include "map/include/map_stats.hpp"
#include "map/include/commonFunc.hpp"

//External includes
#include "common/clipp.h"

namespace skch
{

  /**
   * @brief                   Parse the file which has list of reference or query files
   * @param[in]   fileToRead  File containing list of ref/query files 
   * @param[out]  fileList    List of files will be saved in this vector   
   */
  template <typename VEC>
    void parseFileList(std::string &fileToRead, VEC &fileList)
    {
      std::string line;

      std::ifstream in(fileToRead);

      if (in.fail())
      {
        std::cerr << "ERROR, skch::parseFileList, Could not open " << fileToRead << "\n";
        exit(1);
      }

      while (std::getline(in, line))
      {
        //trim whitespaces
        skch::CommonFunc::trim (line);

        if (line.length() > 0)        //avoid empty strings
          fileList.push_back(line);
      }
    }

  /**
   * @brief                     validate the reference and query file(s)
   * @param[in] querySequences  vector containing query file names
   * @param[in] refSequences    vector containing reference file names
   */
  template <typename VEC>
    void validateInputFiles(VEC &querySequences, VEC &refSequences)
    {
      if (querySequences.size() == 0 || refSequences.size() == 0)
      {
        std::cerr << "ERROR, skch::validateInputFiles, Count of query and ref genomes should be non-zero" << std::endl;
        exit(1);
      }

      //Open file one by one
      for(auto &e : querySequences)
      {
        std::ifstream in(e);

        if (in.fail())
        {
          std::cerr << "ERROR, skch::validateInputFiles, Could not open " << e << std::endl;
          exit(1);
        }
      }

      for(auto &e : refSequences)
      {
        std::ifstream in(e);

        if (in.fail())
        {
          std::cerr << "ERROR, skch::validateInputFiles, Could not open " << e << std::endl;
          exit(1);
        }
      }
    }

  /**
   * @brief                   Print the parsed cmd line options
   * @param[in]  parameters   parameters parsed from command line
   */
  void printCmdOptions(skch::Parameters &parameters)
  {
    std::cerr << ">>>>>>>>>>>>>>>>>>" << std::endl;
    std::cerr << "Reference = " << parameters.refSequences << std::endl;
    std::cerr << "Query = " << parameters.querySequences << std::endl;
    std::cerr << "Kmer size = " << parameters.kmerSize << std::endl;
    std::cerr << "Fragment length = " << parameters.minReadLength << std::endl;
    std::cerr << "Threads = " << parameters.threads << std::endl;
    std::cerr << "ANI output file = " << parameters.outFileName << std::endl;
    std::cerr << ">>>>>>>>>>>>>>>>>>" << std::endl;
  }

  /**
   * @brief                   Parse the cmd line options
   * @param[in]   cmd
   * @param[out]  parameters  sketch parameters are saved here
   */
  void parseandSave(int argc, char** argv, 
      skch::Parameters &parameters)
  {
    //defaults
    parameters.kmerSize = 16;
    parameters.minReadLength = 3000;
    parameters.alphabetSize = 4;
    parameters.minFraction = 0.2;
    parameters.threads = 1;
    parameters.p_value = 1e-03;
    parameters.percentageIdentity = 80;
    parameters.visualize = false;
    parameters.matrixOutput = false;
    parameters.referenceSize = 5000000;
    parameters.reportAll = true; //we need all mappings per fragment, not just best 1% as in mashmap


    std::string refName, refList;
    std::string qryName, qryList;
    bool versioncheck = false;
    bool help = false;

    auto help_cmd = clipp::option("-h", "--help").set(help).doc("print this help page");
    auto ref_cmd = (clipp::option("-r", "--ref") & clipp::value("value", refName)) % "reference genome (fasta/fastq)[.gz]";
    auto refList_cmd = (clipp::option("--rl", "--refList") & clipp::value("value", refList)) % "a file containing list of reference genome files, one genome per line";
    auto qry_cmd = (clipp::option("-q", "--query") & clipp::value("value", qryName)) % "query genome (fasta/fastq)[.gz]";
    auto qryList_cmd = (clipp::option("--ql", "--queryList") & clipp::value("value", qryList)) % "a file containing list of query genome files, one genome per line";
    auto kmer_cmd = (clipp::option("-k", "--kmer") & clipp::value("value", parameters.kmerSize)) % "kmer size <= 16 [default : 16]";
    auto thread_cmd = (clipp::option("-t", "--threads") & clipp::value("value", parameters.threads)) % "thread count for parallel execution [default : 1]";
    auto fraglen_cmd = (clipp::option("--fragLen") & clipp::value("value", parameters.minReadLength)) % "fragment length [default : 3,000]";
    auto minfraction_cmd = (clipp::option("--minFraction") & clipp::value("value", parameters.minFraction)) % "minimum fraction of genome that must be shared for trusting ANI. If reference and query genome size differ, smaller one among the two is considered. [default : 0.2]";
    auto visualize_cmd = clipp::option("--visualize").set(parameters.visualize).doc("output mappings for visualization, can be enabled for single genome to single genome comparison only [disabled by default]");
    auto matrix_cmd = clipp::option("--matrix").set(parameters.matrixOutput).doc("also output ANI values as lower triangular matrix (format inspired from phylip). If enabled, you should expect an output file with .matrix extension [disabled by default]");
    auto output_cmd = (clipp::option("-o", "--output") & clipp::value("value", parameters.outFileName)) % "output file name";
    auto version_cmd = clipp::option("-v", "--version").set(versioncheck).doc("show version");

    auto cli =
      (
       help_cmd,
       ref_cmd,
       refList_cmd,
       qry_cmd,
       qryList_cmd,
       kmer_cmd,
       thread_cmd,
       fraglen_cmd,
       minfraction_cmd,
       visualize_cmd,
       matrix_cmd,
       output_cmd,
       version_cmd
      );

    //with formatting options
    auto fmt = clipp::doc_formatting{}
    .first_column(0)
      .doc_column(5)
      .last_column(80);

    std::string description = "fastANI is a fast alignment-free implementation for computing whole-genome Average Nucleotide Identity (ANI) between genomes\n-----------------\nExample usage:\n$ fastANI -q genome1.fa -r genome2.fa -o output.txt\n$ fastANI -q genome1.fa --rl genome_list.txt -o output.txt";

    if(!clipp::parse(argc, argv, cli))
    {
      //print help page
      clipp::operator<<(std::cout, clipp::make_man_page(cli, argv[0], fmt).prepend_section("-----------------", description)) << std::endl;
      exit(1);
    }

    if (help)
    {
      clipp::operator<<(std::cout, clipp::make_man_page(cli, argv[0], fmt).prepend_section("-----------------", description)) << std::endl;
      exit(0);
    }

    if (versioncheck)
    {
      std::cerr << "version 1.33\n\n";
      exit(0);
    }

    if (refName == "" && refList == "")
    {
      std::cerr << "Provide reference file (s)\n";
      exit(1);
    }

    if (qryName == "" && qryList == "")
    {
      std::cerr << "Provide query file (s)\n";
      exit(1);
    }

    if (refName != "")
      parameters.refSequences.push_back(refName);
    else
      parseFileList(refList, parameters.refSequences);

    if (qryName != "")
      parameters.querySequences.push_back(qryName);
    else
      parseFileList(qryList, parameters.querySequences);

    assert(parameters.minFraction >= 0.0 && parameters.minFraction <= 1.0);

    //Compute optimal window size
    parameters.windowSize = skch::Stat::recommendedWindowSize(parameters.p_value,
        parameters.kmerSize, parameters.alphabetSize,
        parameters.percentageIdentity,
        parameters.minReadLength, parameters.referenceSize);

    printCmdOptions(parameters);

    //Check if files are valid
    validateInputFiles(parameters.querySequences, parameters.refSequences);
  }
}


#endif
