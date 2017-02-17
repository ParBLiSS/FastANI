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

//Own includes
#include "map_parameters.hpp"
#include "map_stats.hpp"
#include "commonFunc.hpp"

//External includes
#include "argvparser.hpp"

namespace skch
{

  /**
   * @brief           Initialize the command line argument parser 
   * @param[out] cmd  command line parser object
   */
  void initCmdParser(CommandLineProcessing::ArgvParser &cmd)
  {
    cmd.setIntroductoryDescription("Approximate read mapper based on Jaccard similarity");

    cmd.setHelpOption("h", "help", "Print this help page");

    cmd.defineOption("subject", "an input reference file (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("subject","s");

    cmd.defineOption("subjectList", "a file containing list of reference files, one per line", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("subjectList","sl");

    cmd.defineOption("query", "an input query file (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("query","q");

    cmd.defineOption("kmer", "kmer size <= 16 [default 16 (DNA), 5 (AA)]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("kmer","k");

    cmd.defineOption("pval", "p-value cutoff, used to determine window/sketch sizes [default e-03]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("pval","p");

    cmd.defineOption("window", "window size [default : computed using pvalue cutoff]\n\
P-value is not considered if a window value is provided. Lower window size implies denser sketch", 
        ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("window","w");

    cmd.defineOption("fragLen", "fragment length [default : 3,000]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("fragLen","m");

    cmd.defineOption("perc_identity", "threshold for identity during mapping [default : 80]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("perc_identity","pi");

    cmd.defineOption("protein", "set alphabet type to proteins, default is nucleotides");
    cmd.defineOptionAlternative("protein","a");

    cmd.defineOption("output", "output file name", ArgvParser::OptionRequired | ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("output","o");
  }

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
      //Open file one by one
      for(auto &e : querySequences)
      {
        std::ifstream in(e);

        if (in.fail())
        {
          std::cerr << "ERROR, skch::validateInputFiles, Could not open " << e << "\n";
          exit(1);
        }
      }

      for(auto &e : refSequences)
      {
        std::ifstream in(e);

        if (in.fail())
        {
          std::cerr << "ERROR, skch::validateInputFiles, Could not open " << e << "\n";
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
    std::cerr << "Window size = " << parameters.windowSize << std::endl;
    std::cerr << "Fragment length = " << parameters.minReadLength << std::endl;
    std::cerr << "Alphabet = " << (parameters.alphabetSize == 4 ? "DNA" : "AA") << std::endl;
    std::cerr << "P-value = " << parameters.p_value << std::endl;
    std::cerr << "Percentage identity threshold = " << parameters.percentageIdentity << std::endl;
    std::cerr << "Mapping output file = " << parameters.outFileName << std::endl;
    std::cerr << ">>>>>>>>>>>>>>>>>>" << std::endl;
  }

  /**
   * @brief                   Parse the cmd line options
   * @param[in]   cmd
   * @param[out]  parameters  sketch parameters are saved here
   */
  void parseandSave(int argc, char** argv, 
      CommandLineProcessing::ArgvParser &cmd, 
      skch::Parameters &parameters)
  {
    int result = cmd.parse(argc, argv);

    //Make sure we get the right command line args
    if (result != ArgvParser::NoParserError)
    {
      std::cerr << cmd.parseErrorDescription(result) << "\n";
      exit(1);
    }
    else if (!cmd.foundOption("subject") && !cmd.foundOption("subjectList"))
    {
      std::cerr << "Provide reference file (s)\n";
      exit(1);
    }
    else if (!cmd.foundOption("query"))
    {
      std::cerr << "Provide reference file (s)\n";
      exit(1);
    }

    std::stringstream str;

    //Parse reference files
    if(cmd.foundOption("subject"))
    {
      std::string ref;

      str << cmd.optionValue("subject");
      str >> ref;

      parameters.refSequences.push_back(ref);
    }
    else //list of files
    {
      std::string listFile;

      str << cmd.optionValue("subjectList");
      str >> listFile;

      parseFileList(listFile, parameters.refSequences);
    }

    //Size of reference
    parameters.referenceSize = skch::CommonFunc::getReferenceSize(parameters.refSequences); 
    str.clear();

    //Parse query files
    {
      std::string query;

      str << cmd.optionValue("query");
      str >> query;

      parameters.querySequences.push_back(query);
    }
    
    str.clear();

    if(cmd.foundOption("protein"))
    {
      parameters.alphabetSize = 20;
    }
    else
      parameters.alphabetSize = 4;

    parameters.reportAll = true;

    //Parse algorithm parameters
    if(cmd.foundOption("kmer"))
    {
      str << cmd.optionValue("kmer");
      str >> parameters.kmerSize;
      str.clear();
    }
    else
    {
      if(parameters.alphabetSize == 4)
        parameters.kmerSize = 16;
      else
        parameters.kmerSize = 5;
    }

    if(cmd.foundOption("pval"))
    {
      str << cmd.optionValue("pval");
      str >> parameters.p_value;
      str.clear();
    }
    else
      parameters.p_value = 1e-03;

    if(cmd.foundOption("fragLen"))
    {
      str << cmd.optionValue("fragLen");
      str >> parameters.minReadLength;
      str.clear();
    }
    else
      parameters.minReadLength = 3000;

    if(cmd.foundOption("perc_identity"))
    {
      str << cmd.optionValue("perc_identity");
      str >> parameters.percentageIdentity;
      str.clear();
    }
    else
      parameters.percentageIdentity = 80;

    /*
     * Compute window size for sketching
     */

    if(cmd.foundOption("window"))
    {
      str << cmd.optionValue("window");
      str >> parameters.windowSize;
      str.clear();

      //Re-estimate p value
      int s = parameters.minReadLength * 2 / parameters.windowSize; 
      parameters.p_value = skch::Stat::estimate_pvalue (s, parameters.kmerSize, parameters.alphabetSize, 
          parameters.percentageIdentity, 
          parameters.minReadLength, parameters.referenceSize);
    }
    else
    {
      //Compute optimal window size
      parameters.windowSize = skch::Stat::recommendedWindowSize(parameters.p_value,
          parameters.kmerSize, parameters.alphabetSize,
          parameters.percentageIdentity,
          parameters.minReadLength, parameters.referenceSize);
    }

    str << cmd.optionValue("output");
    str >> parameters.outFileName;
    str.clear();

    printCmdOptions(parameters);

    //Check if files are valid
    validateInputFiles(parameters.querySequences, parameters.refSequences);
  }
}


#endif
