#include "CommandIndex.h"
#include "Index.h"

using namespace::std;

CommandIndex::CommandIndex()
{
    name = "index";
    description = "Create a reference index";
    argumentString = "fast(a|q)[.gz] ...";
    
    addOption("help", Option(Option::Boolean, "h", "Help", ""));
    addOption("kmer", Option(Option::Number, "k", "Kmer size. Hashes will be based on strings of this many nucleotides.", "11"));
    addOption("factor", Option(Option::Number, "c", "Compression factor. The number of min-hashes kept for each sequence will be its length divided by this number.", "100"));
    addOption("prefix", Option(Option::File, "p", "Output prefix (first input file used if unspecified). The suffix '.mash' will be appended.", ""));
}

int CommandIndex::run() const
{
    if ( arguments.size() == 0 || options.at("help").active )
    {
        print();
        return 0;
    }
    
    int kmer = options.at("kmer").getArgumentAsNumber();
    float factor = options.at("factor").getArgumentAsNumber();
    
    Index index;
    
    index.initFromSequence(arguments, kmer, factor);
    
    string prefix = options.at("prefix").argument.length() > 0 ? options.at("prefix").argument : arguments[0];
    
    index.writeToCapnp((prefix + ".mash").c_str());
    
    return 0;
}