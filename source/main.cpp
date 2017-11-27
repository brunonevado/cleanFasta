
#include <iostream>

#include <Bpp/Seq/Alphabet.all> /* this includes all alphabets in one shot */
#include <Bpp/Seq/Container.all> /* this includes all containers */
#include <Bpp/Seq/AlphabetIndex/DefaultNucleotideScore.h>
#include <Bpp/Seq/Io.all> /* this includes all sequence readers and writers */
#include <Bpp/Seq/Site.h>
#include <Bpp/Seq/SiteTools.h>

#include "args.h"

// 15.11.2016: added protein option, only works with missing threshold


void help(){
    std::cout << "###################\n  cleanFasta v 15112016 \n###################" << std::endl;;
    std::cout << "Clean fasta file according to different criteria." << std::endl;;
    std::cout << "Usage: bpp_cleanFasta -infile in.fas -outfile out.fas -allowed acgt|iupac|missing [-thr 30] [ -prot 1 ]" << std::endl;
    std::cout << " -acgt: any site with ambiguous symbol (not acgt) is removed." << std::endl;
    std::cout << " -iupac: any site with non-iupac symbol (gaps) is removed." << std::endl;
    std::cout << " -missing: any site with more than -thr % non-iupac symbols (gaps) is removed." << std::endl;
    std::cout << " -prot: denotes a protein alignment input (default is DNA). Only available with -allowed missing." << std::endl;
    std::cout << "Requires Bio++ core & seq." << std::endl;

    
};


int main(int argc, const char * argv[])
{

    
    sargs myargs;
    try{
        myargs = args::getargs(argc, argv, std::vector<std::string> {"infile", "outfile", "allowed"}, std::vector<std::string> {}, std::vector<std::string>  {}, std::string  {"prot"}, std::string  {"thr"});
    }
    catch (std::string e){
        std::cerr << "ERROR: " << e << std::endl;
        help();
        exit(1);
    }
    // check that allowed is either acgt or iupac
    if( myargs.args_string.at(2) != "acgt" && myargs.args_string.at(2) != "iupac" && myargs.args_string.at(2) != "missing" ){
        std::cerr << "ERROR: option -allowed must be 'acgt' or 'iupac' or 'missing'" << std::endl;
        help();
        exit(1);
    }
    int aMaxMissing = 0;
    if ( myargs.args_string.at(2) == "missing" ){
        try{
            aMaxMissing = myargs.args_int_optional.at(0);
        }
        catch(...){
            std::cerr << "ERROR: for option 'missing', optional int argument 'thr' must be defined" << std::endl;
            exit(1);
        }
        std::clog << "Removing sites with more than  " << aMaxMissing << "% missing data (N/-)" << std::endl;
    }
    
    std::string aInfile =  myargs.args_string.at(0);
    std::string aOutfile = myargs.args_string.at(1);
    std::string aAllowed = myargs.args_string.at(2);
    
    bool isDNA = (myargs.args_string_optional.size() != 1) ? true : false;
    
    if (!isDNA && myargs.args_string.at(2) != "missing" ){
        std::cerr << "ERROR: for protein alignments only the 'missing' option is available." << std::endl;
        exit(1);
    }
    std::cout << "Assuming " << std::string ((isDNA) ? "DNA" : "PROTEIN") << " alignment, reading file " << aInfile << std::endl;
    
    bpp::Fasta fasReader(-1);
    bpp::OrderedSequenceContainer *sequences = (isDNA) ? fasReader.readSequences(aInfile, &bpp::AlphabetTools::DNA_ALPHABET) : fasReader.readSequences(aInfile, &bpp::AlphabetTools::PROTEIN_ALPHABET);
    bpp::SiteContainer *vsc2 = new bpp::VectorSiteContainer(*sequences);
    
    std::clog << "Read fasta file " << aInfile << " with " << vsc2->getNumberOfSequences() << " sequences "
    << vsc2->getNumberOfSites() << " bp long." << std::endl;

    if( aAllowed == "iupac" ){
        for ( int isite = int (vsc2->getNumberOfSites() - 1); isite >= 0; isite-- ){
            bpp::ApplicationTools::displayGauge(isite, int (vsc2->getNumberOfSites() - 1), '.', "Sites missing " );
            for ( int iline = 0; iline < vsc2->getNumberOfSequences() ; iline++ ) {
                if( vsc2->getSite(isite).getValue(iline) > 13
                   || vsc2->getSite(isite).getValue(iline) < 0 ){
                    vsc2->deleteSite(isite);
                    break;
                }
            }
        }
    }
    else if ( aAllowed == "acgt"  ){
        for ( int isite = int (vsc2->getNumberOfSites() - 1); isite >= 0; isite-- ){
            bpp::ApplicationTools::displayGauge(isite, int (vsc2->getNumberOfSites() - 1),  '.', "Sites missing " );
            for ( int iline = 0; iline < vsc2->getNumberOfSequences() ; iline++ ) {
                if( vsc2->getSite(isite).getValue(iline) > 3
                   || vsc2->getSite(isite).getValue(iline) < 0 ){
                    vsc2->deleteSite(isite);
                    break;
                }
            }
        }
    }
    if( aAllowed == "missing" ){
        for ( int isite = int (vsc2->getNumberOfSites() - 1); isite >= 0; isite-- ){
            bpp::ApplicationTools::displayGauge(isite, int (vsc2->getNumberOfSites() - 1), '.', "Sites missing " );
            int miss = 0;
            for ( int iline = 0; iline < vsc2->getNumberOfSequences() ; iline++ ) {
                if( vsc2->getSite(isite).getValue(iline) > 13
                   || vsc2->getSite(isite).getValue(iline) < 0 ){
                    miss++;
                    if ( 100 * float(miss) / float(vsc2->getNumberOfSequences()) > aMaxMissing   ){
                        break;
                    }
                }
            }
            if ( 100 * float(miss) / float(vsc2->getNumberOfSequences()) > aMaxMissing   ){
                vsc2->deleteSite(isite);
            }
        }
    }
    std::clog << "Fasta file now has " << vsc2->getNumberOfSequences() << " sequences "
    << vsc2->getNumberOfSites() << " bp long." << std::endl;
    fasReader.writeSequences(aOutfile, *vsc2);
    std::clog << "Fasta file written to " <<  aOutfile << std::endl;
    
}
