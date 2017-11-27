cleanFasta:  Clean fasta file according to different criteria.  
  
Author: B. Nevado  
  
Usage: cleanFasta -infile in.fas -outfile out.fas -allowed acgt|iupac|missing [-thr 30] [ -prot 1 ]  
    -acgt: any site with ambiguous symbol (not acgt) is removed.  
    -iupac: any site with non-iupac symbol is removed.  
    -missing: any site with more than -thr % non-iupac symbols is removed.  
    -prot: denotes a protein alignment input (default is DNA). Only available with -allowed missing.  
 
Output: 1 filtered fasta file to -outfile.  
  
Requirements:  
    BIO++ (bpp-core and bpp-seq libraries)  
  
  
Installation (Linux):  
git clone https://github.com/brunonevado/cleanFasta  
cd cleanFasta  
make  
./cleanFasta  

