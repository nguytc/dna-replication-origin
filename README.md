# DNA Replication Origin
This program determines the origin of DNA replication for prokaryotes and
automatically creates a new FASTA file with the sequence rearranged to begin
at the predicted origin.

Tim Nguyen (nguy1877@umn.edu)  
Minnesota Supercomputing Institute

## How to Use
    Usage: python origin.py -f <file.fasta> [options] 

    Options:
        -h, --help      
                show this help message and exit
        -f IN_FA_FILENAME, --fasta IN_FA_FILENAME
                Required input FASTA (.fasta/.fa) file
        -o OUT_FILENAME, --out OUT_FILENAME
                Basename of rearranged FASTA output file [adjusted]
        -s SIZE, --size SIZE  
                Starting window size. Decreases by 2 every iteration [524288]
        -d, --debug           
                Prints out all intermediate files obtained [False]
        -m, --multiscale      
                Use the multiscaling method (default; recommended)
        -w, --wavelet         
                Use the wavelet transform method instead

## Get
The main source code is contained in origin.py.  
Python 2.7+ is required for this program.

Modules:

- SeqIO from BioPython (Required)
    - parses the input FASTA file
    - http://biopython.org/wiki/SeqIO
- PyWavelets (Required only for wavelet transform method) 
    - performs wavelet transform decomposition
    - http://www.pybytes.com/pywavelets/
- robjects from rpy2 (Required only for wavelet transform method)
    - runs R code within Python
    - http://rpy.sourceforge.net/rpy2/doc-2.1/html/robjects.html

## Logic
Asymmetry in base composition between the leading and lagging strands in
bacterial chromosomes can be used to determine the origin of replication. 

A "multiscaling" method was implemented to take advantage of this asymmetry
to predict the origin. A second method using wavelet transforms was also
implemented as an alternative way to find the origin, but it is highly 
recommended that the multiscaling method be used.

##### Multiscaling
This method obtains the origin of DNA replication using a “multiscaling” 
approach. GC skews are considered and calculated for the indices of the 
sequence by focusing on windows of base pairs throughout the entirety of the 
sequence. The default starting window size is 2^19 = 524,288 base pairs.

- The GC skew is calculated for each window of base pairs throughout the 
sequence, and the skew for each given window is assigned to the index at the 
middle of the window. 
- The linear regression slopes of the GC skews are calculated for the same 
window size. 
- As given by the difference in GC skews between the leading and lagging 
strands in DNA replication, the largest slope is selected to be the origin 
of DNA replication, indicating the jump from the negative GC skews in the 
lagging strand to the positive GC skews in the leading strand. 
- After estimating the origin of replication for the current subsequence, 
the window of base pairs around that estimated origin index is used as the 
sequence for further consideration with the window size cut in half. 
- The method of GC skew and slope calculation is repeated for this narrowed
sequence, with the origin being estimated and the window of interest being 
narrowed in each iteration. 
- This continues until the window size decreases to 128 base pairs, which at 
this point, the origin obtained from this final iteration will be taken to be 
the final origin of DNA replication of the genome.

##### Wavelet Transform
This method obtains the origin of DNA replication using wavelet transforms.

- Uses multilevel decomposition Haar discrete wavelet transform method from 
the PyWavelets library to obtain the details coefficients for each level.
- R is used with rpy2 to obtain the quantitative data for the wavelet
transform regarding the segments in the levels and their detail coefficients.
- The R data is analyzed to determine the segments that are connected with
the segments in the higher levels to help visualize a tree structure, where the
level segments are the nodes.
- With the information stored in a data structure representing a tree, the
path within the tree with the largest sum of details coefficients will be
chosen as the path containing the origin. The largest sum path is chosen to
indicate the switch from lagging strand to leading strand.
- Depth-first traversal is used to keep track of the coefficent sums and to
determine the path with the largest sum.
- The origin of replication is chosen as the location at the bottom level of
the tree on the path with the highest sum.

## Notes
- It is highly recommended that the multiscaling approach is used over the
wavelet transform method, as it is more accurate and precise.
- The multiscaling approach will only require a few seconds to run, while
the wavelet transform method will require a few minutes. The wavelet transform
method also requires a significant amount of memory.
- Functions in the code that have names starting with an underscore indicate
that it is not used in the main execution. They are kept for testing purposes.
- The script does not support FASTA files with more than one sequence

Contact Tim with any questions, comments, or suggestions at nguy1877@umn.edu