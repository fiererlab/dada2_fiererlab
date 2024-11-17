
# dada2 tutorial with MiSeq dataset for Fierer Lab

*This tutorial was created by Angela Oliverio and Hannah Holland-Moritz,
and is maintained by Cliff Bueno de Mesquita and other current members
of the Fierer Lab*  
*Updated November 16th, 2024*

This pipeline runs the dada2 workflow for Big Data (paired-end) from
RStudio on the microbe server.

We suggest opening the dada2 tutorial online to understand more about
each step. The original pipeline on which this tutorial is based can be
found here: <https://benjjneb.github.io/dada2/bigdata_paired.html>

| <span>                                                                                                                                                                                                                                                       |
|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **NOTE:** There is a slightly different pipeline for ITS and non-“Big data” workflows. The non-“Big data” pipeline, in particular, has very nice detailed explanations for each step and can be found here: <https://benjjneb.github.io/dada2/tutorial.html> |
| <span>                                                                                                                                                                                                                                                       |

## Preliminary Checklist (part 0) - Before You Begin

1.  Check to make sure you know what your target ‘AMPLICON’ length is.
    This can vary between primer sets, as well as WITHIN primer sets.
    For example, ITS (internal transcribed spacer) amplicons can vary
    from \~100 bps to 300 bps.

    For examples regarding commonly used primer sets (515f/806r, Fungal
    ITS2, 1391f/EukBr) see protocols on the Earth Microbiome Project
    website:
    <http://press.igsb.anl.gov/earthmicrobiome/protocols-and-standards/>

2.  Check to make sure you know how long your reads should be (i.e., how
    long should the reads be coming off the sequencer?) This is not the
    same as fragment length, as many times, especially with longer
    fragments, the entire fragment is not being sequenced in one
    direction. When long *amplicons* are not sequenced with a *read
    length* that allows for substantial overlap between the forward and
    reverse read, you can potentially insert biases into the data. If
    you intend to merge your paired end reads, ensure that your read
    length is appropriate. For example, with a MiSeq 2 x 150, 300 cycle
    kit, you will get bidirectional reads of 150 base pairs.

3.  Make note of which sequencing platform was used, as this can impact
    both read quality and downstream analysis.

4.  Decide which database is best suited for your analysis needs. Note
    that DADA2 requires databases be in a custom format! If a custom
    database is required, further formatting will be needed to ensure
    that it can run correctly in dada2.

    See the following link for details regarding database formatting:
    <https://benjjneb.github.io/dada2/training.html#formatting-custom-databases>

5.  For additional tutorials and reporting issues, please see link
    below:  
    dada2 tutorial: <https://benjjneb.github.io/dada2/tutorial.html>  
    dada2 pipeline issues\*:
    <https://github.com/fiererlab/dada2_fiererlab/issues>

    \*Note by default, only ‘OPEN’ issues are shown. You can look at all
    issues by removing “is:open” in the search bar at the top.

## Set up (part 1) - Steps before starting pipeline

If you are running it through Fierer Lab “microbe” server:

#### Logging in:

To use the microbe server, open a terminal window, and type and hit
return:

``` bash
ssh <your microbe user name>@microbe.colorado.edu 
```

#### For your first time:

*Setting up your password:*

When you log in for the first time, you will have to set a new password.
First, log in using your temporary password. The command prompt will
then ask you to write a new password. Type it in once and hit return,
then type it in again to verify. When you set a new password, make sure
that it is something secure (i.e. has at least letters and numbers in
it). Note, nothing will show up on the screen when you enter passwords.

Important: DO NOT give out your password to anyone!

#### Downloading this tutorial from github

Once you have logged in, you can download a copy of the tutorial into
your directory on the server. To retrieve the folder with this tutorial
from github directly to the server, type the following into your
terminal and hit return after each line.

``` bash
wget https://github.com/fiererlab/dada2_fiererlab/archive/master.zip
unzip master.zip
```

If there are ever updates to the tutorial on github, you can update the
contents of this folder by downloading the new version from the same
link as above.

#### Login to RStudio on the server

1.  Open a your web browser and start a new empty tab
2.  Type `microbe.colorado.edu:8787` in the address bar
3.  Use your server login credentials to log into rstudio server

If you are running it on your own computer (runs slower!):

1.  Download this tutorial from github. Go to [the
    homepage](https://github.com/fiererlab/dada2_fiererlab/dada2_fiererlab),
    and click the green “Clone or download” button. Then click “Download
    ZIP”, to save it to your computer. Unzip the file to access the
    R-script.
2.  Download the tutorial data from here
    <http://cme.colorado.edu/projects/bioinformatics-tutorials>
3.  Install idemp and cutadapt.
    - idemp can be found here: <https://github.com/yhwu/idemp>
    - cutadapt can be installed from here:
      <https://cutadapt.readthedocs.io/en/stable/installation.html>
4.  Download the dada2-formatted reference database of your choice. Link
    to download here: <https://benjjneb.github.io/dada2/training.html>

## Set up (part 2) - You are logged in to Rstudio on server (or have it open on your computer)

First open the R script in Rstudio. The R script is located in the
tutorial folder you downloaded in the first step. You can navigate to
the proper folder in Rstudio by clicking on the files tab and navigating
to the location where you downloaded the github folder. Then click
dada2_fiererlab and dada2_tutorial_16S.R to open the R script.

Now, install DADA2 & other necessary packages. If this is your first
time on Rstudio server, when you install a package you might get a
prompt asking if you want to create your own library. Answer ‘yes’ twice
in the console to continue.

| <span>                                                                                                                  |
|:------------------------------------------------------------------------------------------------------------------------|
| **WARNING:** This installation may take a long time, so only run this code if these packages are not already installed! |
| <span>                                                                                                                  |

``` r
install.packages("BiocManager")
BiocManager::install("dada2")
BiocManager::install("ShortRead")

install.packages("dplyr")
install.packages("tidyr")
install.packages("Hmisc")
install.packages("ggplot2")
install.packages("plotly")
```

Load DADA2 and required packages

``` r
library(dada2); packageVersion("dada2") # the dada2 pipeline
## [1] '1.32.0'
library(ShortRead); packageVersion("ShortRead") # dada2 depends on this
## [1] '1.62.0'
library(dplyr); packageVersion("dplyr") # for manipulating data
## [1] '1.1.4'
library(tidyr); packageVersion("tidyr") # for creating the final graph at the end of the pipeline
## [1] '1.3.1'
library(Hmisc); packageVersion("Hmisc") # for creating the final graph at the end of the pipeline
## [1] '5.1.3'
library(ggplot2); packageVersion("ggplot2") # for creating the final graph at the end of the pipeline
## [1] '3.5.1'
library(plotly); packageVersion("plotly") # enables creation of interactive graphs, especially helpful for quality plots
## [1] '4.10.4'
```

Once the packages are installed, you can check to make sure the
auxiliary software is working and set up some of the variables that you
will need along the way.

| <span>                                                                                                                                                                |
|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **NOTE:** If you are not working from microbe server, you will need to change the file paths for idemp and cutadapt to where they are stored on your computer/server. |
| <span>                                                                                                                                                                |

For this tutorial we will be working with some samples that we obtained
16S amplicon data for, from a Illumina Miseq run. The data for these
samples can be found on the CME website.
<http://cme.colorado.edu/projects/bioinformatics-tutorials>

``` r
# Set up pathway to idemp (demultiplexing tool) and test
idemp <- "/usr/bin/idemp" # CHANGE ME if not on microbe
system2(idemp) # Check that idemp is in your path and you can run shell commands from R

# Set up pathway to cutadapt (primer trimming tool) and test
cutadapt <- "/usr/local/Python27/bin/cutadapt" # CHANGE ME if not on microbe
system2(cutadapt, args = "--version") # Check by running shell command from R

# Set path to shared data folder and contents
data.fp <- "/data/shared/2019_02_20_MicrMethods_tutorial"

# List all files in shared folder to check path
list.files(data.fp)
##  [1] "barcode_demultiplex_short.txt"                                                
##  [2] "CompletedJobInfo.xml"                                                         
##  [3] "cvanderburgh"                                                                 
##  [4] "GenerateFASTQRunStatistics.xml"                                               
##  [5] "Molecular_Methods_18_515fBC_16S_Mapping_File_SHORT_vFinal_Fierer_10252018.txt"
##  [6] "other_files"                                                                  
##  [7] "RunInfo.xml"                                                                  
##  [8] "runParameters.xml"                                                            
##  [9] "Undetermined_S0_L001_I1_001.fastq.gz"                                         
## [10] "Undetermined_S0_L001_R1_001.fastq.gz"                                         
## [11] "Undetermined_S0_L001_R2_001.fastq.gz"

# Set file paths for barcodes file, map file, and fastqs
    # Barcodes may need to have 'N' on the end of each 12bp sequence for compatibility
    # Barcode file must have column 1 with barcodes and column 2 with sample IDs
barcode.fp <- file.path(data.fp, "barcode_demultiplex_short.txt") # .txt file: barcode </t> sampleID
I1.fp <- file.path(data.fp, "Undetermined_S0_L001_I1_001.fastq.gz") 
R1.fp <- file.path(data.fp, "Undetermined_S0_L001_R1_001.fastq.gz") 
R2.fp <- file.path(data.fp, "Undetermined_S0_L001_R2_001.fastq.gz") 
```

| <span>                                                                                                                                                                                                                                                                                                                                                                                                                               |
|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **NOTE:** idemp relies on having a match in length between the index file and and the barcode sequences. Since the index file usually includes a extra linker basepair (making it 13bp long), you should append the barcode sequences with “N” to make sure each is 13bp long. If you are not sure of the length of index reads, check with the sequencing center. If your index reads are 12bp long, you do NOT need to add an “N”. |
| <span>                                                                                                                                                                                                                                                                                                                                                                                                                               |

**For ITS, 18S, and Univeral Barcode Primers:** Depending on how your
sequences were run and how your primers were designed your barcodes may
need to be reverse-complemented. A good clue that you might need to
reverse complement your barcodes is if your first attempt at
demultiplexing results in the wrong number of demultiplexed files, or if
your demultiplexed files are really small. Here is a link to a handy
tool, that can help you reverse complement your barcodes:
<http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html>

Set up file paths in YOUR directory where you want data; you do not need
to create the subdirectories but they are nice to have for
organizational purposes.

``` r
project.fp <- "/data/cliffb/MicroMethods_dada2_tutorial" # CHANGE ME to project directory; don't append with a "/"

# Set up names of sub directories to stay organized
preprocess.fp <- file.path(project.fp, "01_preprocess")
    demultiplex.fp <- file.path(preprocess.fp, "demultiplexed")
    filtN.fp <- file.path(preprocess.fp, "filtN")
    trimmed.fp <- file.path(preprocess.fp, "trimmed")
filter.fp <- file.path(project.fp, "02_filter") 
table.fp <- file.path(project.fp, "03_tabletax") 
```

## Pre-processing data for dada2 - demultiplex, remove sequences with Ns, cutadapt

#### Call the demultiplexing script

Demultiplexing splits your reads out into separate files based on the
barcodes associated with each sample. If you already have demultiplexed
files, you can skip this step. Just put the demultiplexed fastq files in
your demultiplex file path “demultiplex.fp”, and make sure the forward
read file names start with R1\_ and the reverse read file names start
with R2\_.

``` r
flags <- paste("-b", barcode.fp, "-I1", I1.fp, "-R1", R1.fp, "-R2", R2.fp, "-o", demultiplex.fp) 
system2(idemp, args = flags) 

# Look at output of demultiplexing
list.files(demultiplex.fp)
##  [1] "Undetermined_S0_L001_I1_001.fastq.gz.decode"                 
##  [2] "Undetermined_S0_L001_I1_001.fastq.gz.decode.stat"            
##  [3] "Undetermined_S0_L001_R1_001.fastq.gz_ANT7.fastq.gz"          
##  [4] "Undetermined_S0_L001_R1_001.fastq.gz_ANT8.fastq.gz"          
##  [5] "Undetermined_S0_L001_R1_001.fastq.gz_BA1A5566_10_1D.fastq.gz"
##  [6] "Undetermined_S0_L001_R1_001.fastq.gz_BA1A5566_22_1C.fastq.gz"
##  [7] "Undetermined_S0_L001_R1_001.fastq.gz_BA1A5566_45_1J.fastq.gz"
##  [8] "Undetermined_S0_L001_R1_001.fastq.gz_BB1S.fastq.gz"          
##  [9] "Undetermined_S0_L001_R1_001.fastq.gz_BB1W.fastq.gz"          
## [10] "Undetermined_S0_L001_R1_001.fastq.gz_BNS1.fastq.gz"          
## [11] "Undetermined_S0_L001_R1_001.fastq.gz_BNS2.fastq.gz"          
## [12] "Undetermined_S0_L001_R1_001.fastq.gz_C2S.fastq.gz"           
## [13] "Undetermined_S0_L001_R1_001.fastq.gz_C2W.fastq.gz"           
## [14] "Undetermined_S0_L001_R1_001.fastq.gz_CM2A_10_9C.fastq.gz"    
## [15] "Undetermined_S0_L001_R1_001.fastq.gz_CM2A_22_9J.fastq.gz"    
## [16] "Undetermined_S0_L001_R1_001.fastq.gz_COL11.fastq.gz"         
## [17] "Undetermined_S0_L001_R1_001.fastq.gz_COL12.fastq.gz"         
## [18] "Undetermined_S0_L001_R1_001.fastq.gz_COL3.fastq.gz"          
## [19] "Undetermined_S0_L001_R1_001.fastq.gz_COL4.fastq.gz"          
## [20] "Undetermined_S0_L001_R1_001.fastq.gz_DT3S.fastq.gz"          
## [21] "Undetermined_S0_L001_R1_001.fastq.gz_DT3W.fastq.gz"          
## [22] "Undetermined_S0_L001_R1_001.fastq.gz_OM18_BC.fastq.gz"       
## [23] "Undetermined_S0_L001_R1_001.fastq.gz_OM18_BJ.fastq.gz"       
## [24] "Undetermined_S0_L001_R1_001.fastq.gz_unsigned.fastq.gz"      
## [25] "Undetermined_S0_L001_R1_001.fastq.gz_WAB105_22_6J.fastq.gz"  
## [26] "Undetermined_S0_L001_R1_001.fastq.gz_WAB105_45_6D.fastq.gz"  
## [27] "Undetermined_S0_L001_R1_001.fastq.gz_WAB188_10_4D.fastq.gz"  
## [28] "Undetermined_S0_L001_R1_001.fastq.gz_WAB71_45_3D.fastq.gz"   
## [29] "Undetermined_S0_L001_R2_001.fastq.gz_ANT7.fastq.gz"          
## [30] "Undetermined_S0_L001_R2_001.fastq.gz_ANT8.fastq.gz"          
## [31] "Undetermined_S0_L001_R2_001.fastq.gz_BA1A5566_10_1D.fastq.gz"
## [32] "Undetermined_S0_L001_R2_001.fastq.gz_BA1A5566_22_1C.fastq.gz"
## [33] "Undetermined_S0_L001_R2_001.fastq.gz_BA1A5566_45_1J.fastq.gz"
## [34] "Undetermined_S0_L001_R2_001.fastq.gz_BB1S.fastq.gz"          
## [35] "Undetermined_S0_L001_R2_001.fastq.gz_BB1W.fastq.gz"          
## [36] "Undetermined_S0_L001_R2_001.fastq.gz_BNS1.fastq.gz"          
## [37] "Undetermined_S0_L001_R2_001.fastq.gz_BNS2.fastq.gz"          
## [38] "Undetermined_S0_L001_R2_001.fastq.gz_C2S.fastq.gz"           
## [39] "Undetermined_S0_L001_R2_001.fastq.gz_C2W.fastq.gz"           
## [40] "Undetermined_S0_L001_R2_001.fastq.gz_CM2A_10_9C.fastq.gz"    
## [41] "Undetermined_S0_L001_R2_001.fastq.gz_CM2A_22_9J.fastq.gz"    
## [42] "Undetermined_S0_L001_R2_001.fastq.gz_COL11.fastq.gz"         
## [43] "Undetermined_S0_L001_R2_001.fastq.gz_COL12.fastq.gz"         
## [44] "Undetermined_S0_L001_R2_001.fastq.gz_COL3.fastq.gz"          
## [45] "Undetermined_S0_L001_R2_001.fastq.gz_COL4.fastq.gz"          
## [46] "Undetermined_S0_L001_R2_001.fastq.gz_DT3S.fastq.gz"          
## [47] "Undetermined_S0_L001_R2_001.fastq.gz_DT3W.fastq.gz"          
## [48] "Undetermined_S0_L001_R2_001.fastq.gz_OM18_BC.fastq.gz"       
## [49] "Undetermined_S0_L001_R2_001.fastq.gz_OM18_BJ.fastq.gz"       
## [50] "Undetermined_S0_L001_R2_001.fastq.gz_unsigned.fastq.gz"      
## [51] "Undetermined_S0_L001_R2_001.fastq.gz_WAB105_22_6J.fastq.gz"  
## [52] "Undetermined_S0_L001_R2_001.fastq.gz_WAB105_45_6D.fastq.gz"  
## [53] "Undetermined_S0_L001_R2_001.fastq.gz_WAB188_10_4D.fastq.gz"  
## [54] "Undetermined_S0_L001_R2_001.fastq.gz_WAB71_45_3D.fastq.gz"
```

| <span>                                                                                                                                                                                                                                                                                           |
|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **WARNING:** The demultiplexing step may take a while. If it takes too long you can safely close RStudio on the server and the demultiplexing will run in the background. You should be able to resume the pipeline after demultiplexing is complete by logging back into RStudio on the server. |
| <span>                                                                                                                                                                                                                                                                                           |

#### Clean up the output from idemp

``` r
# Change names of unassignable reads so they are not included in downstream processing
unassigned_1 <- paste0("mv", " ", demultiplex.fp, "/Undetermined_S0_L001_R1_001.fastq.gz_unsigned.fastq.gz",
                       " ", demultiplex.fp, "/Unassigned_reads1.fastq.gz")
unassigned_2 <- paste0("mv", " ", demultiplex.fp, "/Undetermined_S0_L001_R2_001.fastq.gz_unsigned.fastq.gz", 
                       " ", demultiplex.fp, "/Unassigned_reads2.fastq.gz")
system(unassigned_1)
system(unassigned_2)

# Rename files - use gsub to get names in order!
R1_names <- gsub(paste0(demultiplex.fp, "/Undetermined_S0_L001_R1_001.fastq.gz_"), "", 
                 list.files(demultiplex.fp, pattern="Undetermined_S0_L001_R1_001", full.names = TRUE))
file.rename(list.files(demultiplex.fp, pattern="Undetermined_S0_L001_R1_001", full.names = TRUE), 
            paste0(demultiplex.fp, "/R1_", R1_names))
##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
## [16] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

R2_names <- gsub(paste0(demultiplex.fp, "/Undetermined_S0_L001_R2_001.fastq.gz_"), "", 
                 list.files(demultiplex.fp, pattern="Undetermined_S0_L001_R2_001", full.names = TRUE))
file.rename(list.files(demultiplex.fp, pattern="Undetermined_S0_L001_R2_001", full.names = TRUE),
            paste0(demultiplex.fp, "/R2_", R2_names))
##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
## [16] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

# Get full paths for all files and save them for downstream analyses
# Forward and reverse fastq filenames have format: 
fnFs <- sort(list.files(demultiplex.fp, pattern="R1_", full.names = TRUE))
fnRs <- sort(list.files(demultiplex.fp, pattern="R2_", full.names = TRUE))
```

#### Pre-filter to remove sequence reads with Ns

Ambiguous bases will make it hard for cutadapt to find short primer
sequences in the reads. To solve this problem, we will remove sequences
with ambiguous bases (Ns)

``` r
# Name the N-filtered files to put them in filtN/ subdirectory
fnFs.filtN <- file.path(preprocess.fp, "filtN", basename(fnFs))
fnRs.filtN <- file.path(preprocess.fp, "filtN", basename(fnRs))

# Filter Ns from reads and put them into the filtN directory
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) 
# CHANGE multithread to FALSE on Windows (here and elsewhere in the program)
```

| <span>                                                                                                                                                                                       |
|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Note:** The `multithread = TRUE` setting can sometimes generate an error (names not equal). If this occurs, try rerunning the function. The error normally does not occur the second time. |
| <span>                                                                                                                                                                                       |

#### Prepare the primers sequences and custom functions for analyzing the results from cutadapt

Assign the primers you used to “FWD” and “REV” below. Note primers
should be not be reverse complemented ahead of time. Our tutorial data
uses 515f and 806br. Those are the primers below. Change if you
sequenced with other primers. See the Earth Microbiome project for
standard 16S, 18S, and ITS primer sequences.

**For ITS data:** `CTTGGTCATTTAGAGGAAGTAA` is the ITS forward primer
sequence (ITS1F) and `GCTGCGTTCTTCATCGATGC` is ITS reverse primer
sequence (ITS2)

``` r
# Set up the primer sequences to pass along to cutadapt
FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE ME # this is 515f
REV <- "GGACTACNVGGGTWTCTAAT"  ## CHANGE ME # this is 806Br

# Write a function that creates a list of all orientations of the primers
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
                 RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}

# Save the primer orientations to pass to cutadapt
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
##               Forward            Complement               Reverse 
## "GTGYCAGCMGCCGCGGTAA" "CACRGTCGKCGGCGCCATT" "AATGGCGCCGMCGACYGTG" 
##               RevComp 
## "TTACCGCGGCKGCTGRCAC"
REV.orients
##                Forward             Complement                Reverse 
## "GGACTACNVGGGTWTCTAAT" "CCTGATGNBCCCAWAGATTA" "TAATCTWTGGGVNCATCAGG" 
##                RevComp 
## "ATTAGAWACCCBNGTAGTCC"

# Write a function that counts how many time primers appear in a sequence
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}
```

Before running cutadapt, we will look at primer detection for the first
couple samples, as a check. There may be some primers here, we will
remove them below using cutadapt. If your first couple samples happen to
be blanks, change the number to a real sample.

``` r
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
##                  Forward Complement Reverse RevComp
## FWD.ForwardReads       2          0       0       0
## FWD.ReverseReads       0          0       0      35
## REV.ForwardReads       0          0       0      56
## REV.ReverseReads       0          0       0       0
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[2]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[2]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[2]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[2]]))
##                  Forward Complement Reverse RevComp
## FWD.ForwardReads       0          0       0       0
## FWD.ReverseReads       0          0       0      39
## REV.ForwardReads       0          0       0      62
## REV.ReverseReads       0          0       0       0
```

#### Remove primers with cutadapt and assess the output

``` r
# Create directory to hold the output from cutadapt
if (!dir.exists(trimmed.fp)) dir.create(trimmed.fp)
fnFs.cut <- file.path(trimmed.fp, basename(fnFs.filtN))
fnRs.cut <- file.path(trimmed.fp, basename(fnRs.filtN))

# Save the reverse complements of the primers to variables
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

##  Create the cutadapt flags ##
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC, "--minimum-length 50") 

# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC, "--minimum-length 50") 

# Run Cutadapt
for (i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# As a sanity check, we will check for primers in the first cutadapt-ed sample:
    ## should all be zero!
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
##                  Forward Complement Reverse RevComp
## FWD.ForwardReads       0          0       0       0
## FWD.ReverseReads       0          0       0       0
## REV.ForwardReads       0          0       0       0
## REV.ReverseReads       0          0       0       0
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[2]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[2]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[2]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[2]]))
##                  Forward Complement Reverse RevComp
## FWD.ForwardReads       0          0       0       0
## FWD.ReverseReads       0          0       0       0
## REV.ForwardReads       0          0       0       0
## REV.ReverseReads       0          0       0       0
```

# Now start DADA2 pipeline

``` r
# Put filtered reads into separate sub-directories for big data workflow
dir.create(filter.fp)
    subF.fp <- file.path(filter.fp, "preprocessed_F") 
    subR.fp <- file.path(filter.fp, "preprocessed_R") 
dir.create(subF.fp)
dir.create(subR.fp)

# Move R1 and R2 from trimmed to separate forward/reverse sub-directories
fnFs.Q <- file.path(subF.fp,  basename(fnFs.filtN)) 
fnRs.Q <- file.path(subR.fp,  basename(fnRs.filtN))
file.rename(from = fnFs.cut, to = fnFs.Q)
##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
## [16] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
file.rename(from = fnRs.cut, to = fnRs.Q)
##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
## [16] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE

# File parsing; create file names and make sure that forward and reverse files match
filtpathF <- file.path(subF.fp, "filtered") # files go into preprocessed_F/filtered/
filtpathR <- file.path(subR.fp, "filtered") # files go into preprocessed_R/filtered/
fastqFs <- sort(list.files(subF.fp, pattern="fastq.gz"))
fastqRs <- sort(list.files(subR.fp, pattern="fastq.gz"))

# Check to make sure your forward and reverse read files are in order.
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
```

### 1. FILTER AND TRIM FOR QUALITY

Before choosing sequence variants, we want to trim reads where their
quality scores begin to drop (the `truncLen` and `truncQ` values) and
remove any low-quality reads that are left over after we have finished
trimming (the `maxEE` value).

**You will want to change this depending on run chemistry and quality:**
For 2x250 bp runs you can try `truncLen=c(240,160)` (as per the [dada2
tutorial](https://benjjneb.github.io/dada2/tutorial.html#inspect-read-quality-profiles))
if your reverse reads drop off in quality. Or you may want to choose a
higher value, for example, `truncLen=c(240,200)`, if they do not. In
`truncLen=c(xxx,yyy)`, `xxx` refers to the forward read truncation
length, `yyy` refers to the reverse read truncation length.

**For ITS data:** Due to the expected variable read lengths in ITS data
you should run this command without the `trunclen` parameter. See here
for more information and appropriate parameters for ITS data:
<https://benjjneb.github.io/dada2/ITS_workflow.html>.

*From dada2 tutorial:* \>If there is only one part of any amplicon
bioinformatics workflow on which you spend time considering the
parameters, it should be filtering! The parameters … are not set in
stone, and should be changed if they don’t work for your data. If too
few reads are passing the filter, increase maxEE and/or reduce truncQ.
If quality drops sharply at the end of your reads, reduce truncLen. If
your reads are high quality and you want to reduce computation time in
the sample inference step, reduce maxEE.

#### Inspect read quality profiles

It’s important to get a feel for the quality of the data that we are
using. To do this, we will plot the quality of some of the samples.

*From the dada2 tutorial:* \>In gray-scale is a heat map of the
frequency of each quality score at each base position. The median
quality score at each position is shown by the green line, and the
quartiles of the quality score distribution by the orange lines. The red
line shows the scaled proportion of reads that extend to at least that
position (this is more useful for other sequencing technologies, as
Illumina reads are typically all the same lenghth, hence the flat red
line).

``` r
# If the number of samples is 20 or less, plot them all, otherwise, just plot 20 randomly selected samples
if( length(fastqFs) <= 20) {
  plotQualityProfile(paste0(subF.fp, "/", fastqFs))
  plotQualityProfile(paste0(subR.fp, "/", fastqRs))
} else {
  rand_samples <- sample(size = 20, 1:length(fastqFs)) # grab 20 random samples to plot
  fwd_qual_plots <- plotQualityProfile(paste0(subF.fp, "/", fastqFs[rand_samples]))
  rev_qual_plots <- plotQualityProfile(paste0(subR.fp, "/", fastqRs[rand_samples]))
}

fwd_qual_plots
```

<img src="figure/unnamed-chunk-11-1.png" width="98%" height="98%" />

``` r
rev_qual_plots
```

<img src="figure/unnamed-chunk-11-2.png" width="98%" height="98%" />

``` r
# Optional: To make these quality plots interactive, call the plots through plotly
ggplotly(fwd_qual_plots)
ggplotly(rev_qual_plots)
```

``` r
# write plots to disk
saveRDS(fwd_qual_plots, paste0(filter.fp, "/fwd_qual_plots.rds"))
saveRDS(rev_qual_plots, paste0(filter.fp, "/rev_qual_plots.rds"))

ggsave(plot = fwd_qual_plots, filename = paste0(filter.fp, "/fwd_qual_plots.png"), 
       width = 10, height = 10, dpi = "retina")
ggsave(plot = rev_qual_plots, filename = paste0(filter.fp, "/rev_qual_plots.png"), 
       width = 10, height = 10, dpi = "retina")
```

#### Filter the data

| <span>                                                                                                                                                                                                                                                                                                                        |
|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **WARNING:** THESE PARAMETERS ARE NOT OPTIMAL FOR ALL DATASETS. Make sure you determine the trim and filtering parameters for your data. The following settings are generally appropriate for MiSeq runs that are 2x150 bp. These are the recommended default parameters from the dada2 pipeline. See above for more details. |
| <span>                                                                                                                                                                                                                                                                                                                        |

``` r
filt_out <- filterAndTrim(fwd = file.path(subF.fp, fastqFs), 
                          filt = file.path(filtpathF, fastqFs),
                          rev = file.path(subR.fp, fastqRs), 
                          filt.rev = file.path(filtpathR, fastqRs),
                          truncLen = c(150,140), 
                          maxEE = c(2,2), 
                          truncQ = 2, 
                          maxN = 0,
                          rm.phix = TRUE,
                          compress = TRUE, 
                          verbose = TRUE, 
                          multithread = TRUE)

# look at how many reads were kept
filt_out
##                            reads.in reads.out
## R1_ANT7.fastq.gz              40165     38746
## R1_ANT8.fastq.gz              25526     24723
## R1_BA1A5566_10_1D.fastq.gz    50515     48763
## R1_BA1A5566_22_1C.fastq.gz    42238     40623
## R1_BA1A5566_45_1J.fastq.gz    31057     29830
## R1_BB1S.fastq.gz              34353     33063
## R1_BB1W.fastq.gz              39824     38337
## R1_BNS1.fastq.gz               3792      3684
## R1_BNS2.fastq.gz              16915     16362
## R1_C2S.fastq.gz               38469     37037
## R1_C2W.fastq.gz               41979     40359
## R1_CM2A_10_9C.fastq.gz        34309     33206
## R1_CM2A_22_9J.fastq.gz        31079     29844
## R1_COL11.fastq.gz             31560     30423
## R1_COL12.fastq.gz             24413     23608
## R1_COL3.fastq.gz              40397     38933
## R1_COL4.fastq.gz              32665     31406
## R1_DT3S.fastq.gz              36902     35316
## R1_DT3W.fastq.gz              44790     42982
## R1_OM18_BC.fastq.gz            4803      4603
## R1_OM18_BJ.fastq.gz             452       391
## R1_WAB105_22_6J.fastq.gz      44993     43423
## R1_WAB105_45_6D.fastq.gz      37348     36026
## R1_WAB188_10_4D.fastq.gz      31134     29793
## R1_WAB71_45_3D.fastq.gz       43340     41726

# summary of samples in filt_out by percentage
filt_out %>% 
  data.frame() %>% 
  mutate(Samples = rownames(.),
         percent_kept = 100*(reads.out/reads.in)) %>%
  select(Samples, everything()) %>%
  summarise(min_remaining = paste0(round(min(percent_kept), 2), "%"), 
            median_remaining = paste0(round(median(percent_kept), 2), "%"),
            mean_remaining = paste0(round(mean(percent_kept), 2), "%"), 
            max_remaining = paste0(round(max(percent_kept), 2), "%"))
##   min_remaining median_remaining mean_remaining max_remaining
## 1         86.5%           96.28%         95.93%        97.15%

# If very few reads made it through, go back and tweak your parameters.
```

Plot the quality of the filtered fastq files.

``` r
# figure out which samples, if any, have been filtered out
if( length(fastqFs) <= 20) {
  remaining_samplesF <-  fastqFs[
  which(fastqFs %in% list.files(filtpathF))] # keep only samples that haven't been filtered out
remaining_samplesR <-  fastqRs[
  which(fastqRs %in% list.files(filtpathR))] # keep only samples that haven't been filtered out
fwd_qual_plots_filt <- plotQualityProfile(paste0(filtpathF, "/", remaining_samplesF))
rev_qual_plots_filt <- plotQualityProfile(paste0(filtpathR, "/", remaining_samplesR))
} else {
remaining_samplesF <-  fastqFs[rand_samples][
  which(fastqFs[rand_samples] %in% list.files(filtpathF))] # keep only samples that haven't been filtered out
remaining_samplesR <-  fastqRs[rand_samples][
  which(fastqRs[rand_samples] %in% list.files(filtpathR))] # keep only samples that haven't been filtered out
fwd_qual_plots_filt <- plotQualityProfile(paste0(filtpathF, "/", remaining_samplesF))
rev_qual_plots_filt <- plotQualityProfile(paste0(filtpathR, "/", remaining_samplesR))
}

fwd_qual_plots_filt
```

<img src="figure/unnamed-chunk-14-1.png" width="98%" height="98%" />

``` r
rev_qual_plots_filt
```

<img src="figure/unnamed-chunk-14-2.png" width="98%" height="98%" />

``` r

# write plots to disk
saveRDS(fwd_qual_plots_filt, paste0(filter.fp, "/fwd_qual_plots_filt.rds"))
saveRDS(rev_qual_plots_filt, paste0(filter.fp, "/rev_qual_plots_filt.rds"))

ggsave(plot = fwd_qual_plots_filt, filename = paste0(filter.fp, "/fwd_qual_plots_filt.png"), 
       width = 10, height = 10, dpi = "retina")
ggsave(plot = rev_qual_plots_filt, filename = paste0(filter.fp, "/rev_qual_plots_filt.png"), 
       width = 10, height = 10, dpi = "retina")
```

### 2. INFER sequence variants

In this part of the pipeline dada2 will learn to distinguish error from
biological differences using a subset of our data as a training set.
After it understands the error rates, we will reduce the size of the
dataset by combining all identical sequence reads into “unique
sequences”. Then, using the dereplicated data and error rates, dada2
will infer the amplicon sequence variants (ASVs) in our data. These are
essentially OTUs at 100% similarity. Finally, we will merge the
corresponding forward and reverse reads to create a list of the fully
denoised sequences and create a sequence table from the result.

#### Housekeeping step - set up and verify the file names for the output:

``` r
# File parsing
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)

# Sample names in order
sample.names <- substring(basename(filtFs), regexpr("_", basename(filtFs)) + 1) # doesn't drop fastq.gz
sample.names <- gsub(".fastq.gz", "", sample.names)
sample.namesR <- substring(basename(filtRs), regexpr("_", basename(filtRs)) + 1) # doesn't drop fastq.gz
sample.namesR <- gsub(".fastq.gz", "", sample.namesR)

# Double check
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

#### Learn the error rates

``` r
set.seed(100) # set seed to ensure that randomized steps are repeatable

# Learn forward error rates (Notes: randomize default is FALSE)
errF <- learnErrors(filtFs, nbases = 1e8, multithread = TRUE, randomize = TRUE)
## 101455950 total bases in 676373 reads from 22 samples will be used for learning the error rates.

# Learn reverse error rates
set.seed(100) # reset seed to ensure that randomized steps are repeatable
errR <- learnErrors(filtRs, nbases = 1e8, multithread = TRUE, randomize = TRUE)
## 102407340 total bases in 731481 reads from 24 samples will be used for learning the error rates.
```

#### Plot Error Rates

We want to make sure that the machine learning algorithm is learning the
error rates properly. In the plots below, the red line represents what
we should expect the learned error rates to look like for each of the 16
possible base transitions (A-\>A, A-\>C, A-\>G, etc.) and the black line
and grey dots represent what the observed error rates are. If the black
line and the red lines are very far off from each other, it may be a
good idea to increase the `nbases` parameter. This allows the machine
learning algorithm to train on a larger portion of your data and may
help improve the fit.

``` r
errF_plot <- plotErrors(errF, nominalQ = TRUE)
errR_plot <- plotErrors(errR, nominalQ = TRUE)

errF_plot
```

<img src="figure/unnamed-chunk-17-1.png" width="98%" height="98%" />

``` r
errR_plot
```

<img src="figure/unnamed-chunk-17-2.png" width="98%" height="98%" />

``` r
# write error objects and plots to disk
saveRDS(errF, paste0(filtpathF, "/errF.rds"))
saveRDS(errR, paste0(filtpathR, "/errR.rds"))
saveRDS(errF_plot, paste0(filtpathF, "/errF_plot.rds"))
saveRDS(errR_plot, paste0(filtpathR, "/errR_plot.rds"))
```

#### Dereplication, sequence inference, and merging of paired-end reads

In this part of the pipeline, dada2 will make decisions about assigning
sequences to ASVs (called “sequence inference”). There is a major
parameter option in the core function dada() that changes how samples
are handled during sequence inference. The parameter `pool =` can be set
to: `pool = FALSE` (default), `pool = TRUE`, or `pool = psuedo`. For
details on parameter choice, please see below, and further information
on this blogpost
<http://fiererlab.org/2020/02/17/whats-in-a-number-estimating-microbial-richness-using-dada2/>,
and explanation on the dada2 tutorial
<https://benjjneb.github.io/dada2/pool.html>. Note that some amplicons
like ITS have variable read lengths which in some cases means that reads
may not overlap and cannot be merged. In such cases you may want to opt
to just use the forward or reverse reads rather than merged reads. For
more information see
<https://benjjneb.github.io/dada2/ITS_workflow.html>.

**Details**  
`pool = FALSE`: Sequence information is not shared between samples. Fast
processing time, less sensitivity to rare taxa.  
`pool = psuedo`: Sequence information is shared in a separate “prior”
step. Intermediate processing time, intermediate sensitivity to rare
taxa.  
`pool = TRUE`: Sequence information from all samples is pooled together.
Slow processing time, most sensitivity to rare taxa.

#### Default: SAMPLES NOT POOLED

For simple communities or when you do not need high sensitivity for rare
taxa

``` r
# make lists to hold the loop output
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
ddF <- vector("list", length(sample.names))
names(ddF) <- sample.names
ddR <- vector("list", length(sample.names))
names(ddR) <- sample.names

# For each sample, get a list of merged and denoised sequences
for(sam in sample.names) {
    cat("Processing:", sam, "\n")
    # Dereplicate forward reads
    derepF <- derepFastq(filtFs[[sam]])
    # Infer sequences for forward reads
    dadaF <- dada(derepF, err = errF, multithread = TRUE)
    ddF[[sam]] <- dadaF
    # Dereplicate reverse reads
    derepR <- derepFastq(filtRs[[sam]])
    # Infer sequences for reverse reads
    dadaR <- dada(derepR, err = errR, multithread = TRUE)
    ddR[[sam]] <- dadaR
    # Merge reads together
    merger <- mergePairs(ddF[[sam]], derepF, ddR[[sam]], derepR)
    mergers[[sam]] <- merger
}
## Processing: ANT7 
## Sample 1 - 38746 reads in 16505 unique sequences.
## Sample 1 - 38746 reads in 14810 unique sequences.
## Processing: ANT8 
## Sample 1 - 24723 reads in 10738 unique sequences.
## Sample 1 - 24723 reads in 9620 unique sequences.
## Processing: BA1A5566_10_1D 
## Sample 1 - 48763 reads in 12161 unique sequences.
## Sample 1 - 48763 reads in 11339 unique sequences.
## Processing: BA1A5566_22_1C 
## Sample 1 - 40623 reads in 9985 unique sequences.
## Sample 1 - 40623 reads in 8732 unique sequences.
## Processing: BA1A5566_45_1J 
## Sample 1 - 29830 reads in 7601 unique sequences.
## Sample 1 - 29830 reads in 6225 unique sequences.
## Processing: BB1S 
## Sample 1 - 33063 reads in 11503 unique sequences.
## Sample 1 - 33063 reads in 9938 unique sequences.
## Processing: BB1W 
## Sample 1 - 38337 reads in 10213 unique sequences.
## Sample 1 - 38337 reads in 7858 unique sequences.
## Processing: BNS1 
## Sample 1 - 3684 reads in 886 unique sequences.
## Sample 1 - 3684 reads in 809 unique sequences.
## Processing: BNS2 
## Sample 1 - 16362 reads in 3931 unique sequences.
## Sample 1 - 16362 reads in 3576 unique sequences.
## Processing: C2S 
## Sample 1 - 37037 reads in 14537 unique sequences.
## Sample 1 - 37037 reads in 12349 unique sequences.
## Processing: C2W 
## Sample 1 - 40359 reads in 9705 unique sequences.
## Sample 1 - 40359 reads in 6851 unique sequences.
## Processing: CM2A_10_9C 
## Sample 1 - 33206 reads in 5265 unique sequences.
## Sample 1 - 33206 reads in 4317 unique sequences.
## Processing: CM2A_22_9J 
## Sample 1 - 29844 reads in 5688 unique sequences.
## Sample 1 - 29844 reads in 4663 unique sequences.
## Processing: COL11 
## Sample 1 - 30423 reads in 15072 unique sequences.
## Sample 1 - 30423 reads in 13813 unique sequences.
## Processing: COL12 
## Sample 1 - 23608 reads in 11851 unique sequences.
## Sample 1 - 23608 reads in 10657 unique sequences.
## Processing: COL3 
## Sample 1 - 38933 reads in 19379 unique sequences.
## Sample 1 - 38933 reads in 17624 unique sequences.
## Processing: COL4 
## Sample 1 - 31406 reads in 15661 unique sequences.
## Sample 1 - 31406 reads in 14126 unique sequences.
## Processing: DT3S 
## Sample 1 - 35316 reads in 12977 unique sequences.
## Sample 1 - 35316 reads in 11359 unique sequences.
## Processing: DT3W 
## Sample 1 - 42982 reads in 11632 unique sequences.
## Sample 1 - 42982 reads in 8518 unique sequences.
## Processing: OM18_BC 
## Sample 1 - 4603 reads in 1285 unique sequences.
## Sample 1 - 4603 reads in 1101 unique sequences.
## Processing: OM18_BJ 
## Sample 1 - 391 reads in 130 unique sequences.
## Sample 1 - 391 reads in 112 unique sequences.
## Processing: WAB105_22_6J 
## Sample 1 - 43423 reads in 9283 unique sequences.
## Sample 1 - 43423 reads in 8167 unique sequences.
## Processing: WAB105_45_6D 
## Sample 1 - 36026 reads in 7663 unique sequences.
## Sample 1 - 36026 reads in 7031 unique sequences.
## Processing: WAB188_10_4D 
## Sample 1 - 29793 reads in 6088 unique sequences.
## Sample 1 - 29793 reads in 4871 unique sequences.
## Processing: WAB71_45_3D 
## Sample 1 - 41726 reads in 10190 unique sequences.
## Sample 1 - 41726 reads in 9586 unique sequences.

rm(derepF); rm(derepR)
```

#### Alternative: SAMPLES POOLED

For complex communities when you want to preserve rare taxa. Note: This
adds substantially more computational time!! Learn more about pooling
here: <https://benjjneb.github.io/dada2/pool.html>. The DADA2 author
Benjamin Callahan states there that: “The large majority of ASVs, and
the vast majority of reads (\>99%) are commonly assigned between these
methods. Either method provides an accurate reconstruction of your
sampled communities.” alternative: swap `pool = TRUE` with
`pool = "pseudo"`

``` r
# same steps, but not in loop because samples are now pooled

# Dereplicate forward reads
derepF.p <- derepFastq(filtFs)
names(derepF.p) <- sample.names
# Infer sequences for forward reads
dadaF.p <- dada(derepF.p, err = errF, multithread = TRUE, pool = TRUE)
names(dadaF.p) <- sample.names

# Dereplicate reverse reads
derepR.p <- derepFastq(filtRs)
names(derepR.p) <- sample.names
# Infer sequences for reverse reads
dadaR.p <- dada(derepR.p, err = errR, multithread = TRUE, pool = TRUE)
names(dadaR.p) <- sample.names

# Merge reads together
mergers <- mergePairs(dadaF.p, derepF.p, dadaR.p, derepR.p)
```

#### Construct sequence table

``` r
seqtab <- makeSequenceTable(mergers)

# Save table as an r data object file
dir.create(table.fp)
saveRDS(seqtab, paste0(table.fp, "/seqtab.rds"))
```

#### Optional: Merge sequence tables

If you have multiple sequencing runs, if you have processed them each to
this point, you can now load all of your different seqtab.rds files in
here and merge them. Then you can remove chimeras and assign taxonomy
just on the one merged table.

``` r
# seqtab1 <- readRDS(paste0(table.fp, "/seqtab.rds"))
# seqtab2 <- readRDS("path/to/other/seqtab.rds")
# seqtab <- mergeSequenceTables(seqtab1, seqtab2)
```

### 3. REMOVE Chimeras and ASSIGN Taxonomy

Although dada2 has searched for indel errors and substitutions, there
may still be chimeric sequences in our dataset (sequences that are
derived from forward and reverse sequences from two different organisms
becoming fused together during PCR and/or sequencing). To identify
chimeras, we will search for rare sequence variants that can be
reconstructed by combining left-hand and right-hand segments from two
more abundant “parent” sequences. After removing chimeras, we will use a
taxonomy database to train a naive Bayesian classifier-algorithm to
assign names to our sequence variants.

For the tutorial 16S, we will assign taxonomy with SILVA db v138.1, but
you might want to use other databases for your data. This one includes
species but be wary of species level assignments. It’s always a good
idea to BLAST your top ASVs of interest to confirm the taxonomy. Other
options for 16S include Greengenes, GTDB, and RDP. Below are paths to
some of the databases we use often. If you want to use SILVA for 18S
instead of PR2, make sure you use the one that is suitable for 18S,
which is different than the 16S one and is from release 132. Databases
are periodically updated so make sure you check for new releases and use
the latest versions. (If you are on your own computer you can download
the database you need from this link
<https://benjjneb.github.io/dada2/training.html>:)

- 16S bacteria and archaea (SILVA db):
  /db_files/dada2/silva_nr99_v138.1_train_set.fa

- ITS fungi (UNITE db):
  /db_files/dada2/sh_general_release_dynamic_25.07.2023.fasta

- 18S protists (PR2 db):
  /db_files/dada2/pr2_version_5.0.0_SSU_dada2.fasta.gz

``` r
# Read in RDS 
st.all <- readRDS(paste0(table.fp, "/seqtab.rds"))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(st.all, 
                                    method = "consensus", 
                                    multithread = TRUE)

# Print percentage of our sequences that were not chimeric.
100*sum(seqtab.nochim)/sum(seqtab)
## [1] 98.34444

# Print starting number of ASVs, number of chimeras, and non-chimeric ASVs
ncol(st.all) # starting number of ASVs
## [1] 6283
ncol(st.all) - ncol(seqtab.nochim) # number of chimeras. there are removed.
## [1] 443
ncol(seqtab.nochim) # number of non-chimeric ASVs
## [1] 5840

# Assign taxonomy
# Note: assignTaxonomy implements the RDP Naive Bayesian Classifier algorithm described in Wang et al. Applied and Environmental Microbiology 2007, with kmer size 8 and 100 bootstrap replicates.
tax <- assignTaxonomy(seqs = seqtab.nochim, 
                      refFasta = "/db_files/dada2/silva_nr99_v138.1_train_set.fa", 
                      minBoot = 50,
                      tryRC = TRUE,
                      outputBootstraps = FALSE,
                      taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
                      multithread = TRUE,
                      verbose = FALSE)

# If you want to add species level, you must use this separate function that uses exact matching
tax <- addSpecies(taxtab = tax, 
                  refFasta = "/db_files/dada2/silva_species_assignment_v138.1.fa.gz",
                  allowMultiple = FALSE,
                  tryRC = TRUE,
                  n = 1e5,
                  verbose = FALSE)

# Write results to disk
saveRDS(seqtab.nochim, paste0(table.fp, "/seqtab_final.rds"))
saveRDS(tax, paste0(table.fp, "/tax_final.rds"))
```

### 4. Optional - FORMAT OUTPUT to obtain ASV IDs and repset, and input for mctoolsr and phyloseq

For convenience, we will now rename our ASVs with numbers, output our
results as a traditional taxa table, and create a matrix with the
representative sequences for each ASV. Having the actual sequences for
each ASV is important if you want to build a phylogenetic tree or do
additional taxonomic investigation.

``` r
# Flip table
seqtab.t <- as.data.frame(t(seqtab.nochim))

# Pull out ASV repset
rep_set_ASVs <- as.data.frame(rownames(seqtab.t))
rep_set_ASVs <- mutate(rep_set_ASVs, ASV_ID = 1:n())
rep_set_ASVs$ASV_ID <- sub("^", "ASV_", rep_set_ASVs$ASV_ID)
rep_set_ASVs$ASV <- rep_set_ASVs$`rownames(seqtab.t)` 
rep_set_ASVs$`rownames(seqtab.t)` <- NULL

# Add ASV numbers to table
rownames(seqtab.t) <- rep_set_ASVs$ASV_ID

# Add ASV numbers to taxonomy
taxonomy <- as.data.frame(tax)
taxonomy$ASV <- as.factor(rownames(taxonomy))
taxonomy <- merge(rep_set_ASVs, taxonomy, by = "ASV")
rownames(taxonomy) <- taxonomy$ASV_ID
# Note: If you used a database that only goes to genus level, delete the "Species" level here.
# Note: In mctoolsr, these will correspond to taxonomy1-8 (or 1-7 for genus level)
taxonomy_for_mctoolsr <- unite(taxonomy, 
                               "taxonomy", 
                                c("Kingdom", "Phylum", "Class", "Order", 
                                  "Family", "Genus", "Species", "ASV_ID"),
                                sep = ";")

# Write repset to fasta file
# create a function that writes fasta sequences
writeRepSetFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

# Arrange the taxonomy dataframe for the writeRepSetFasta function
# Note: Remove Species according to your database.
taxonomy_for_fasta <- taxonomy %>%
  unite("TaxString", 
        c("Kingdom", "Phylum", "Class", "Order", 
          "Family", "Genus", "Species", "ASV_ID"), 
        sep = ";", 
        remove = FALSE) %>%
  unite("name", 
        c("ASV_ID", "TaxString"), 
        sep = " ", 
        remove = TRUE) %>%
  select(ASV, name) %>%
  rename(seq = ASV)

# write fasta file
writeRepSetFasta(taxonomy_for_fasta, paste0(table.fp, "/repset.fasta"))

# Merge taxonomy and table
seqtab_wTax <- merge(seqtab.t, taxonomy_for_mctoolsr, by = 0)
seqtab_wTax$ASV <- NULL 

# Set name of table in mctoolsr format and save
out_fp <- paste0(table.fp, "/seqtab_wTax_mctoolsr.txt")
names(seqtab_wTax)[1] = "#ASV_ID"
write("#Exported for mctoolsr", out_fp)
suppressWarnings(write.table(seqtab_wTax, out_fp, sep = "\t", row.names = FALSE, append = TRUE))

# Also export files as .txt
write.table(seqtab.t, file = paste0(table.fp, "/seqtab_final.txt"),
            sep = "\t", row.names = TRUE, col.names = NA)
write.table(tax, file = paste0(table.fp, "/tax_final.txt"), 
            sep = "\t", row.names = TRUE, col.names = NA)


# Optional: Output for phyloseq, if you use phyloseq instead of mctoolsr
# library(BiocManager)
# BiocManager::install("phyloseq")
# library(phyloseq)
# Use seqtab.nochim for the otu_table
# Use the tax object from assignTaxonomy() for the tax_table
# Load in your sample metadata and call it samdf. Sample IDs much match the row names of seqtab.nochim
# ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE), 
#                sample_data(samdf), 
#                tax_table(tax))
```

### Summary of output files:

1.  seqtab_final.txt - A tab-delimited sequence-by-sample (i.e. OTU)
    table
2.  tax_final.txt - a tab-delimited file showing the relationship
    between ASVs, ASV IDs, and their taxonomy
3.  seqtab_wTax_mctoolsr.txt - a tab-delimited file with ASVs as rows,
    samples as columns and the final column showing the taxonomy of the
    ASV ID
4.  repset.fasta - a fasta file with the representative sequence of each
    ASV. Fasta headers are the ASV ID and taxonomy string.

### 5. Summary of reads throughout pipeline

Here we track the reads throughout the pipeline to see if any step is
resulting in a greater-than-expected loss of reads. If a step is showing
a greater than expected loss of reads, it is a good idea to go back to
that step and troubleshoot why reads are dropping out. The dada2
tutorial has more details about what can be changed at each step.

``` r
getN <- function(x) sum(getUniques(x)) # function to grab sequence counts from output objects
# tracking reads by counts
filt_out_track <- filt_out %>%
  data.frame() %>%
  mutate(Sample = gsub("(R1\\_)(.{1,})(\\.fastq\\.gz)","\\2",rownames(.))) %>%
  rename(input = reads.in, filtered = reads.out)
rownames(filt_out_track) <- filt_out_track$Sample

# If you did not pool samples, use this
ddF_track <- data.frame(denoisedF = sapply(ddF[sample.names], getN)) %>%
  mutate(Sample = row.names(.))
ddR_track <- data.frame(denoisedR = sapply(ddR[sample.names], getN)) %>%
  mutate(Sample = row.names(.))

# If you used samples = pooled, use this
# ddF_track <- data.frame(denoisedF = sapply(derepF.p[sample.names], getN)) %>% mutate(Sample = row.names(.)) 
# ddR_track <- data.frame(denoisedR = sapply(derepR.p[sample.names], getN)) %>% mutate(Sample = row.names(.))

merge_track <- data.frame(merged = sapply(mergers, getN)) %>%
  mutate(Sample = row.names(.))
chim_track <- data.frame(nonchim = rowSums(seqtab.nochim)) %>%
  mutate(Sample = row.names(.))

track <- left_join(filt_out_track, ddF_track, by = "Sample") %>%
  left_join(ddR_track, by = "Sample") %>%
  left_join(merge_track, by = "Sample") %>%
  left_join(chim_track, by = "Sample") %>%
  replace(., is.na(.), 0) %>%
  select(Sample, everything())
row.names(track) <- track$Sample
track
##                        Sample input filtered denoisedF denoisedR merged nonchim
## ANT7                     ANT7 40165    38746     36988     37217  30862   30726
## ANT8                     ANT8 25526    24723     23307     23554  18873   18762
## BA1A5566_10_1D BA1A5566_10_1D 50515    48763     48233     48059  46759   45399
## BA1A5566_22_1C BA1A5566_22_1C 42238    40623     40162     40225  38940   37117
## BA1A5566_45_1J BA1A5566_45_1J 31057    29830     29398     29519  28618   25794
## BB1S                     BB1S 34353    33063     31892     32190  28326   28191
## BB1W                     BB1W 39824    38337     38014     38028  33867   33359
## BNS1                     BNS1  3792     3684      3668      3667   3655    3655
## BNS2                     BNS2 16915    16362     16308     16295  16170   16170
## C2S                       C2S 38469    37037     35623     35976  31203   31087
## C2W                       C2W 41979    40359     40093     40185  39706   39477
## CM2A_10_9C         CM2A_10_9C 34309    33206     33116     33124  32979   32751
## CM2A_22_9J         CM2A_22_9J 31079    29844     29785     29786  29418   29418
## COL11                   COL11 31560    30423     27972     28588  20030   19875
## COL12                   COL12 24413    23608     21511     22048  15347   15055
## COL3                     COL3 40397    38933     35689     36566  25853   25666
## COL4                     COL4 32665    31406     28888     29523  21331   21070
## DT3S                     DT3S 36902    35316     34167     34553  30666   30663
## DT3W                     DT3W 44790    42982     42575     42697  41470   40453
## OM18_BC               OM18_BC  4803     4603      4564      4575   4526    4526
## OM18_BJ               OM18_BJ   452      391       377       382    376     376
## WAB105_22_6J     WAB105_22_6J 44993    43423     42958     43102  41792   41121
## WAB105_45_6D     WAB105_45_6D 37348    36026     35477     35596  34455   33207
## WAB188_10_4D     WAB188_10_4D 31134    29793     29587     29656  29107   29101
## WAB71_45_3D       WAB71_45_3D 43340    41726     41478     41484  40637   40607

# tracking reads by percentage
track_pct <- track %>% 
  data.frame() %>%
  mutate(Sample = rownames(.),
         filtered_pct = ifelse(filtered == 0, 0, 100 * (filtered/input)),
         denoisedF_pct = ifelse(denoisedF == 0, 0, 100 * (denoisedF/filtered)),
         denoisedR_pct = ifelse(denoisedR == 0, 0, 100 * (denoisedR/filtered)),
         merged_pct = ifelse(merged == 0, 0, 100 * merged/((denoisedF + denoisedR)/2)),
         nonchim_pct = ifelse(nonchim == 0, 0, 100 * (nonchim/merged)),
         total_pct = ifelse(nonchim == 0, 0, 100 * nonchim/input)) %>%
  select(Sample, ends_with("_pct"))

# summary stats of tracked reads averaged across samples
track_pct_avg <- track_pct %>% summarize_at(vars(ends_with("_pct")), 
                                            list(avg = mean))
track_pct_med <- track_pct %>% summarize_at(vars(ends_with("_pct")), 
                                            list(med = stats::median))
track_pct_avg
##   filtered_pct_avg denoisedF_pct_avg denoisedR_pct_avg merged_pct_avg
## 1         95.93076          97.20694          97.85535       91.05881
##   nonchim_pct_avg total_pct_avg
## 1        98.61956      84.18239
track_pct_med
##   filtered_pct_med denoisedF_pct_med denoisedR_pct_med merged_pct_med
## 1         96.27596          98.86517          98.95743       97.11918
##   nonchim_pct_med total_pct_med
## 1        99.41186      87.87585

# Plotting each sample's reads through the pipeline
track_plot <- track %>% 
  data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  gather(key = "Step", value = "Reads", -Sample) %>%
  mutate(Step = factor(Step, 
                       levels = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim"))) %>%
  ggplot(aes(x = Step, y = Reads)) +
  geom_line(aes(group = Sample), alpha = 0.2) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0)) + 
  stat_summary(fun = median, geom = "line", group = 1, color = "steelblue", size = 1, alpha = 0.5) +
  stat_summary(fun = median, geom = "point", group = 1, color = "steelblue", size = 2, alpha = 0.5) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5), 
               geom = "ribbon", group = 1, fill = "steelblue", alpha = 0.2) +
  geom_label(data = t(track_pct_avg[1:5]) %>% data.frame() %>% 
               rename(Percent = 1) %>%
               mutate(Step = c("filtered", "denoisedF", "denoisedR", "merged", "nonchim"),
                      Percent = paste(round(Percent, 2), "%")),
             aes(label = Percent), y = 1.1 * max(track[,2])) +
  geom_label(data = track_pct_avg[6] %>% data.frame() %>%
               rename(total = 1),
             aes(label = paste("Total\nRemaining:\n", round(track_pct_avg[1,6], 2), "%")), 
             y = mean(track[,6]), x = 6.5) +
  expand_limits(y = 1.1 * max(track[,2]), x = 7) +
  theme_classic()

track_plot
```

<img src="figure/unnamed-chunk-25-1.png" width="98%" height="98%" />

``` r
# Write results to disk
saveRDS(track, paste0(project.fp, "/tracking_reads.rds"))
saveRDS(track_pct, paste0(project.fp, "/tracking_reads_percentage.rds"))
saveRDS(track_plot, paste0(project.fp, "/tracking_reads_summary_plot.rds"))
```

## Next Steps

You can now transfer over the output files onto your local computer. The
table and taxonomy can be read into R with ‘mctoolsr’ package or another
R package of your choosing.

### Post-pipeline considerations

After following this pipeline, you will need to think about the
following in downstream applications:

1.  Remove ASVs assigned as mitochondrial or chloroplast
2.  Remove ASVs assigned as eukaryotes
3.  Remove ASVs that are unassigned at the domain level
4.  Remove ASVs that only have 1 or 2 sequences in the entire dataset
5.  Normalize or rarefy your ASV table

For more information on these steps, see
<https://github.com/cliffbueno/Alpha-Beta-Tutorial>
