
# WHAT, WHO, WHEN

- Folder containing the executable files for T-number prediction based on Lee et al. Predicting the capsid architecture of tailed phages from metagenomic data https://doi.org/10.1016/j.csbj.2021.12.032

- Diana Y. Lee, Luque Lab, SDSU / dlee@sdsu.edu
- 2021-09-24, revised 2022-05-01



# FOLDERS & FILES

**FILE:** G2T.py\
--> This the current (as of 2021-11-10) G2T model for predicting the capsid architecture (as measured by the T-number) of a tailed phage from the genome\
--> Requires: \ 
    phageFunctions.py  :  Functions for calculating T based on genome size\
--> Input: \
    This code requires genome length that can be entered in one of three ways:\
    1) a single genome length\
    2) a fasta file containing up to 5000 items\
    3) a csv file containing at least two columns: "Virus_ID" and "genome_length"\ 
--> Output: \
    Single genome option results will display on screen.\ 
    All others will be provided in the file G2TResults.csv

## TO INSTALL:
--> save both the G2T.py and phageFunctions.py files locally.

## TO USE:
--> Run G2T.py via command line and follow the prompts. G2T.py requires input data in the form of:
    1) a single, integer genome length value in basepairs (bp)
    2) a fasta file with up to 1000 records
    3) a csv file that includes the columns "Virus_ID" and "genome_length". Any additional columns are left in place in the final result, but are not used.
    NOTE: We advise you to add a leading column before the data, ie: do not make "Virus_ID" or "genome_length" the first column. 
    We've encountered issues with text parsing in the first column header, and that bug is currently under investigation.



**FILE:** MCP2TRF.py\ 
--> Base random forest model for predicting the capsid architecture (as measured by the T-number) of a tailed phage from the MCP sequence \
--> Requires: \ 
    phageFunctions.py  :  Functions for calculating T based on genome size\
    MCP2T_RF_state.db  :  Trained random forest model database\
--> Input: \ 
    Two options: \ 
      1) a fasta format file for up to 2000 MCPs (protein format) \
      2) an MCP phage data .csv to analyze. Must include the following columns:\
        'Virus_ID'\
        'MCP_Sequence'\
        'MCP_len'\
        'IPC'\
--> Output: \
    MCP2TResults.csv  :  results of the random forest prediction

## TO INSTALL:
--> save the MCP2TRF.py, phageFunctions.py, MCP2T_RF_state.db files locally.

## TO USE:





**FILE:** MCP2TRF_state.db \
--> Trained random forest model database

**FILE:** phageFunctions.py \
--> Contains the basic functions that are reused across all jupyter notebooks associated with the G2T/MCP2T project
