
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


**FILE:** MCP2TRF_state.db \
--> Trained random forest model database

**FILE:** phageFunctions.py \
--> Contains the basic functions that are reused across all jupyter notebooks associated with the G2T/MCP2T project
