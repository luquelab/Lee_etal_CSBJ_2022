""" Author:
   Diana Y. Lee, Luque Lab, SDSU
   dlee@sdsu.edu
   
 Purpose:
    This the current (as of 2021-11-10) G2T model for predicting the capsid architecture (as measured by the 
    T-number) of a tailed phage from the genome
 Requires: 
    phageFunctions.py  :  Functions for calculating T based on genome size
 Input: 
    This code requires genome length that can be entered in one of three ways:
    1) a single genome length
    2) a fasta file containing up to 5000 items
    3) a csv file containing at least two columns: "Virus_ID" and "genome_length" 
      
 Output:
    
    Single genome option results will display on screen. 
    All others will be provided in the file G2TResults.csv    """



# imports
import pandas as pd
import numpy as np
np.random.seed(42)
import os

from Bio import SeqIO

from phageFunctions import tNum

# set the error margin
errMar = 0.09

genome_Type = input("Get in the VAN! Enter input type (1: genome(bp), 2: fasta file, 3: csv file): ")
if genome_Type=="1":
    test_genome = float(input("Enter genome length in bp: "))

    T_raw_test = round(tNum(test_genome/1000,0),4)
    T_nearest_test = tNum(test_genome/1000,1)
    T_errMar_test = tNum(test_genome/1000,2,errMar)
    
    print("Success! ")
    print("T raw from model: ", T_raw_test)
    print("Nearest T, no restriction: ", T_nearest_test)
    print("Nearest T within specified error margin (",errMar,"): ", T_errMar_test)

elif genome_Type=="2":
    print("Vroom! Let's go!")
    genome_Fasta_Loc = input("Enter file location: ")
    
    assert os.path.exists(genome_Fasta_Loc), "Error: file does not exist at "+str(genome_Fasta_Loc)
    genome_fasta_file = open(genome_Fasta_Loc,'r+')
    print("File found. Processing.")
    # import a single fasta for a whole genome
    genome_fasta_data = []
    for record in SeqIO.parse(genome_fasta_file, "fasta"):
        #print(record.description)
        #print("genome length ",len(record.seq))
        genome_fasta_data.append(str(record.description[0:75]))
        genome_fasta_data.append(len(record.seq))
    genome_fasta_file.close()

    genome_fasta_data = np.reshape(np.ravel(genome_fasta_data), (-1, 2))
    genome_fasta_data = pd.DataFrame(genome_fasta_data)
    genome_fasta_data = genome_fasta_data.rename(columns={0: 'Virus_ID', 1: 'genome_length'})
    genome_fasta_data["genome_length"] = genome_fasta_data["genome_length"].astype('int64')
    
    n = len(genome_fasta_data["Virus_ID"])
    # calculate T numbers
    ny = genome_fasta_data.shape[0]
    Y_T = []
    
    for i in range(ny):
        Y_T.append(genome_fasta_data.iloc[i]["Virus_ID"])
        Y_T.append(round(tNum(genome_fasta_data.iloc[i]["genome_length"]/1000,0),4))
        Y_T.append(tNum(genome_fasta_data.iloc[i]["genome_length"]/1000,1))
        Y_T.append(tNum(genome_fasta_data.iloc[i]["genome_length"]/1000,2,errMar))
        
    Y = np.asarray(Y_T)
    Y = np.reshape(np.ravel(Y), (ny, 4));
    Y = np.asarray(Y)
    
    df_T = pd.DataFrame(Y)
    df_T = df_T.rename(columns={0: 'Virus_ID', 1: 'T_raw', 2: 'T_nearest', 3: 'T_nearest_errMar'})
    
    df_T["T_raw"] = df_T["T_raw"].astype('float64')
    df_T["T_nearest"] = df_T["T_nearest"].astype('float64')
    df_T["T_nearest_errMar"] = df_T["T_nearest_errMar"].astype('float64')
    
    # add T predictions to the phage data
    phageG2TResult = genome_fasta_data.merge(df_T, how='left', on='Virus_ID')
    phageG2TResult.to_csv(r'G2TResults.csv', index=False)
    print("Success! see G2TResults.csv")

elif genome_Type=="3":
    print("I can do that. Your .csv file will require two columns: Virus_ID and genome_length")
    genome_data_Loc = input("Enter file location: ")
    
    assert os.path.exists(genome_data_Loc), "Error: file does not exist at "+str(genome_Fasta_Loc)
    genome_data_file = open(genome_data_Loc,'r+')
    print("File found. Processing.")
    phageData = pd.read_csv(genome_data_file)
    genome_data_file.close()
    #print(phageData.dtypes)
    
    n = len(phageData["Virus_ID"])
    # calculate T numbers
    ny = phageData.shape[0]
    Y_T = []
    
    for i in range(ny):
        Y_T.append(phageData.iloc[i]["Virus_ID"])
        Y_T.append(round(tNum(phageData.iloc[i]["genome_length"]/1000,0),4))
        Y_T.append(tNum(phageData.iloc[i]["genome_length"]/1000,1))
        Y_T.append(tNum(phageData.iloc[i]["genome_length"]/1000,2,errMar))
        
    Y = np.asarray(Y_T)
    Y = np.reshape(np.ravel(Y), (ny, 4));
    Y = np.asarray(Y)
    
    df_T = pd.DataFrame(Y)
    df_T = df_T.rename(columns={0: 'Virus_ID', 1: 'T_raw', 2: 'T_nearest', 3: 'T_nearest_errMar'})
    
    df_T["T_raw"] = df_T["T_raw"].astype('float64')
    df_T["T_nearest"] = df_T["T_nearest"].astype('float64')
    df_T["T_nearest_errMar"] = df_T["T_nearest_errMar"].astype('float64')
    # add T predictions to the phage data
    phageG2TResult = phageData.merge(df_T, how='left', on='Virus_ID')    
    phageG2TResult.to_csv(r'G2TResults.csv', index=False)    #stuff you do with the file goes here
    print("Success! see G2TResults.csv")

else:
    print("Sad honk. Invalid option.")