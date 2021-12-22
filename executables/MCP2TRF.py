# -*- coding: utf-8 -*-
"""
 Author:
    Diana Y. Lee, Luque Lab, SDSU
    dlee@sdsu.edu
 Purpose:
    Base random forest model for predicting the capsid architecture (as measured by the 
    T-number) of a tailed phage from the MCP sequence
 Requires: 
    phageFunctions.py  :  Functions for calculating T based on genome size
    MCP2T_RF_state.db  :  Trained random forest model database
Input: 
    
    Two options: 
      1) a fasta format file for up to 2000 MCPs (protein format)
      2) an MCP phage data to analyze. Must include the following columns:
        'Virus_ID'
        'MCP_Sequence'
        'MCP_len'
        'IPC'
Output:
    MCP2TResults.csv  :  results of the random forest prediction
"""

# basic imports
import pandas as pd
import numpy as np
np.random.seed(42)
import os
import copy

# biopython imports
from Bio import SeqIO
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP

# custom function to count amino acids
# amino acids are hardcoded to avoid broken dependencies, since they do not change
def createFreq(acidSeq, normF=None):
    normF = normF or 0
    if (normF > 1):
        print("Valid tTypes are 0 (raw) or 1 (normalized by sample). Defaults to 0.")
        return
    AA = []
    aaList = np.asarray(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'])
    aaLen=len(aaList)
    n = len(acidSeq)
    for i in range(n):
        for acid in (aaList):
            trans1=str(acidSeq[i])
            a = trans1.count(str(acid))
            AA.append(a)
    rFreq = np.asarray(AA).reshape((n, aaLen))
    if (normF == 0):
#        print("Success! Created an nx20 array, where n is the length of the list provided:",n)
#        print("Columns are frequency totals for each amino acid:",aaList)
        return rFreq
    if (normF == 1):
        nFreq = copy.copy(rFreq).astype(float)
        fff3 = copy.copy(rFreq).astype(float)
        nf = rFreq.shape[1]
        for i in range(nf):
            nFreq[:,i] = fff3[:,i]/fff3.sum(axis=1)
#        print("Success! Created an nx20 array, where n is the length of the list provided:",n)
#        print("Columns are frequency percentages for each amino acid:",aaList)
        return nFreq

## load kernel state
## rerun imports cell as well
import dill
dill.load_session('MCP2T_RF_state.db')

# set the error margin
errMar = 0.09

MCP_Type = input("Get in the VAN! Enter MCP input type (1: fasta file, 2: csv file): ")
if MCP_Type=="1":
    print("Vroom! Let's go!")
    MCP_Fasta_Loc = input("Enter file location: ")
    
    assert os.path.exists(MCP_Fasta_Loc), "Error: file does not exist at "+str(MCP_Fasta_Loc)
    MCP_fasta_file = open(MCP_Fasta_Loc,'r+')
    print("File found. Processing.")
    # import a single fasta for a whole genome
    MCP_fasta_data = []
    for record in SeqIO.parse(MCP_fasta_file, "fasta"):
        #print(record.description)
        #print("genome length ",len(record.seq))
        MCP_fasta_data.append(str(record.description[0:75]))
        MCP_fasta_data.append(str(record.seq))
        MCP_fasta_data.append(len(record.seq))
        MCP_fasta_data.append(IP(str(record.seq)).pi())       
    MCP_fasta_file.close()

    MCP_fasta_data = np.reshape(np.ravel(MCP_fasta_data), (-1, 4))
    MCP_fasta_data = pd.DataFrame(MCP_fasta_data)
    MCP_fasta_data = MCP_fasta_data.rename(columns={0: 'Virus_ID', 1: 'MCP_Sequence', 2: 'MCP_len', 3: 'IPC'})
    MCP_fasta_data["MCP_len"] = MCP_fasta_data["MCP_len"].astype('int64')
    MCP_fasta_data["IPC"] = MCP_fasta_data["IPC"].astype('float64')
    
    n = len(MCP_fasta_data["Virus_ID"])

    # create the dataset for the random forest
    freq = createFreq(MCP_fasta_data["MCP_Sequence"], 1)
    AAT = []
    for i in range(n):
        AAT.append(MCP_fasta_data.iloc[i]["Virus_ID"])
        AAT.append(MCP_fasta_data.iloc[i]["IPC"])
        AAT.append(MCP_fasta_data.iloc[i]["MCP_len"])
        for j in range(20):
            AAT.append(freq[i][j])
    AAT = np.asarray(np.reshape(np.ravel(AAT), (n, 23)))
    x_FaaPhage = pd.DataFrame(AAT)

    # assign the features
    x_FaaActual = (x_FaaPhage.iloc[:,1:23]).astype(float)
    
    # predict T-numbers (ignore the error; the trained forest is imported in line 104)
    y_Pred = rfBest_clf.predict(x_FaaActual)  

    # create an output file
    y_PredTemp = []
    
    for i in range(n):
        y_PredTemp.append(x_FaaPhage[i][0])
        y_PredTemp.append(tdict2rev[y_Pred[i]])
    y_PredTemp = pd.DataFrame(np.reshape(np.ravel(y_PredTemp),(n,2)))
    y_PredTemp = y_PredTemp.rename(columns={0: 'Virus_ID', 1: 'T_pred'})
    phageMCP2TResult = MCP_fasta_data.merge(y_PredTemp, how='left', on='Virus_ID')
    phageMCP2TResult.to_csv(r'MCP2TResults.csv', index=False)
    print("Success! see MCP2TResults.csv")
    
elif MCP_Type=="2":
    print("Vroom! Let's go!")
    print("Your .csv file will require four columns: Virus ID, MCP_Sequence, MCP_len, and IPC")
    MCP_File_Loc = input("Enter file location: ")
    
    assert os.path.exists(MCP_File_Loc), "Error: file does not exist at "+str(MCP_File_Loc)
    MCP_data_file = open(MCP_File_Loc,'r+')
    MCPData = pd.read_csv(MCP_data_file)
    MCP_data_file.close()

    # count number of records
    n=MCPData.shape[0]

    # create the whole dataset including Virus ID
    freq = createFreq(MCPData["MCP_Sequence"], 1)
    AAT = []
    for i in range(n):
        AAT.append(MCPData.iloc[i]["Virus_ID"])
        AAT.append(MCPData.iloc[i]["IPC"])
        AAT.append(MCPData.iloc[i]["MCP_len"])
        for j in range(20):
            AAT.append(freq[i][j])
    AAT = np.reshape(np.ravel(AAT), (n, 23))
    x_Phage = np.asarray(AAT)
    # grab the subset of the dataset with just the features
    x_actual = (x_Phage[0:n,1:23]).astype(float)

    # predict T-numbers (ignore the error; the trained forest is imported in line 104)
    y_Pred = rfBest_clf.predict(x_actual)  

    # create an output file
    y_PredTemp = []
    for i in range(n):
        y_PredTemp.append(x_Phage[i,0])
        y_PredTemp.append(tdict2rev[y_Pred[i]])
    y_PredTemp = pd.DataFrame(np.reshape(np.ravel(y_PredTemp),(n,2)))
    y_PredTemp = y_PredTemp.rename(columns={0: 'Virus_ID', 1: 'T_pred'})
    phageMCP2TResult = MCPData.merge(y_PredTemp, how='left', on='Virus_ID')
    phageMCP2TResult.to_csv(r'MCP2TResults.csv', index=False)
    print("Success! see MCP2TResults.csv")
    
else:
    print("Sad honk. Invalid input.")
