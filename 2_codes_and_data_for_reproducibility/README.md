# WHAT, WHO, WHEN

- Folder containing the Jupyter notebook files associated with Lee et al. Predicting the capsid architecture of phages from metagenomic data

- Diana Y. Lee, Luque Lab, SDSU / dlee@sdsu.edu
- 2021-09-24, revised 2021-11-12


# FOLDERS & FILES

**FILE:** phage_functions.ipynb\
--> Jupyter notebook for Python 3.0\ 
Creates the basic functions that are reused across all jupyter notebooks associated with the G2T/MCP2T project

**FILE:** Figure_3_Statistical_Model_Updated.ipynb\
--> Jupyter notebook for Python 3.0 \
Creates a linear regression model for relating the genome length of a tailed phage to the capsid architecture \
as measured by the T-number.\
Evaluates the accuracy of that model.\
--> Requires:\
    data\HRphagedata.csv\
--> Creates (results folder):\
    Fig3_DNA_vs_T-number(ticks).svg\
    Fig3_DNA_vs_T-number.svg .png\
    Fig3_DNA_vs_T_model_Rarification.svg .png

**FILE:** Figure_4_model_application_to_MCP_Genome_DB_Updated.ipynb\
--> Jupyter notebook for Python 3.0 \
Creates a kernel distribution for analyzing the MCP database\
Applies the G2T model to the MCP database\
--> Requires:\
    phage_functions.ipynb\
    data\AllDNA.csv\
    data\PHAGE_TABLE3.xlsx\
--> Creates (results folder):\
    Fig4a_MCP_Kernel_Density.svg / .png\
    Fig4b_Predicted_T-numbers_for_MCP.svg / .png\
    Fig4c_Percent_architectures_predicted_by_T-number.svg / .png\
    Fig4_tWhole.csv\
    FigS2_MCP_Genome_length_histogram.svg / .png

**FILE:** Figure_5_MCP2T_Similarity_Model_Updated.ipynb
--> Jupyter notebook for Python 3.0 \
Creates an algorithm to assess T-number of a phage based on similarity of MCP\
Assesses the accuracy of that algorithm\
--> Requires:\
    phage_functions.ipynb \
    data\PHAGE_TABLE4.xlsx : phage data with indexes, genome size, and translations\
    data\MCPdbblastpe001.csv : Similarity table from blastp (e value = 0.001)\
--> Creates (results folder):\
    Fig5a_Tdist_vs_Similarity(e.001)_violin.png\
    Fig5c_Similarity_vs_accuracy_and_predictions.svg

**FILE:** Figure_6_Phage_Manuscript_Updated.ipynb\
--> Jupyter notebook for Python 3.0 \
Creates a random forest model for predicting the capsid architecture (as measured by the T-number) \
    of a tailed phage from the MCP sequence\
Evaluates the accuracy of that model.\
--> Requires:\
    phage_functions.ipynb\
    data\PHAGE_TABLE4.xlsx : phage data with indexes, genome size, and translations\
    data\clade1.csv : the list of phages identified based on the tree\
--> Creates (results folder):\
    Fig6b_RF_Predicted_T-numbers_for_MCP.png / .svg\
    Fig6c_RF_model_Rarification.png / .svg\
    FigS4_RF_confusion_Matrix.png / svg\
    FigS5a_RF_impurity-based_feature_importance.png / .svg\
    FigS5b_RF_Permutation_feature_importance.png / .svg\
    FigS5c_RF_Dropout_feature_importance.png / .svg

**FILE:** Figure_7_NCBI data_Updated.ipynb\
--> Jupyter notebook for Python 3.0 \
Processes PHANNs output for the NCBI data from Benler et al \
--> Requires:\
    phage_functions.ipynb\
    data\All_PHANNs_results.csv\
    data\NCBI_Genome_Len.csv\
    data\ac11.faa, data\ac13.faa, data\ac14.faa\
    data\MCP2T_RF_state.db\
--> Creates (results folder):\
    NCBI_all.csv\
    NCBI_predicted.csv

**FILE:** Figure_S6_Sim_v_RF_Timing.ipynb\
--> Jupyter notebook for Python 3.0 \
    Times the Similarity algorithms\
    Implements hypothetical large datasets\
    Compares the Similarity timing with (imported) Random Forest timing \
--> Requires:\
    phage_functions.ipynb \
    data\PHAGE_TABLE1.xlsx : phage data with indexes, genome size, and translations\
    data\MCPdbblastpe001.csv : Similarity table from blastp (e value = 0.001)\
--> Creates (results folder):\
    FigS6_sim_v_RF_time_train.svg \ .png\
    FigS6_sim_v_RF_time_predict.svg \ .png\
    FigS6_sim_v_RF_time_all.svg \ .png\
    FigS6_sim_v_RF_time_train_to_1M \ .png\
    FigS6_sim_v_RF_time_predict_to_1M.svg \ .png
