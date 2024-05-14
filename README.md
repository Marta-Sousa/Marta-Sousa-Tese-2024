# Marta-Sousa-Thesis_2024
 This is the deliverables for my master thesis.
 In this repository are 2 files in python. One to compute the acoustic parametrics using Parselmouth and the MFCC, wich is in "All_features.ipynb", while "Models.ipynb" is the code for the feature selection, and to train, test and validate the SVM, NN and XGBoost models.
The other file, the "Pearson_r.mat" is the code developed to compute the Pearson's correlation coefficient, also used as an input, that requires some data from "All_features.ipynb" such as the funcamental frequency of each siganl.
As for the pastes, _Database_ contains the data from VOiCeD (sub-pastes _waves1s_ and _waves2s_) and SVD, the paste _Dados_ has the products of the file "All_features.ipynb" ans "Pearson_r.mat", while _Modelos_ has the results of "Models.ipynb"
