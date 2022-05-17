This repository contains the R code of the algorithm Magma, presented in the paper 'MAGMA: Inference and Prediction with Multi-Task Gaussian Processes' by Leroy, Latouche, Guedj and Gey (https://link.springer.com/article/10.1007/s10994-022-06172-1). If you are considering citing this work, you are welcome to use the file `reference.bib`.

The folder 'Simulations' contains the synthetic datasets, used for the experimental study in the paper, as well as the corresponding trained models and tables of results.

The file 'Simulation_study.R' contains the code used to generate synthetic datasets, conduct the experiments, evaluate and display the results. 

The file 'Computing_functions.R' contains the code of many useful functions used in MAGMA such as kernels, log-likelihoods, gradients, E-step and M-step.

The file 'MAGMA.R' contains the core of the code for the MAGMA algorithm. Fonctions implementing the learning, hyperposterior computing, and prediction steps are provided.
Several ploting function and testing examples are provided as well. 

An extended and up-to-date R package, called MagmaClustR, implementing Magma with additional features can be found at https://github.com/ArthurLeroy/MagmaClustR and on the CRAN.

Corresponding author and main contributor to the code:  
Arthur Leroy  
arthur.leroy.pro@gmail.com