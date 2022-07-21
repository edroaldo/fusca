# FUSCA
Framework for Unified Single-Cell Analysis (FUSCA) is a software package for single-cell data analysis developed in R. At the moment, it contains two modules: 1) CellRouter, for reconstruction of complex single-cell trajectories, and 2) CellComm, to infer intercellular communication networks from scRNA-seq data.

:telephone_receiver: Contact: edroaldo.lummertz@ufsc.br / edroaldo@gmail.com

## :notebook_with_decorative_cover: Tutorial for preprocessing and processing scRNA-Seq data using FUSCA (data: 3k PBMC) [Link](https://github.com/edroaldo/fusca/blob/main/tutorial/FUSCA_tutorial_3k_PBMC_ds.ipynb).


# Important note:
When running the CellComm tutorial, please make sure that your cell type names, or any other annotation used for inference of cell communication, do not contain underscores. For instance, intead of using T_cell or B_cell, please, use Tcell and Bcell.

## :notebook_with_decorative_cover: Tutorial to apply CellComm to spatial transcriptomics data [Link](https://github.com/edroaldo/fusca/blob/main/tutorial/CellComm_tutorial.ipynb).
