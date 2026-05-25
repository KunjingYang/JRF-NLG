
1. Overview
This package contains the MATLAB implementation of the paper:
``A Coarse-to-Fine Hybrid Registration and Fusion Framework for Hyperspectral Superresolution via Batch Image Alignment''
SIAM Journal on Imaging Sciences, 2026

Coarse stage: Joint Registration and Fusion (JRF) model with batch image alignment
Fine stage: Nonconvex Low-rank and Group-sparse (NLG) fusion model
Simultaneously handles misalignment, spectral variability, and detail loss
Output: High-spatial-resolution hyperspectral image (HR-HSI)


2. File Descriptions
├── demo_cave : Main script for Cave  dataset
├── demo_pavia: Main script for Pavia dataset
├── JRF_outer : Registration framework
├── JRF_inner : Generalized Gauss-Newton algorithm and ADMM
├── NLRGS_fus : the algorithm for solving the NLG model 
├── NLRGS_cave: Main script for Cave  dataset for Registed HSI and MSI
└── NLRGS_Pavia: Main script for Pavia dataset for Registed HSI and MSI


3. If you find this code helpful, please consider citing our work.
@article{yang2026coarse,
  title={A Coarse-to-Fine Hybrid Registration and Fusion Framework for Hyperspectral Superresolution via Batch Image Alignment},
  author={Yang, Kunjing and Bai, Minru and Lu, Ting and Chen, Liang},
  journal={SIAM Journal on Imaging Sciences},
  volume={19},
  number={2},
  pages={1207--1243},
  year={2026},
  publisher={SIAM}
}
