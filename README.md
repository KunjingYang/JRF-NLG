
# A Coarse-to-Fine Hybrid Registration & Fusion framework for Hyperspectral Image Super-Resolution



## 📝 Overview
This repository contains the official implementation of the paper:
> A Coarse-to-Fine Hybrid Registration and Fusion Framework for Hyperspectral Superresolution via Batch Image Alignment, Kunjing Yang, Minru Bai, Ting Lu, Liang Chen, SIAM Journal on Imaging Sciences, 19(2), 1207–1243, 2026

Paper: [https://epubs.siam.org/doi/full/10.1137/25M1758738](https://epubs.siam.org/doi/full/10.1137/25M1758738)  
arXiv: [https://arxiv.org/html/2407.05279v1](https://arxiv.org/html/2407.05279v1)

---

## 📂 File Descriptions
```matlab
├── demo_cave.m     % Main script for CAVE dataset experiments
├── demo_pavia.m    % Main script for Pavia dataset experiments
├── JRF_outer.m     % Outer registration framework
├── JRF_inner.m     % Generalized Gauss-Newton algorithm & ADMM solver
├── NLRGS_fus.m     % Core solver for NLG fusion model
├── NLRGS_cave.m    % Fusion script for registered HSI & MSI (CAVE)
└── NLRGS_Pavia.m   % Fusion script for registered HSI & MSI (Pavia)
```

## 📌 How to Run
1. Place your dataset in the data folder
2. Run the demo script corresponding to your dataset
```matlab
demo_cave    % for CAVE dataset
demo_pavia   % for Pavia dataset
```

## ✉️ Contact
Kunjing Yang
Email: kunjing-yang@hnu.edu.cn

## 📚 Citation
If you find this code helpful or use it in your research, please cite our paper:
```bibtex
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
```