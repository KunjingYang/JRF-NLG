
# A Coarse-to-Fine Hybrid Registration & Fusion framework for Hyperspectral Image Super-Resolution via Batch Image Alignment

---

## 📝 Overview
This package provides the **MATLAB implementation** of our paper:

**A Coarse-to-Fine Hybrid Registration and Fusion Framework for Hyperspectral Superresolution via Batch Image Alignment**  
*SIAM Journal on Imaging Sciences, 2026*

### Key Features
- **Coarse Stage**: Joint Registration and Fusion (JRF) model with batch image alignment
- **Fine Stage**: Nonconvex Low-rank and Group-sparse (NLG) fusion model
- Simultaneously addresses **misalignment**, **spectral variability**, and **detail loss**
- Output: High-spatial-resolution hyperspectral image (HR-HSI)

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