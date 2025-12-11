# Fisher Randomization Tests with Conformal Selective Borrowing

This repository contains replication code for the ICML 2025 paper:

Zhu, K., Yang, S.\*, & Wang, X. (2025). [Enhancing statistical validity and power in hybrid controlled trials: A randomization inference approach with conformal selective borrowing.](https://arxiv.org/abs/2410.11713) *Proceedings of the 42nd International Conference on Machine Learning (ICML)*, PMLR, in press. [[Slides]](https://drive.google.com/file/d/1LkTDY12CjUL0BAGQAt4JmOBlgk-MN3Ix/view?usp=sharing) [[Poster]](https://drive.google.com/file/d/1g5vFT6irtPWFQWwh6AGe-iYCvMF4z0B2/view?usp=share_link) [[Package]](https://github.com/ke-zhu/intFRT)


## Folder Structure

- `sim_main/`: Main simulation study 
- `sim_power_curve/`: Simulations for generating power curves
- `sim_corr_X/`: Additional simulations with correlated covariates, reported in the appendix
- `C9633+NCDB`: Code for real data analysis. Data set is not included due to data privacy restrictions

## Citation

If you use this code or refer to our methodology, please cite:

```bibtex
@inproceedings{zhu2025enhancing,
  title={Enhancing Statistical Validity and Power in Hybrid Controlled Trials: A Randomization Inference Approach with Conformal Selective Borrowing},
  author={Zhu, Ke and Yang, Shu and Wang, Xiaofei},
  booktitle={Proceedings of the 42nd International Conference on Machine Learning},
  year={2025},
  publisher={PMLR}
}
```
