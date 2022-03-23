# An efficient 146-line 3D sensitivity analysis code of stress-based topology optimization written in MATLAB

Please run the MATLAB program `run_stress_min.m` for the first numerical example in the paper.

**Recommended environment:** MATLAB 2020a

## Introduction

The repository here is for the paper: [An efficient 146-line 3D sensitivity analysis code of stress-based topology optimization written in MATLAB](https://link.springer.com/article/10.1007/s11081-021-09675-3). The 146 lines code includes the finite element analysis and p-norm stress sensitivity based on the adjoint method, verified by finite difference approximation. The code uses the [Method of Moving Asymptotes (MMA) Optimizer](http://www.smoptit.se/) as the nonlinear optimizer. The code can be extended for different stress related 3D topology optimization problems and is intended for educational purposes.

## Included Files

- MMA Optimization (These files can be downloaded from [here](http://www.smoptit.se/). Please email them if you are using the code)
  - `asymp.m`
  - `gcmmasub.m`
  - `kktcheck.m`
  - `subsolv.m`
- `prepare_filter.m`: Prepare the density filter to smooth the density variables and sensitivities.
- `run_stress_min.m`: Example to produce the first cantilever beam example from the paper.
- `Stress_3D_Sensitivity_Comp.m`: 146-line code to determine stress and sensitivities.
- `stress_minimize.m`: Computes filtered density variable for 146-line code and plots figures.

## Citation

If you have used this code, please cite our paper:
```
@Article{Deng2021,
author={Deng, Hao
and Vulimiri, Praveen S.
and To, Albert C.},
title={An efficient 146-line 3D sensitivity analysis code of stress-based topology optimization written in MATLAB},
journal={Optimization and Engineering},
year={2021},
month={Aug},
day={26},
issn={1573-2924},
doi={10.1007/s11081-021-09675-3},
url={https://doi.org/10.1007/s11081-021-09675-3}
}
```
