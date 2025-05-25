# Optimized Dimensionality Reduction for Moment-based Distributionally Robust Optimization

The software and data in this repository are a snapshot of the software and data that were used in the research reported in the paper Optimized Dimensionality Reduction for Moment-based Distributionally Robust Optimization by Shiyi Jiang, Jianqiang Cheng, Kai Pan, and Zuo-Jun Max Shen.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2024.0719

https://doi.org/10.1287/ijoc.2024.0719.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{Jiang2025,
  author =        {Jiang, Shiyi AND Cheng, Jianqiang AND Pan, Kai AND Shen, Zuo-Jun Max},
  publisher =     {Operations Research},
  title =         {Optimized Dimensionality Reduction for Moment-based Distributionally Robust Optimization},
  year =          {2025},
  doi =           {10.1287/ijoc.2024.0719.cd},
  url =           {https://github.com/jsy1164014200/ODR-MDRO},
  note =          {Available for download at https://github.com/jsy1164014200/ODR-MDRO},
}  
```

## Description

The goal of this repository is to share code and data related to the moment-based distributionally robust optimization (DRO) problem, solved using our proposed optimized dimensionality reduction (ODR) approach. 

## Building 

We propose an optimized dimensionality reduction (ODR) approach to solve the moment-based distributionally robust optimization (DRO) problem, which leads to two outer and one inner approximations of the original DRO problem. As these approximations are nonconvex  semidefinite programming (SDP) problems, we develop modified alternating direction method of multipliers (ADMM) algorithms to solve them efficiently. We demonstrate the effectiveness of our proposed ODR approach and algorithm in solving multiproduct newsvendor and production-transportation problems.

## Prerequisites

- MATLAB R2022a or higher
- CVX 
- Mosek 
- BMIBNB
 

## Structure

The source code is available in the [src](src) directory, which includes six folders.
- `Example1`: This folder includes the code used in Example 1. 
- `Example2`: This folder includes the code used in Example 2. 
- `MultiproductNewsvendor`: This folder includes the code used to solve multiproduct newsvendor problems.
- `ProductionTransportation`: This folder includes the code used to solve production transportation problems.
- `Insight`: This folder includes the code used in Section 7.2.3.
- `BMIBNB`: This folder includes the code used in Appendix F.1.

### Example1

This folder includes the code used in Example 1. 
- `CVaR.m`: Defines the procedure to solve Problem (8) and its PCA approximations in Example 1. 
- `CVaRPrimary.m`: Defines the `CVaRPrimary` function, which solves the SDP reformulation in Example 1. 

### Example2

This folder includes the code used in Example 2. 
- `SDP_fixX.m`: Defines the `SDP_fixX` function, which solves Problems (55) and (56) in Example 2. 

### MultiproductNewsvendor

This folder includes the code used to solve multiproduct newsvendor problems.
- `GenerateData.m`: Defines the procedure to generate data of multiproduct newsvendor problems.
- `ADMMnewsvendor.m`: Defines the procedure to use Algorithm 1 to solve two outer and one inner approximations under the ODR approach. 
- `PrimaryandPCA.m`: Defines the procedure to use the Mosek solver with default settings to solve the original DRO problem and its PCA approximations. 
- `covTransformerDecomposer.m`: Defines the `covTransformerDecomposer` function, which does eigenvalue decomposition for $\Sigma$.
- `ADMMnewsvendorGivenB.m`: Defines the `ADMMnewsvendorGivenB` function, which solves the first subproblem (i.e., Problem (33) given $\mathbf{B}$ and $\boldsymbol{\beta}$) of the ADMM algorithm of the first outer approximation. 
- `upperboundADMMGivenB.m`: Defines the `upperboundADMMGivenB` function, which solves the first subproblem of the ADMM algorithm of the inner approximation. 
- `lowerboundADMMGivenB.m`: Defines the `lowerboundADMMGivenB` function, which solves the first subproblem of the ADMM algorithm of the second outer approximation. 
- `PrimaryNewsvendor.m`: Defines the `PrimaryNewsvendor` function, which solves the SDP reformulation of the original problem.
- `BurerAlgorithm.m`: Defines the `BurerAlgorithm` function, which uses the low-rank algorithm to solve the SDP reformulation of the original problem. 
- `PrimaryNewsvendorGivenB.m`, and `PrimaryNewsvendorUpperGivenB.m`: Define two functions, which solve the PCA approximations of the original problem.


### ProductionTransportation 

This folder includes the code used to solve production transportation problems.
- `generate_data.m`: Defines the procedure to generate data of production transportation problems.
- `exe_admm.m`: Defines the procedure to use Algorithm 1 to solve our proposed outer and inner approximations under the ODR approach. 
- `exe_primary_dual_pca_vs.m`: Defines the procedure to use the Mosek solver with default settings to solve the original problem and its PCA approximations. 
- `covTransformerDecomposer.m`: Defines the `covTransformerDecomposer` function, which does eigenvalue decomposition for $\Sigma$.
- `admm_lowerbound.m`: Defines the `admm_lowerbound` function, which solves the first subproblem (i.e., Problem (33) given $\mathbf{B}$ and $\boldsymbol{\beta}$) of the ADMM algorithm of the first outer approximation. 
- `admm_upperbound.m`: Defines the `admm_lowerbound` function, which solves the first subproblem of the ADMM algorithm of the inner approximation. 
- `admm_revisted_lowerbound.m`: Defines the `admm_revisted_lowerbound` function, which solves the first subproblem of the ADMM algorithm of the second outer approximation. 
- `primary_sdp.m`: Defines the `primary_sdp` function, which solves the SDP reformulation of the original problem.
- `low_rank.m`: Defines the `low_rank` function, which uses the low-rank algorithm to solve the SDP reformulation of the original problem.
- `pca_lowerbound.m`, and `pca_upperbound.m`: Define two functions, which solve the PCA approximations of the original problem, which generate PCA-based lower bound and upper bound, respectively.

### Insight 

This folder includes the code used in Section 7.2.3.
- `Insight.m`: Defines the `Insight` function, which solves Problem (43). 

### BMIBNB

This folder includes the code used in Appendix F.1.
- `yalmipBiSDP.m`: Defines the `yalmipBiSDP` function, which uses the BMIBNB solver to solve the bilinear SDP problem. 


## Data 

The instance data is available in the [data](data) directory, which includes two folders.
- `MultiproductNewsvendor`: This folder includes all instances of multiproduct newsvendor problems. Each instance (denoted by "m_i") represents a multiproduct newsvendor problem with m products. Note that i is the instance id. 
- `ProductionTransportation`: This folder includes all instances of production transportation problems. Each instance (denoted by "K_G_H_i") represents a production transportation problem with G suppliers and H customers. Note that K is the number of pieces in the objective function and i is the instance id. 

Note that the data used in Appendix F.1 is also included in the folder `MultiproductNewsvendor`. The six instances in Table F2 are `100-1`, `100-2`, `100-3`, `200-1`, `200-2`, and `200-3`. 

Note that the data used in Example 1, Example 2, and Insight (Section 7.2.3) is included in the corresponding code. 


## Results

The results are available in the [results](results) directory, which includes three folders.
- `MultiproductNewsvendor`: This folder includes all results of multiproduct newsvendor problems. 
- `ProductionTransportation`: This folder includes all results of production transportation problems. 
- `BMIBNB`: This folder includes the results of using the BMIBNB solver to solve the six instances in Table F2. This folder includes six files. Each file records the log information given by the BMIBNB solver. For example, the file `100-1` records the log information of using the BMIBNB solver to solve the instance `data/MultiproductNewsvendor/100-1`.

### MultiproductNewsvendor

This folder includes two folders.
- `PCA`: This folder includes the results of three benchmark approaches. Each file (named "m") records the results of using three benchmark approaches to solve multiproduct newsvendor problems with m products. For example, the file "100" records the results of solving 5 instances with the prefix "100-" in the folder "data/MultiproductNewsvendor." Each file includes a table and each row records the results of an instance. The first and second columns are the objective value and computational time of the first benchmark approach (Mosek). The third and fourth columns are the objective value and computational time of the second benchmark approach (Low-rank Algorithm). The remaining columns are the objective values and computational times of PCA approximations with reduced dimension $m_1$ equals $100\% \times {\rm{dim}}(\boldsymbol{\xi})$, 
$80\%\times {\rm{dim}}(\boldsymbol{\xi})$, 
$60\%\times {\rm{dim}}(\boldsymbol{\xi})$, 
$40\%\times {\rm{dim}}(\boldsymbol{\xi})$, $20\%\times {\rm{dim}}(\boldsymbol{\xi})$, and $ K $, respectively.
- `ADMM`: This folder includes the results of our proposed ODR approach. Each file (named "m") records the results of using the ODR approach to solve multiproduct newsvendor problems with m products. For example, the file "100" records the results of solving 5 instances with the prefix "100-" in the folder "data/MultiproductNewsvendor." Each file includes a table and each row records the results of an instance. The first, second, and third columns are the lower bound, computational time, and theoretical gap of the ODR-LB. The fourth, fifth, and sixth columns are the upper bound, computational time, and theoretical gap of the ODR-UB. The remaining columns are the lower bound, computational time, and theoretical gap of the ODR-RLB.


### ProductionTransportation

This folder includes three folders: `K5`, `K10`, and `K15`. They include the results of production transportation problems when $K=5$, $K=10$, and $K=15$, respectively. Each folder includes two folders. The `PCA` folder includes the results of three benchmark approaches. The `ADMM` folder includes the results of our proposed ODR approach. 
- `K5/PCA`, `K10/PCA`, and `K15/PCA`: These folders include the results of three benchmark approaches. Each file (named "K-G-H") records the results of using three benchmark approaches to solve production transportation problems with G suppliers and H customers. For example, the file "5-4-25" records the results of solving 5 instances with the prefix "5-4-25-" in the folder "data/ProductionTransportation."
- `K5/ADMM/m1=3`, `K5/ADMM/m1=5`, `K5/ADMM/m1=7`: These folders include the results of our proposed ODR approach. Note that $m1=3$ means that we set $m_1=3$ in our proposed ODR approach although $K=5$. These results provide sensitivity analyses for our ODR approach. 
- `K10/ADMM/m1=8`, `K10/ADMM/m1=10`, `K10/ADMM/m1=12`: These folders provide results related to sensitivity analyses. 
- `K15/ADMM/m1=13`, `K15/ADMM/m1=15`, `K15/ADMM/m1=17`: These folders provide results related to sensitivity analyses. 
- `K10/ADMM/10-5-20-m1=2`, `K10/ADMM/10-5-20-m1=4`, ..., and `K10/ADMM/10-5-20-m1=40`: These files include results of sensitivity analyses presented in Figure 1. 



## Run

To execute the code, we just need to run the corresponding code file in MATLAB. For example, if we would like to use our ODR approach to solve the multiproduct newsvendor problem, we can 
```bash
cd src/MultiproductNewsvendor
./ADMMnewsvendor
```
