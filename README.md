# ATE of Parental IQ: Unconditional and Conditional Bounds Estimation

This repository contains the code accompanying our paper, where we estimate both **unconditional** and **conditional** bounds using **Chernozhukov et al.'s half-median unbiased estimator**.

## Repository Structure

- **`code.ado`**: Prepares the data for analysis.
- **`omega.ado`**: Estimates the **unconditional bounds**.
- **`omega_conditional.ado`**: Estimates the **conditional bounds**.

## Requirements

This code is written in **Stata** and requires the following dependencies:
- Stata 15 or later  
- `moremata` package (for matrix operations)  
- `gtools` (for fast data processing)  

You can install these by running:  
```stata
ssc install moremata, replace
ssc install gtools, replace
