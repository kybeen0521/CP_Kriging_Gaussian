# Pump Efficiency Analysis with cubic spline + Gaussian Kriging
## 📌 Overview
This project analyzes pump efficiency data using:
Cubic Spline Interpolation (direct implementation)
Gaussian Variogram Model 
1D Ordinary Kriging with adjustable noise

The workflow integrates hydraulic efficiency computation with spatial statistics (kriging) to generate smooth, reliable efficiency curves from experimental pump test data.

---

## ✨ Features
- Direct computation of pump efficiency and Best Efficiency Point.
- Cubic Spline interpolation implemented manually.
- Experimental variogram estimation from data.
- Gaussian variogram fitting with nugget, sill, and range extraction.
- Ordinary Kriging (1D) implementation with noise stabilization.
- Model evaluation with RMSE and R² score.
- Visualization of raw data, spline approximation, and kriging prediction.

---

## 📂 Code Structure
```
pump_kriging.py

├── Constants               # Physical constants & noise parameter
│
├── Input Data              # Pump test data (pressure, flow, torque, etc.)
│
├── compute_efficiency()    # Pump efficiency calculation
│
├── cubic_spline_interpolation() # Natural cubic spline
│
├── gaussian_variogram()    # Gaussian variogram model
│
├── ordinary_kriging_1d()   # 1D Ordinary Kriging solver
│
├── semivariogram()         # Experimental variogram computation
│
├── compute_rmse(), compute_r^2() # Evaluation metrics
│
└── Main Execution          # Workflow: efficiency → variogram → kriging → plot
```
---

## 🔑 Kriging Workflow in This Project

1) Efficiency Data Preparation
- Computes pump efficiency (η) as function of flow rate (Q).
- Selects unique Q values for modeling.

2) Experimental Variogram
- Calculates semivariogram from efficiency differences.
- Extracts nugget, sill, and range parameters.

3) Gaussian Variogram Model
- Fits a parametric variogram using Gaussian kernel.

4) Ordinary Kriging (1D)
- Builds covariance matrix from variogram.
- Adds small white noise (KRIGING_NOISE) to diagonal for numerical stability.
- Solves augmented system (with Lagrange multiplier) for weights.
- Produces predictions and variance estimates.

5) Model Evaluation
- Computes RMSE and R² between spline ground truth and kriging predictions.


---


## 📊 Variogram Parameters

- Nugget (γ(0)) → measurement noise / micro-scale variation.
- Sill → variance plateau (long-range variability).
- Range → effective distance of spatial correlation.

Dependencies:

- numpy
- matplotlib


These parameters are automatically extracted from experimental variogram in the code.


---

## 📖 Key Insights

- Kriging is applied to efficiency vs flow rate curve, treating flow rate as spatial coordinate.
- White noise ensures stable inversion of covariance matrix.
- Variogram is the conceptual bridge: extracted from data → used in kriging (model covariance).
- Combined spline + kriging workflow allows smooth and statistically consistent efficiency prediction.


---


## 👤 Author
**Yongbeen Kim**  
Researcher, Intelligent Mechatronics Research Center, KETI
Email address: ybin521@keti.re.kr


📅 Document last updated 2025.10.09
