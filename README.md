# Pump Efficiency Analysis with B-spline + Gaussian Kriging
## 📌 Overview
This project analyzes pump efficiency data using:
Natural Cubic Spline Interpolation (direct implementation, no SciPy)
Gaussian Variogram Modeling
1D Ordinary Kriging with Adjustable Noise

The workflow integrates hydraulic efficiency computation with spatial statistics (kriging) to generate smooth, reliable efficiency curves from experimental pump test data.
---

## ✨ Features
- Direct computation of pump efficiency and Best Efficiency Point (BEP).
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
├── Input Data              # Pump test data (pressure, flow, torque, etc.)
├── compute_efficiency()    # Pump efficiency calculation
├── cubic_spline_interpolation() # Natural cubic spline
├── gaussian_variogram()    # Gaussian variogram model
├── ordinary_kriging_1d()   # 1D Ordinary Kriging solver
├── semivariogram()         # Experimental variogram computation
├── compute_rmse(), compute_r2() # Evaluation metrics
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

## 📊 Variogram Parameters

Nugget (γ(0)) → measurement noise / micro-scale variation.
Sill → variance plateau (long-range variability).
Range → effective distance of spatial correlation.

These parameters are automatically extracted from experimental variogram in the code.

🧩 Example Output

Red dots → Original pump efficiency test data.

Spline curve → Smooth reference (Cubic Spline).

Colored curves → Kriging predictions with different sampling resolutions (10, 15, 25, 30 pts).

BEP (Best Efficiency Point) marked with red circle and label.

Logs provide estimated variogram parameters and performance metrics for each kriging run:
