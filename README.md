#Pump Efficiency Analysis with B-spline + Gaussian Kriging
##📌 Overview

This project analyzes pump efficiency data using:

Natural Cubic Spline Interpolation (direct implementation, no SciPy)

Gaussian Variogram Modeling

1D Ordinary Kriging with Adjustable Noise

The workflow integrates hydraulic efficiency computation with spatial statistics (kriging) to generate smooth, reliable efficiency curves from experimental pump test data.

✨ Features

Direct computation of pump efficiency and Best Efficiency Point (BEP).

Cubic Spline interpolation implemented manually.

Experimental variogram estimation from data.

Gaussian variogram fitting with nugget, sill, and range extraction.

Ordinary Kriging (1D) implementation with noise stabilization.

Model evaluation with RMSE and R² score.

Visualization of raw data, spline approximation, and kriging prediction.

📂 Code Structure
