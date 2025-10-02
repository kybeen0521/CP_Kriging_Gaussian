#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pump Efficiency Analysis with B-spline + Gaussian Kriging using Covariance
(실제 데이터 사용, Cubic Spline 직접 구현)

Author: Yongbeen Kim
Date: 2025-09-30
"""

import logging
from typing import Tuple
import matplotlib.pyplot as plt
import numpy as np

# -----------------------------
# Logging Configuration
# -----------------------------
logging.basicConfig(level=logging.INFO, format="%(message)s")
logger = logging.getLogger(__name__)

# -----------------------------
# Constants
# -----------------------------
RHO: float = 1000.0
G: float = 9.81
KRIGING_NOISE: float = 1e-12

# -----------------------------
# Input Data (동일)
# -----------------------------
RPM = np.full(20, 900)
P_IN = np.array([1.262, 1.262, 1.212, 0.858, 0.454, 0.0, -0.303, -0.555, -0.909, -1.262,
                 -1.464, -1.767, -1.969, -2.020, -2.322, -2.423, -2.474, -2.474, -2.575, -2.575]) * 1000.0
P_OUT = np.array([21.48, 20.78, 19.64, 18.15, 17.17, 15.45, 14.54, 13.91, 12.77, 11.86,
                  11.16, 10.25, 10.02, 9.74, 9.16, 9.04, 9.24, 9.14, 9.06, 9.06]) * 1000.0
Q = np.array([0.0527, 0.1191, 0.2793, 0.4258, 0.5449, 0.6641, 0.7168, 0.7695, 0.8242, 0.9023,
              0.9160, 0.9570, 0.9824, 1.0098, 1.0352, 1.0762, 1.0625, 1.0625, 1.0762, 1.0625]) / 1000.0
V_IN = np.array([0.1216, 0.2747, 0.6439, 0.9817, 1.2563, 1.5310, 1.6526, 1.7742, 1.9003, 2.0804,
                 2.1119, 2.2065, 2.2650, 2.3281, 2.3866, 2.4812, 2.4496, 2.4496, 2.4812, 2.4496])
V_OUT = np.array([0.2192, 0.4953, 1.1612, 1.7702, 2.2655, 2.7609, 2.9801, 3.1993, 3.4267, 3.7515,
                  3.8084, 3.9789, 4.0844, 4.1981, 4.3037, 4.4742, 4.4174, 4.4174, 4.4742, 4.4174])
HE = np.full(20, 0.075)
TORQUE = np.array([0.0402, 0.1098, 0.1345, 0.1484, 0.1561, 0.2041, 0.2041, 0.2242, 0.1994, 0.2535,
                   0.2473, 0.2597, 0.2674, 0.2891, 0.2736, 0.2922, 0.3061, 0.2953, 0.3138, 0.3308])

# -----------------------------
# Functions
# -----------------------------
def compute_efficiency(rpm, torque, p_in, p_out, q, v_in, v_out, he):
    omega = rpm * np.pi / 30.0
    w_shaft = torque * omega
    ha = (p_out - p_in) / (RHO * G) + he + (v_out**2 - v_in**2) / (2 * G)
    w_hydraulic = RHO * G * q * ha
    eta = np.clip((w_hydraulic / w_shaft) * 100.0, 0.0, 100.0)
    bep_idx = np.argmax(eta)
    return eta, q[bep_idx], eta[bep_idx]


def cubic_spline_interpolation(x, y, x_new):
    n = len(x)
    h = np.diff(x)
    alpha = np.zeros(n)
    for i in range(1, n-1):
        alpha[i] = (3/h[i])*(y[i+1]-y[i]) - (3/h[i-1])*(y[i]-y[i-1])
    l = np.ones(n)
    mu = np.zeros(n)
    z = np.zeros(n)
    c = np.zeros(n)
    b = np.zeros(n-1)
    d = np.zeros(n-1)
    for i in range(1, n-1):
        l[i] = 2*(x[i+1]-x[i-1]) - h[i-1]*mu[i-1]
        mu[i] = h[i]/l[i]
        z[i] = (alpha[i] - h[i-1]*z[i-1]) / l[i]
    for j in range(n-2, -1, -1):
        c[j] = z[j] - mu[j]*c[j+1]
        b[j] = ((y[j+1]-y[j])/h[j]) - h[j]*(c[j+1]+2*c[j])/3
        d[j] = (c[j+1]-c[j]) / (3*h[j])
    y_new = np.zeros_like(x_new)
    for idx, xi in enumerate(x_new):
        j = np.searchsorted(x, xi) - 1
        j = np.clip(j, 0, n-2)
        dx = xi - x[j]
        y_new[idx] = y[j] + b[j] * dx + c[j] * dx ** 2 + d[j] * dx ** 3
    return y_new


def gaussian_covariance(h, sill, vrange):
    """Gaussian covariance function"""
    return sill * np.exp(-(h / vrange) ** 2)


def ordinary_kriging_covariance_1d(x_data, y_data, x_pred, sill, vrange):
    """1D ordinary Kriging using covariance directly"""
    n = len(x_data)
    # Covariance matrix
    cov_matrix = np.array([[gaussian_covariance(abs(xi - xj), sill, vrange) for xj in x_data] for xi in x_data])
    cov_matrix += KRIGING_NOISE * np.eye(n)

    # Augmented system for Lagrange multiplier
    cov_aug = np.zeros((n + 1, n + 1))
    cov_aug[:n, :n] = cov_matrix
    cov_aug[-1, :-1] = 1.0
    cov_aug[:-1, -1] = 1.0
    cov_aug[-1, -1] = 0.0

    y_pred = np.zeros(len(x_pred))
    for i, xp in enumerate(x_pred):
        k = np.array([gaussian_covariance(abs(xp - xi), sill, vrange) for xi in x_data])
        sol = np.linalg.solve(cov_aug, np.append(k, 1.0))
        weights = sol[:-1]
        y_pred[i] = np.dot(weights, y_data)
    return y_pred


# -----------------------------
# Main Execution
# -----------------------------
if __name__ == "__main__":
    eta, bep_q, bep_eta = compute_efficiency(RPM, TORQUE, P_IN, P_OUT, Q, V_IN, V_OUT, HE)

    q_unique, unique_idx = np.unique(Q, return_index=True)
    eta_unique = eta[unique_idx]
    q_dense = np.linspace(q_unique.min(), q_unique.max(), 50)
    eta_spline_dense = cubic_spline_interpolation(q_unique, eta_unique, q_dense)

    # Covariance parameters (sill and range)
    sill_fit = np.var(eta_unique)
    vrange_fit = (q_unique.max() - q_unique.min()) / 2.0

    logger.info("=== Covariance Kriging Parameters ===")
    logger.info(f"sill={sill_fit:.4f}, range={vrange_fit:.4f}")

    # Kriging + plotting
    points_list = [10, 15, 25, 30]
    colors = ["purple", "blue", "green", "orange"]

    plt.figure(figsize=(12, 6))
    plt.plot(Q, eta, "ro-", label="Original Data")
    plt.scatter(bep_q, bep_eta, color="red", s=100, zorder=5)
    plt.text(bep_q, bep_eta + 1, f"BEP\n(Q={bep_q:.4f}, η={bep_eta:.2f}%)",
             color="red", fontweight="bold", ha="center")

    for n_pts, color in zip(points_list, colors):
        q_spline = np.linspace(q_unique.min(), q_unique.max(), n_pts)
        eta_spline = cubic_spline_interpolation(q_unique, eta_unique, q_spline)
        eta_pred = ordinary_kriging_covariance_1d(q_spline, eta_spline, q_dense, sill_fit, vrange_fit)
        plt.plot(q_dense, eta_pred, color=color, lw=2, label=f"B-spline + Covariance Kriging ({n_pts} pts)")

    plt.xlabel("Flow Rate Q [m³/s]")
    plt.ylabel("Efficiency η [%]")
    plt.title("Pump Efficiency: Cubic Spline + Gaussian Kriging (Covariance)")
    plt.ylim(0, 100)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
