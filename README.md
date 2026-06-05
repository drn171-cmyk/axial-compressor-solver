# Axial Compressor Solver

1D mean-line and radial equilibrium aerothermal analysis, velocity triangles calculation, and CAD geometry generation for axial compressors in Python. 

This repository contains a purely analytical thermodynamic and kinematic solver designed for multi-stage axial compressors. It bridges the gap between initial theoretical thermodynamic cycle requirements and 3D CAD modeling for CFD analysis.

## 🚀 Features

- **Multi-Spool Capability:** Solves for Low-Pressure Compressor (LPC) and High-Pressure Compressor (HPC) sequentially.
- **Radial Equilibrium:** Extends mean-line analysis to 5 radial stations (Hub, 25%, Mean, 75%, Tip) using the **Free Vortex** design assumption ($r \cdot V_{\theta} = constant$).
- **Thermodynamic Tracking:** Calculates stage-by-stage total/static pressures, temperatures, and densities using the Euler Turbomachinery Equation and Ideal Gas relations.
- **Mass Balance & Blockage:** Incorporates annulus blockage factors to realistically estimate the tip radius and account for boundary layer growth.
- **CAD Geometry Export:** Automatically outputs leading edge ($\beta_1$), trailing edge ($\beta_2$), and camber angles across the radial span. These coordinates can be directly used for lofting 3D blade geometries in software like SolidWorks, Siemens NX, or CATIA.

## ⚙️ Theory & Assumptions

The solver operates under the following aerodynamic assumptions:
* Constant axial velocity ($V_x$) across the stages.
* Free Vortex radial distribution for work equilibrium.
* Fluid is treated as a calorically perfect gas ($\gamma = 1.4$).

## 🛠️ Usage

### Prerequisites
The only required external library is `numpy`.


## 📝 License
This project is licensed under the MIT License.




