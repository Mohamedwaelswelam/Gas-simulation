# Gas-simulation

This repository contains the MATLAB implementation of a steady-state natural gas transmission network simulator. The code is part of a larger integrated energy systems project involving coupled electricity and gas networks.

## 📌 Overview

The simulation models pressure and flow across a multi-node pipeline network using steady-state gas flow equations. It includes:
- Nodal pressure and mass flow rate calculations
- Compressor station modelling
- Numerical convergence checking
- Support for structural and load-based sensitivity test cases

This tool enables analysis of how gas network changes (e.g., pipeline length, diameter, load variations) affect system feasibility and performance.

## 🧠 Main Features
- Models frictional pressure drops using a simplified Weymouth-based formula
- Handles node balancing and pipeline flow direction
- Iterative convergence via fixed-point update
- Built-in compressor behaviour for key branches

## 🗂️ Files
- `NewGasModel.m` – Main function to simulate the gas network

## 🔍 Applications
- Sensitivity analysis of gas system stability
- Integration with power system models via Power-to-Gas (P2G) and Combined Heat and Power (CHP)
- Educational tool for gas flow and convergence behaviour

## 📊 Output
- Prints pressures at each node
- Displays gas flows per branch
- Shows number of iterations and convergence status

## ⚙️ Requirements
- MATLAB R2021a or later (tested)
- No external toolboxes required

## 👨‍🎓 Author
Mohamed Wael Swelam – Final Year Electrical and Electronic Engineering  
University of Manchester  
Supervisor: Dr. Levi Elbaz

## 📜 License
This project is for academic use only.
