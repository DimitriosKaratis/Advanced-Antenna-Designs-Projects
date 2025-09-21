# Advanced Antenna Designs Projects

This repository contains the implementation of two projects in **Advanced Antenna Designs**, focusing on antenna modeling, simulation, and beamforming techniques.  
The work was developed as part of the *Advanced Antenna Designs* course at the **Department of Electrical and Computer Engineering, Aristotle University of Thessaloniki**.

---

## 📡 Project 1: Dual-Band H-Ring Antenna Design

- **Topic:** Design and simulation of a **compact H-ring microstrip-fed monopole antenna** with dual-band operation for RFID and wireless sensor applications.  
- **Reference Paper:**  
  *Nasser Ojaroudi and Mohammad Ojaroudi, "Compact H-Ring Antenna with Dual-Band Operation for Wireless Sensors and RFID Tag Systems in ISM Frequency Bands,"* Microwave Opt Technol Lett 55:697–700, 2013.  
- **Objective:**  
  - Reproduce the antenna geometry and results described in the paper.  
  - Simulate the antenna using MATLAB Antenna Toolbox.  
  - Validate performance metrics:  
    - Return loss (S11)  
    - Radiation patterns  
    - Resonance frequencies (2.45 GHz and 5.8 GHz)  
- **Implementation:**  
  - MATLAB scripts for antenna modeling and simulation.  
  - Comparison of results with the referenced paper.  
  - PDF report including methodology, code, and diagrams.

---

## 📶 Project 2: Beamforming with Linear Arrays and Optical Phased Arrays

- **Topic:** Design and simulation of **beamformers** for uniform linear antenna arrays and optical phased arrays.  
- **Objectives:**  
  1. **Uniform Linear Array (ULA) Beamforming**  
     - Simulation of 5, 7, and 9 element arrays at 1 GHz.  
     - Beam steering for different angles (0°–180°).  
     - Effect of element spacing (λ/2 vs λ/4).  
     - Application of binomial amplitude distribution for sidelobe suppression.  
  2. **Optical Phased Array Beamforming**  
     - Simulation of beam steering in an **Optical Phased Array (OPA)**.  
     - Modeling of Mach–Zehnder interferometers and phase shifters.  
     - Demonstration of beam steering at 30°, 60°, and 90°.  
     - Qualitative discussion of implementing binomial amplitude distribution in OPA networks.  
- **Implementation:**  
  - MATLAB scripts for beamforming factor and field distribution calculations.  
  - Polar radiation diagrams for different array configurations.  
  - PDF report including methodology, code, and diagrams.
