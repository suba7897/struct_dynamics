# **Structural Dynamics Tools**

*A MATLAB toolbox for response analysis, vibration mode computation, and numerical integration in structural dynamics.*

---

## **1. Overview**

This repository contains a set of MATLAB scripts developed for the study and analysis of structural dynamics, focusing on single-degree-of-freedom (SDOF) systems and simplified multi-degree-of-freedom (MDOF) approximations. The tools implement classical and numerically robust procedures commonly used in earthquake engineering, modal analysis, and structural vibration research.

The repository is designed for researchers, graduate students, and practitioners who require transparent, reproducible implementations of:

* Duhamel’s integral for SDOF response
* Newmark-β numerical integration
* Response spectrum generation
* Rayleigh and Rayleigh–Stodola mode-shape estimation
* Mass-orthogonal multi-mode extraction

All scripts are accompanied by theoretical documentation in LaTeX format.

---

## **2. Included MATLAB Tools**

### **2.1 Duhamel Integration Methods**

* **`duhamel_simpson_flex.m`**
  Computes the time-history response of an SDOF system using Simpson’s rule for evaluating Duhamel’s integral. Supports mixed unit systems and force/acceleration excitation.

* **`duhamel_response_spectrum.m`**
  Computes displacement, velocity, and acceleration response spectra directly from the Duhamel integral for a user-defined period range.

---

### **2.2 Newmark-β Integration**

* **`newmark_beta_response_spectrum.m`**
  Implements the general Newmark-β method using:
  [
  v_{t+\Delta t},; \dot{v}*{t+\Delta t},; \ddot{v}*{t+\Delta t}
  ]
  Based on arbitrary β and γ parameters (default: average-acceleration method).
  Generates response spectra from numerical integration of the equation of motion.

---

### **2.3 Mode Shape Computation**

* **`rayleigh_modes.m`**
  Computes natural mode shapes using:

  * Classical Rayleigh iteration for the first mode
  * Generalized Rayleigh–Stodola procedure for higher modes
  * Mass-weighted Gram–Schmidt orthogonalization:
    [
    \phi_i^T M \phi_j = 0, \quad i \neq j
    ]
    Produces final converged mode shapes side-by-side for interpretability.

---

## **3. Documentation**

Extensive theoretical explanations, derivations, and worked examples are included in the **`docs/`** directory:

* Rayleigh Method
* Orthogonalization and higher-mode extraction
* Newmark-β theory (general β and γ)
* Duhamel integral theory and numerical evaluation
* Worked examples illustrating algorithmic steps

These documents are provided in LaTeX format for academic use, adaptation, and inclusion in theses.

---

## **4. Usage Instructions**

### **4.1 Cloning the Repository**

```bash
git clone [https://github.com/suba7897/structural_dynamics_tools.git]
```

### **4.2 Adding MATLAB Path**

```matlab
addpath('src');
```

### **4.3 Executing Scripts**

From the MATLAB command window:

```matlab
newmark_beta_response_spectrum
duhamel_simpson_flex
duhamel_response_spectrum
rayleigh_modes
```

Scripts provide prompts for units, damping ratio, periods, and input records (force or acceleration).

---

## **5. Research Motivation**

The tools in this repository are motivated by frequent academic needs in structural dynamics research:

* Generating high-resolution response spectra for parametric studies
* Verifying analytical expressions using numerical integration
* Exploring mode shapes in simplified MDOF systems
* Developing teaching material for graduate dynamics and earthquake engineering courses
* Conducting preliminary studies prior to finite element modeling

The scripts prioritize transparency and clarity over optimization, making them suitable for both education and reproducible research workflows.

---

## **6. Citation**

Users of this repository in academic or professional work are kindly requested to cite it as:

```
Subakaran R. (2025). Structural Dynamics Tools (Version 1.0).
GitHub Repository: https://github.com/suba7897/struct_dynamics
Licensed under CC-BY 4.0.
```

Machine-readable citation data is provided via **`CITATION.cff`**.

---

## **7. License**

This project is licensed under the
**Creative Commons Attribution 4.0 International License (CC-BY 4.0)**.

This license permits use, distribution, adaptation, and modification of the materials for any purpose, provided that appropriate credit is given to the original author.

The full license text is provided in the **`LICENSE`** file.

---

8. Acknowledgments

The author gratefully acknowledges:

Professor Stavroula Pantazapoula, whose teaching, guidance, and course material in structural dynamics strongly influenced the formulation and implementation of the methods in this repository.

The classical literature in structural dynamics, including:

Chopra, Dynamics of Structures

Clough & Penzien, Dynamics of Structures

Newmark (1959), “A Method of Computation for Structural Dynamics”

These works form the theoretical foundation for the algorithms implemented here.

