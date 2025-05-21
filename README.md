# 🧮 Barrier Options Pricing using Finite Difference Method

This project implements **finite difference methods (FDM)** for pricing **barrier options**, including single and double barrier types, under both **European** and **American** exercise styles. The numerical models are developed and tested in MATLAB.

---

## 📌 Project Overview

Barrier options are exotic options whose payoff depends on whether the underlying asset price breaches a specified barrier level during the option's life. This project explores the pricing of:

- **Single barrier options**: Up-and-in, Up-and-out, Down-and-in, Down-and-out
- **Double barrier options**: Knock-in and Knock-out types
- **Exercise styles**: European and American
- **Numerical method**: Finite Difference Method (FDM) with Crank-Nicolson scheme

The models are designed to accommodate non-uniform spatial grids and varying time step resolution to handle barrier discontinuities more accurately.

---

## 🛠 Features

- Support for **European and American** barrier options
- **Single and double barrier** option pricing
- Modular and object-oriented MATLAB design
- Comparison with **Black-Scholes** and **Binomial Tree** benchmarks
- Integrated visualization and tabular output
- Configurable grid resolution and Crank-Nicolson parameter `θ`

---

## 📁 Project Structure

```

barrier_options/
├── ame_test_double_barriers.m
├── ame_test_single_barriers.m
├── ame_test_vanilla.m
├── eur_test_double_barriers.m
├── eur_test_single_barriers.m
├── eur_test_vanilla.m
├── gentable.m
├── option.m
├── option_new.m
├── README.md
├── utils.m
└── visualization.m
└── README.md # This file

```

---

## 🚀 How to Run

1. Open MATLAB.
2. Navigate to the `barrier_options/` folder.
3. Run one of the test scripts:

   ```matlab
   ame_test_double_barriers
   ```

4. Results will be visualized and printed in the console.

---

## 📈 Methods Used

- **Finite Difference Method (FDM)**:

  - Crank-Nicolson (θ = 0.5)
  - Non-uniform spatial grids for better accuracy near barriers

- **Binomial Tree**:

  - Used for American vanilla option benchmark

- **Black-Scholes Formula**:

  - Used for European vanilla options

---

## 👥 Team: M.A.N.S

We are a passionate team of students dedicated to computational finance and numerical methods:

- **Arsen**
- **Asyat**
- **Abay**
- **Zhasulan**
- **Esbol**

Project Team Name: **M.A.N.S**

---

## 📚 References

- Hull, J. C. _Options, Futures, and Other Derivatives_
- Wilmott, P. _Quantitative Finance_
- MATLAB documentation for PDEs and numerical schemes

---

## 📌 License

This project is intended for educational and academic use. Feel free to modify and adapt for your learning or research.
