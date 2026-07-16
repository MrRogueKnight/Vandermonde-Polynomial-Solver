# Vandermonde Polynomial Interpolation & Extrapolation Framework

## Overview

A comprehensive scientific computing framework for polynomial interpolation and adaptive extrapolation using shifted and scaled Vandermonde matrices. This implementation provides numerical stability, extensive diagnostics, and professional-grade visualization.

**Key Features:**
- Shifted and scaled Vandermonde basis for numerical stability
- Exact interpolation and least-squares fitting
- Adaptive extrapolation with polynomial capping
- Comprehensive numerical diagnostics
- Publication-quality plots with Matplotlib
- Interactive Jupyter widgets interface

---

## Table of Contents

- [What is Polynomial Interpolation?](#what-is-polynomial-interpolation)
- [The Vandermonde Matrix Method](#the-vandermonde-matrix-method)
- [The Numerical Stability Problem](#the-numerical-stability-problem)
- [Extrapolation Framework](#extrapolation-framework)
- [Visualization Engine](#visualization-engine)
- [Installation](#installation)
- [Usage Guide](#usage-guide)
- [API Reference](#api-reference)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)

---

## What is Polynomial Interpolation?

Polynomial interpolation finds a polynomial that passes through a given set of points.

**The Problem:**
Given points: `(T₁, Y₁), (T₂, Y₂), ..., (Tₙ, Yₙ)`

Find a polynomial of degree `n-1` that passes through all points:

```
Y = aₙTⁿ + aₙ₋₁Tⁿ⁻¹ + ... + a₁T + a₀
```

**Example:**
Given 4 points:
```
(30.75, 867.2295), (30.88, 888.5687), (31.00, 875.7651), (31.12, 885.1544)
```

We find a cubic polynomial:
```
Y = a₃T³ + a₂T² + a₁T + a₀
```

---

## The Vandermonde Matrix Method

### Step 1: Write Equations
For each point, substitute T and Y:

```
Point 1: a₃(30.75)³ + a₂(30.75)² + a₁(30.75) + a₀ = 867.2295
Point 2: a₃(30.88)³ + a₂(30.88)² + a₁(30.88) + a₀ = 888.5687
Point 3: a₃(31.00)³ + a₂(31.00)² + a₁(31.00) + a₀ = 875.7651
Point 4: a₃(31.12)³ + a₂(31.12)² + a₁(31.12) + a₀ = 885.1544
```

### Step 2: Matrix Form
This becomes a matrix equation:

```
[30.75³  30.75²  30.75  1]   [a₃]   [867.2295]
[30.88³  30.88²  30.88  1] × [a₂] = [888.5687]
[31.00³  31.00²  31.00  1]   [a₁]   [875.7651]
[31.12³  31.12²  31.12  1]   [a₀]   [885.1544]
```

The matrix is called the **Vandermonde Matrix**.

### Step 3: Solve
We solve the system using `numpy.linalg.solve()` to find `a₃, a₂, a₁, a₀`.

---

## The Numerical Stability Problem

### Why is it a Problem?
When T values are **large and close together** (like 30.75, 30.88, 31.00), the matrix becomes **ill-conditioned**.

This means:
- Small rounding errors (10⁻¹⁶) get amplified massively
- Coefficients become huge with alternating signs
- Example: `a₃ = 5010`, `a₂ = -465226`, `a₁ = 14398050`, `a₀ = -148530345`

### The Solution: Shifting
Instead of T, we use:
```
x = T - shift
```

Where `shift` is a number close to T (like the mean).

### Example
```
shift = 30.9375 (mean of 30.75, 30.88, 31.00, 31.12)

x₁ = 30.75 - 30.9375 = -0.1875
x₂ = 30.88 - 30.9375 = -0.0575
x₃ = 31.00 - 30.9375 = 0.0625
x₄ = 31.12 - 30.9375 = 0.1825
```

Now x-values are **small and centered around zero**, making the matrix well-conditioned.

### Shift Methods

| Method | Description | When to Use |
|--------|-------------|-------------|
| **Mean** | shift = average of all T values | Recommended for most cases |
| **First** | shift = first T value | When you want T₁ to become 0 |
| **None** | shift = 0 | Only for small T values (like 1,2,3) |

---

## The Expansion Process

### Step 1: Solve Shifted Polynomial
We solve for P(x) where `x = T - shift`:
```
P(x) = a₀ + a₁x + a₂x² + a₃x³
```

### Step 2: Expand Using Binomial Theorem
Substitute `x = T - shift`:
```
P(T) = a₀ + a₁(T-shift) + a₂(T-shift)² + a₃(T-shift)³
```

### Step 3: Expand Each Term
```
(T-shift)² = T² - 2shift·T + shift²
(T-shift)³ = T³ - 3shift·T² + 3shift²·T - shift³
```

### Step 4: Collect Terms
Finally get:
```
P(T) = A₃T³ + A₂T² + A₁T + A₀
```

Where:
```
A₃ = a₃
A₂ = a₂ - 3·shift·a₃
A₁ = a₁ - 2·shift·a₂ + 3·shift²·a₃
A₀ = a₀ - shift·a₁ + shift²·a₂ - shift³·a₃
```

---

## Condition Number

The **condition number** tells us how stable the solution is:

- **Condition Number < 10⁶** → Stable solution ✓
- **Condition Number > 10⁶** → Potentially unstable
- **Condition Number > 10¹²** → Ill-conditioned

Our solver shows the condition number after solving.

---

## Verification

We check the solution by:
1. Evaluating P(T) at each original T value
2. Comparing with original Y values
3. Calculating the error

### Error Types
| Error Range | Status |
|-------------|--------|
| < 10⁻¹⁰ | Perfect reconstruction (machine precision) |
| < 10⁻⁶ | Good reconstruction |
| > 10⁻⁶ | Large error detected |

---

## Extrapolation Framework

The framework includes adaptive extrapolation with multiple methods:

### Extrapolation Methods

| Method | Description | Best For |
|--------|-------------|----------|
| **Linear** | Uses last two points | Large extrapolation distances |
| **Polynomial** | Uses interpolated polynomial | Small extrapolation distances |
| **Rational** | Rational function approximation | Asymptotic behavior |
| **Ensemble** | Blends multiple methods | Moderate extrapolation |

### Method Selection Logic

```
Distance < 0.1x  → Polynomial blend (70% polynomial + 30% linear)
Distance < 0.5x  → Ensemble (weighted average)
Distance > 0.5x  → Linear-rational blend (70% linear + 30% rational)
```

### Reliability Scoring

```
Reliability = 0.7 - 0.04 × Extrapolation Distance
Score > 0.8  → HIGH reliability
Score > 0.5  → MEDIUM reliability
Score > 0.2  → LOW reliability
Score < 0.2  → VERY LOW reliability
```

---

## Visualization Engine

Publication-quality plots with Matplotlib:

### Plot Types

1. **Interpolation Plot**
   - Data points with labels
   - Interpolating polynomial
   - Residual analysis

2. **Extrapolation Plot**
   - Full range visualization
   - Confidence intervals
   - Method comparison

3. **Reliability Gauge**
   - Circular gauge display
   - Reliability percentage
   - Extrapolation distance

4. **Comprehensive Analysis**
   - Combined plots
   - Method comparison table
   - Residual analysis

### Example Plots

```python
# Generate interpolation plot
fig = interpolator.plot_interpolation()
plt.show()

# Generate extrapolation plot
fig = interpolator.plot_extrapolation(x_target=33.31)
plt.show()

# Generate comprehensive plot
fig = interpolator.plot_comprehensive(x_targets=[32.0, 33.31, 35.0])
plt.show()
```

---

## Installation

### Requirements

```bash
Python 3.8+
numpy
matplotlib
ipywidgets
jupyter
```

### Install from Source

```bash
git clone https://github.com/MrRogueKnight/Vandermonde-Polynomial-Solver.git
cd Vandermonde-Polynomial-Solver
pip install -r requirements.txt
```

### Quick Start in Jupyter

```python
from interpolation_framework import InterpolationApplication
app = InterpolationApplication()
```

---

## Usage Guide

### Interactive Jupyter Application

```python
from interpolation_framework import InterpolationApplication
app = InterpolationApplication()
```

The application provides:
- Data input table with up to 10 points
- Shift method selection (mean/first/none)
- Least squares fitting option
- Tabbed results display
- Extrapolation with visualization
- Method comparison

### Programmatic Usage

```python
from interpolation_framework import ShiftedVandermondeInterpolator
import numpy as np

# Create interpolator
interp = ShiftedVandermondeInterpolator()

# Sample data
x = np.array([30.75, 30.88, 31.00, 31.12])
y = np.array([867.2295, 888.5687, 875.7651, 885.1544])

# Interpolate
result = interp.interpolate(x, y)

# Evaluate at points
y_pred = interp.evaluate(np.array([30.95, 31.05]))

# Extrapolate
extrap_result = interp.extrapolate(33.31)

# Generate plots
fig = interp.plot_interpolation()
fig = interp.plot_extrapolation(33.31)
```

### Extrapolation Example

```python
from interpolation_framework import ShiftedVandermondeInterpolator

interp = ShiftedVandermondeInterpolator()
interp.interpolate(x, y)

# Extrapolate at multiple points
results = []
for xt in [32.0, 33.31, 35.0]:
    result = interp.extrapolate(xt)
    results.append(result)
    print(f"X = {xt:.2f}: Y = {result.y_predicted:.4f}")
    print(f"Reliability: {result.reliability_score*100:.1f}%")
    print(f"Method: {result.method_used}")
```

---

## API Reference

### `ShiftedVandermondeInterpolator`

**Parameters:**
- `shift_method` (ShiftMethod): MEAN, FIRST, or NONE
- `scale_data` (bool): Apply scaling after shifting
- `expansion_mode` (ExpansionMode): SHIFTED, EXPANDED, or BOTH

**Methods:**
- `interpolate(x_data, y_data)`: Exact interpolation
- `fit(x_data, y_data, degree)`: Least squares fitting
- `evaluate(x_values, use_shifted)`: Evaluate polynomial
- `extrapolate(x_target)`: Adaptive extrapolation
- `plot_interpolation(x_range, title)`: Generate interpolation plot
- `plot_extrapolation(x_target, x_range, title)`: Generate extrapolation plot
- `plot_comprehensive(x_targets, title)`: Generate comprehensive plot

### `ExtrapolationResult`

**Attributes:**
- `x_target`: Target x value
- `y_predicted`: Predicted y value
- `method_used`: Extrapolation method
- `method_predictions`: All method predictions
- `reliability_score`: 0-1 reliability score
- `reliability_level`: Categorical reliability
- `extrapolation_distance`: Normalized distance
- `confidence_interval`: 95% confidence interval
- `warning`: Warning message if any

### `InterpolationResult`

**Attributes:**
- `shifted_coefficients`: Coefficients in shifted basis
- `expanded_coefficients`: Coefficients in original basis
- `shift_value`: Shift value used
- `scale_value`: Scale value used
- `degree`: Polynomial degree
- `diagnostics`: Numerical diagnostics
- `verification`: Verification metrics

---

## Examples

### Example 1: Basic Interpolation

```python
import numpy as np
from interpolation_framework import ShiftedVandermondeInterpolator

# Data
x = np.array([30.75, 30.88, 31.00, 31.12])
y = np.array([867.2295, 888.5687, 875.7651, 885.1544])

# Interpolate
interp = ShiftedVandermondeInterpolator()
result = interp.interpolate(x, y)

# Print polynomial
print("Shifted Form:")
print(interp.get_expanded())
```

### Example 2: Extrapolation with Visualization

```python
# Extrapolate
result = interp.extrapolate(33.31)

print(f"Prediction: {result.y_predicted:.6f}")
print(f"Reliability: {result.reliability_score*100:.1f}%")
print(f"Method: {result.method_used}")

# Generate plot
fig = interp.plot_extrapolation(33.31)
plt.show()
```

### Example 3: Multiple Extrapolations

```python
x_targets = [32.0, 33.31, 35.0, 39.5]
results = []

for xt in x_targets:
    result = interp.extrapolate(xt)
    results.append(result)
    print(f"X={xt:.2f}: Y={result.y_predicted:.4f} "
          f"(Reliability: {result.reliability_score*100:.1f}%)")

# Generate comprehensive plot
fig = interp.plot_comprehensive(x_targets)
plt.show()
```

### Example 4: Least Squares Fitting

```python
# Generate noisy data
x = np.linspace(0, 10, 20)
y = 2*x + 1 + 0.5*np.random.randn(20)

# Fit with degree 1 (linear regression)
interp = ShiftedVandermondeInterpolator()
result = interp.fit(x, y, degree=1)

print(f"R²: {result.verification.r_squared:.4f}")
print(f"RMSE: {result.verification.rmse:.4f}")
```

---

## Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md).

### Development Setup

```bash
git clone https://github.com/MrRogueKnight/Vandermonde-Polynomial-Solver.git
cd Vandermonde-Polynomial-Solver
pip install -e ".[dev]"
pytest tests/
```

### Pull Request Process

1. Fork the repository
2. Create a feature branch
3. Add your changes with tests
4. Update documentation
5. Submit pull request

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

- NumPy for linear algebra operations
- Matplotlib for visualization
- Jupyter for interactive computing

---

## Contact

**Developer:** MrRogueKnight

**GitHub:** [MrRogueKnight/Vandermonde-Polynomial-Solver](https://github.com/MrRogueKnight/Vandermonde-Polynomial-Solver)

**Kaggle:** [View Notebook](https://www.kaggle.com/code/your-username/your-notebook)

---

## Citation

If you use this framework in your research, please cite:

```bibtex
@software{VandermondeFramework2024,
  author = {MrRogueKnight},
  title = {Vandermonde Polynomial Interpolation and Extrapolation Framework},
  year = {2024},
  url = {https://github.com/MrRogueKnight/Vandermonde-Polynomial-Solver}
}
```

---

## Version History

### v2.0.0 (Current)
- Added adaptive extrapolation framework
- Added visualization engine
- Added confidence intervals
- Added reliability scoring
- Added method comparison
- Added comprehensive diagnostics

### v1.0.0
- Initial release
- Shifted Vandermonde interpolation
- Numerical diagnostics
- Basic evaluation

---

## FAQ

### Q: Why do we use shifting?
**A:** To make the matrix well-conditioned when T values are large and close together. This prevents numerical errors.

### Q: What does "condition number" mean?
**A:** It measures how sensitive the solution is to small changes in input. Lower is better.

### Q: How many points can I use?
**A:** 2 to 10 points in the interactive UI. Programmatically, any number of points can be used.

### Q: Can I use any T values?
**A:** Yes, but for best results use the "Mean" shift method when T values are large.

### Q: Why do coefficients have alternating signs?
**A:** This is normal for ill-conditioned systems. It's why we use shifting.

### Q: What is a "good" error?
**A:** Error < 10⁻⁶ is good. Error < 10⁻¹⁰ is perfect (machine precision).

### Q: How reliable are extrapolations?
**A:** Reliability depends on distance from data. The framework provides quantitative reliability scores and confidence intervals.

### Q: What methods are available for extrapolation?
**A:** Linear, Polynomial, Rational, and Ensemble methods with automatic selection based on extrapolation distance.
