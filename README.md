# Vandermonde Polynomial Solver

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

## Technical Implementation

### Libraries Used

| Library | Purpose |
|---------|---------|
| `numpy` | Matrix operations, solving linear systems |
| `ipywidgets` | Interactive UI elements |
| `plotly` | Interactive plots |
| `MathJax` | Display mathematical equations |

### Key Functions

**1. Vandermonde Matrix Construction**
```python
V = np.vander(X, increasing=True)
```

**2. Solving the System**
```python
coefficients = np.linalg.solve(V, Y)
```

**3. Polynomial Evaluation**
```python
Y = np.polyval(coefficients, T)
```

**4. Polynomial Expansion (Binomial Theorem)**
```python
for i, c in enumerate(coeffs):
    for j in range(i + 1):
        result[j] += c * comb(i,j) * ((-shift) ** (i - j))
```

---

## User Interface Guide

### Buttons & Controls

| Element | Purpose |
|---------|---------|
| **Number of Points** | Set how many data points (2-10) |
| **Generate Table** | Create input table with that many rows |
| **Shift Method** | Choose Mean/First/None for stability |
| **Solve Polynomial** | Run the calculation and show results |
| **Evaluate** | Test the polynomial at any T value |

### Output Sections

| Section | Shows |
|---------|-------|
| **Step 1-2** | Assumed polynomial and equations |
| **Step 3-4** | Vandermonde matrices |
| **Step 5-6** | Shift and shifted matrix |
| **Step 7-8** | Solution coefficients |
| **Step 9-10** | Final polynomial |
| **Step 11** | Verification table |
| **Step 12** | Interactive evaluation |
| **Plot** | Interactive graph |

---

## Common Questions

### Q: Why do we use shifting?
**A:** To make the matrix well-conditioned when T values are large and close together. This prevents numerical errors.

### Q: What does "condition number" mean?
**A:** It measures how sensitive the solution is to small changes in input. Lower is better.

### Q: How many points can I use?
**A:** 2 to 10 points. More points = higher degree polynomial.

### Q: Can I use any T values?
**A:** Yes, but for best results use the "Mean" shift method when T values are large.

### Q: Why do coefficients have alternating signs?
**A:** This is normal for ill-conditioned systems. It's why we use shifting.

### Q: What is a "good" error?
**A:** Error < 10⁻⁶ is good. Error < 10⁻¹⁰ is perfect (machine precision).

---

## Mathematical Summary

### Full Process

```
Input: (T₁,Y₁), (T₂,Y₂), ..., (Tₙ,Yₙ)
        ↓
Apply shift: x = T - shift
        ↓
Build Vandermonde: V = [x⁰, x¹, x², ..., xⁿ⁻¹]
        ↓
Solve: V·a = Y → find coefficients a₀, a₁, ..., aₙ₋₁
        ↓
Expand using binomial theorem
        ↓
Output: P(T) = aₙTⁿ + aₙ₋₁Tⁿ⁻¹ + ... + a₁T + a₀
        ↓
Verify: Calculate Y at each T, check errors
        ↓
Evaluate: Test polynomial at any T value
```

---

## Example Walkthrough

### Input
```
Points: 4
T: [30.75, 30.88, 31.00, 31.12]
Y: [867.2295, 888.5687, 875.7651, 885.1544]
Shift Method: Mean
```

### Output (Final Polynomial)
```
P(T) = 5010.7141660891T³ - 465225.8306407345T² 
       + 14398026.1783565581T - 148530098.2401689589
```

### Verification
```
Maximum Error: 6.19e-08 (very small!)
Condition Number: 8.11e+02 (very stable!)
Status: ✓ Good reconstruction
```

### Evaluation Example
```
At T = 30.95
Y = 881.245672...
```

---

## Credits

**Developed by:** MrRogueKnight

**GitHub Repository:** [MrRogueKnight/Vandermonde-Polynomial-Solver](https://github.com/MrRogueKnight/Vandermonde-Polynomial-Solver)

**Kaggle Notebook:** [View on Kaggle](https://www.kaggle.com/code/your-username/your-notebook)

---
