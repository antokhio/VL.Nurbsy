### (WIP) VL.Nurbsy

A NURBS Bezier library for VL, based on [Stride](https://github.com/stride3d/stride) data types.
Based on [LNLib](https://github.com/BIMCoderLiang/LNLib)


# Transpilation Plan: LNLib to Nurbsy

## Overview
This document outlines the ordered phases for porting the C++ geometry library `LNLib` to C# `Nurbsy`. The order is determined by dependency resolution: implementing low-level mathematical helpers first, followed by core data structures, and finally complex geometric algorithms.

## Phase 1: Infrastructure & Bezier Curves
*Bezier curves are simpler and often serve as building blocks for NURBS operations (e.g., decomposition).*

### 1.1 External Dependencies
Ensure the following helper classes are available or ported:
- `Polynomials` (Specifically `AllBernstein` for Bezier evaluation). [x]
- `ValidationUtils` (For argument checking). [x]
- `XYZ`, `XYZW` (Basic vector types). [x]
- `Intersection`, `Interpolation` (For quadratic curve fitting). [x]/[ ]

### 1.2 BezierCurve Logic (`BezierCurve.h` & `.cpp`)
Port these methods into a `BezierCurveHelper` or directly into the `BezierCurve<T>` struct if appropriate. [x]
1.  **`Check`**: Basic validation of degree vs control points. [x]
2.  **`GetPointOnCurveByBernstein`**: [x]
    - *Dependency:* `Polynomials.AllBernstein`.
    - *Note:* Standard evaluation method.
3.  **`GetPointOnCurveByDeCasteljau`**: [x]
    - *Note:* Numerically stable evaluation method.
4.  **`GetPointOnQuadraticArc`**: [x]
    - *Source:* `BezierCurve.cpp`
    - *Note:* Specific analytic solution for quadratic arcs.
5.  **`ComputerMiddleControlPointsOnQuadraticCurve`**: [x]
    - *Source:* `BezierCurve.cpp`
    - *Dependency:* `Intersection.ComputeRays`, `Interpolation.ComputerWeight...` [x]

---

## Phase 2: NURBS Core Evaluation
*Fundamental operations required to define and evaluate a NURBS curve.*

### 2.1 Basic Evaluation
1.  **`Check`**: Validate Knot vector and Control Point counts ($m = n + p + 1$). [x]
2.  **`GetPointOnCurve`**: [x]
    - *Dependency:* `BsplineCurve` (Basis function evaluation).
3.  **`ComputeRationalCurveDerivatives`**:
    - *Dependency:* `BsplineCurve` derivatives.
4.  **`IsLinear`**: Geometric check for linearity.
5.  **`CreateLine`**: Factory method.
6.  **`CreateTransformed`**: Applies a 4x4 matrix to control points.
7.  **`Reverse`**: Reverses control points and knots.
8.  **`Reparametrize`**: Scales the knot vector to a new domain $[a, b]$.

---

## Phase 3: Knot & Degree Manipulation
*Operations that modify the underlying data structure (Knots/ControlPoints) without changing the curve's shape (mostly).*

1.  **`InsertKnot`**: Fundamental for many algorithms.
2.  **`RefineKnotVector`**: Inserts multiple knots.
3.  **`RemoveKnot`**: Inverse of insertion.
4.  **`RemoveExcessiveKnots`**: Cleans up the knot vector.
5.  **`ElevateDegree`**: Increases degree (adds control points).
6.  **`ReduceDegree`**: Decreases degree (approximation).
7.  **`GetPointOnCurveByCornerCut`**: Visualization/approximation helper.

---

## Phase 4: Topological Operations
*Splitting, merging, and decomposing curves.*

1.  **`SplitAt`**: Cuts the curve at a parameter.
    - *Dependency:* `RefineKnotVector`.
2.  **`Segment`**: Extracts a sub-curve.
3.  **`Merge`**: Joins two curves.
4.  **`DecomposeToBeziers`**: Converts NURBS into a list of Bezier curves.
    - *Critical for:* Length calculation and intersection algorithms.
5.  **`IsClosed`**, **`IsClamp`**, **`IsPeriodic`**.
6.  **`ToClampCurve`** / **`ToUnclampCurve`**.

---

## Phase 5: Geometric Properties
*Analyzing the curve's differential properties.*

1.  **`CanComputeDerivative`**: Safety check.
2.  **`Curvature`**: $k = \frac{|r' \times r''|}{|r'|^3}$
3.  **`Normal`**: Principal normal vector.
4.  **`Torsion`**.
5.  **`ProjectNormal`**: Projects a point onto the curve plane (if planar).

---

## Phase 6: Geometric Algorithms & Interpolation
*Constructing curves from points.*

1.  **`GlobalInterpolation`**: Fits a curve through a set of points.
    - *Dependency:* Linear System Solver (LU Decomposition or SVD).
2.  **`LeastSquaresApproximation`**: Fits a curve to many points (approximate).
3.  **`CreateArc`**, **`CreateOpenConic`**: Analytic construction of conics.
4.  **`Tessellate`**: Discretizes the curve into polyline.
    - *Internal Helper:* `TessellateCore`.

---

## Phase 7: Inverse Evaluation & Length
*Iterative methods for finding parameters or lengths.*

1.  **`EquallyTessellate`**: Generates evenly spaced points.
2.  **`GetParamOnCurve(point)`**: Projection.
    - *Strategy:* Use `EquallyTessellate` for initial guess, then Newton-Raphson.
3.  **`ApproximateLength`**:
    - *Dependency:* `DecomposeToBeziers`, `FirstDerivativeLengthFunction` (Internal class), Simpson's Rule integrator.
4.  **`GetParamOnCurve(length)`**: Finds parameter $t$ such that length from start is $L$.

---

## Phase 8: Advanced Modification (Optional)
*Complex shape editing tools.*

1.  **`Offset`**: Creates a parallel curve.
    - *Internal Helpers:* `OffsetByTillerAndHanson`, `OffsetByPieglAndTiller`.
2.  **`Flattening`**, **`Bending`**, **`Warping`**.
3.  **`FitWithConic`**, **`FitWithCubic`**.