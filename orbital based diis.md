# Orbital-Space DIIS for Hartree–Fock and SCF Methods
## Correct Mathematical Formulation and Implementation Guide

**Author:** Fubin Zhao  
**Date:** \today

## Abstract

The direct inversion in the iterative subspace (DIIS) is a widely used convergence accelerator for self-consistent field (SCF) calculations. While the conventional Fock-matrix DIIS works well for many cases, an orbital-space formulation offers superior stability and a larger radius of convergence, especially for multi-configurational wave functions. This document presents the correct orbital-space DIIS algorithm for Hartree–Fock (HF) and general SCF methods, emphasizing that all rotation angles must be defined with respect to a *fixed reference orbital set*. A common pitfall—storing incremental rotations between successive iterations—is identified and corrected. The proper use of matrix logarithms, error vectors based on total-angle differences, and the removal of redundant intra-subspace rotations are described in detail. Pseudocode and practical recommendations are provided to help implement a robust, efficient, and reference‑consistent DIIS accelerator. The presentation is based on the seminal work of Ionova and Carter [1] and on standard SCF theory.

## 1. Introduction

The self-consistent field (SCF) iteration is at the heart of most electronic structure calculations. Pulay's DIIS [2] accelerates convergence by minimizing a linear combination of error vectors from previous iterations. In single-determinant Hartree–Fock (HF), the error vector is usually constructed from the antisymmetric part of the Fock matrix in the molecular orbital (MO) basis [3,4]. This so-called “Fock-matrix DIIS” (or “CFM-DIIS”) is simple and effective, but it can fail when the initial guess is far from the solution.

Alternatively, the DIIS idea can be applied directly to the orbitals themselves, as first proposed by Ionova and Carter for the generalized valence bond (GVB) method [1]. The orbital-space DIIS has a larger radius of convergence and can be used with any SCF method that updates orbitals iteratively, including HF, MCSCF, and Kohn–Sham DFT. However, its implementation is subtle: a naive linear combination of orbital coefficients breaks orthonormality. The solution is to work with antisymmetric rotation matrices $\kappa$ that parameterize unitary transformations of a fixed reference orbital set.

A common mistake in practical codes is to store the *incremental* rotation from the previous iteration to the current one and then treat these increments as if they were total angles relative to a fixed reference. This leads to inconsistent reference frames and usually prevents convergence. The purpose of this document is to present the correct algorithm in a clear, self-contained manner, focusing on HF but with the understanding that it generalizes to any SCF.

## 2. Orbital Parametrization and Reference Frame

### 2.1 Fixed Reference Orbitals
Let $C_{\text{ref}} \in \mathbb{R}^{N \times N}$ be a set of orthonormal molecular orbitals (coefficient matrix) chosen once at the beginning of the SCF procedure. Typically $C_{\text{ref}}$ is the initial guess (e.g., Hückel orbitals). It satisfies $C_{\text{ref}}^T S C_{\text{ref}} = I$, where $S$ is the atomic orbital overlap matrix. **This reference set remains fixed throughout the DIIS process.**

### 2.2 Total Rotation Angle Matrix
At iteration $n$, the orbitals $C^{(n)}$ are related to $C_{\text{ref}}$ by a unitary transformation:
$$
C^{(n)} = C_{\text{ref}} \, U^{(n)}, \qquad U^{(n)} = \exp\bigl(\kappa^{(n)}\bigr),
$$
where $\kappa^{(n)}$ is an antisymmetric matrix ($\kappa^T = -\kappa$), called the **total rotation angle matrix**. The matrix exponential is defined by the power series and yields an orthogonal (or unitary) matrix. The inverse operation, obtaining $\kappa^{(n)}$ from $C^{(n)}$, is given by the matrix logarithm:
$$
\kappa^{(n)} = \ln\!\bigl(C_{\text{ref}}^T S C^{(n)}\bigr).
$$
Because $C_{\text{ref}}^T S C^{(n)}$ is orthogonal, its logarithm exists (up to branch choices, which are irrelevant for small angles near convergence).

### 2.3 Subspace Partition and Redundant Rotations
For Hartree–Fock, we partition the molecular orbitals into two subspaces: occupied (O) and virtual (V). The antisymmetric matrix $\kappa$ then has the block structure
$$
\kappa = \begin{pmatrix}
0_{OO} & \kappa_{OV} \\
-\kappa_{OV}^T & 0_{VV}
\end{pmatrix},
$$
i.e., the diagonal blocks ($OO$ and $VV$) are zero. Rotations within the occupied subspace (or within the virtual subspace) do not change the total energy or the electron density; they are gauge (redundant) degrees of freedom. Consequently, we can *enforce* the condition that these blocks are zero. In practice, the $\kappa$ computed from Eq. (2) may have small non‑zero diagonal blocks due to numerical noise or different gauge choices; we will project them out.

## 3. Error Vector and DIIS Extrapolation

### 3.1 Proper Error Vector
The error vector should measure the physical change between two consecutive iterations. Because all $\kappa^{(n)}$ are defined relative to the same $C_{\text{ref}}$, the difference $\kappa^{(n)} - \kappa^{(n-1)}$ directly represents the step in Lie algebra. Removing the redundant intra‑subspace components, we define
$$
\mathbf{e}_n = \mathcal{P}_{\text{inter}}\bigl(\kappa^{(n)} - \kappa^{(n-1)}\bigr), \qquad (n \ge 2)
$$
where $\mathcal{P}_{\text{inter}}$ sets the $OO$ and $VV$ blocks to zero, keeping only the off-diagonal $OV$ block. For the first iteration ($n=1$), where no previous $\kappa$ exists, we use
$$
\mathbf{e}_1 = \mathcal{P}_{\text{inter}}\bigl(\kappa^{(1)}\bigr).
$$

### 3.2 Why Not Use Incremental Rotations?
Many SCF codes (including standard HF implementations) produce a rotation that takes the current orbitals to the improved ones:
$$
C^{(n)} = C^{(n-1)} e^{\delta\kappa_n}.
$$
Here $\delta\kappa_n$ is an *incremental* rotation (often obtained from the orbital gradient or from a trust‑radius step). If we store $\delta\kappa_n$ and later try to form a linear combination $\sum c_i \delta\kappa_i$, we are mixing rotations from different reference frames ($C^{(n-1)}$ changes each step). The resulting matrix does not correspond to any meaningful unitary transformation, and the DIIS extrapolation fails. Therefore, **one must always work with total angles $\kappa^{(n)}$ relative to a fixed $C_{\text{ref}}$**.

### 3.3 DIIS Extrapolation
Assume we have stored $d$ history pairs $(\kappa^{(i)}, \mathbf{e}_i)$ for $i = i_0, \dots, i_0+d-1$. The Gram matrix is
$$
B_{ij} = \langle \mathbf{e}_i, \mathbf{e}_j \rangle = \sum_{p=1}^{D_{\text{rot}}} e_i(p) e_j(p), \quad i,j=1,\dots,d,
$$
where $D_{\text{rot}} = N_{\text{occ}} N_{\text{virt}}$ is the number of independent parameters (the size of the $OV$ block). The augmented linear system is
$$
\begin{bmatrix}
B_{11} & \cdots & B_{1d} & -1 \\
\vdots & \ddots & \vdots & \vdots \\
B_{d1} & \cdots & B_{dd} & -1 \\
-1 & \cdots & -1 & 0
\end{bmatrix}
\begin{bmatrix}
c_1 \\ \vdots \\ c_d \\ \lambda
\end{bmatrix}
=
\begin{bmatrix}
0 \\ \vdots \\ 0 \\ -1
\end{bmatrix},
\qquad \text{(1)}
$$
with the constraint $\sum_{i=1}^d c_i = 1$. Solving this system yields the coefficients $c_i$. The extrapolated total angle matrix is then
$$
\kappa_{\text{new}} = \sum_{i=1}^d c_i \, \kappa^{(i)}.
$$

### 3.4 Orbital Update
The corresponding unitary transformation is
$$
U_{\text{new}} = \exp(\kappa_{\text{new}}),
$$
and the new orbitals are obtained by rotating the fixed reference:
$$
C_{\text{new}} = C_{\text{ref}} \, U_{\text{new}}.
$$
Because $\kappa_{\text{new}}$ is antisymmetric, $U_{\text{new}}$ is orthogonal, and $C_{\text{new}}$ automatically satisfies $C_{\text{new}}^T S C_{\text{new}} = I$.

## 4. Algorithm Pseudocode

```python
# Initialization
C_ref = C_initial                # fixed reference orbitals (e.g., Huckel)
history = []                     # list of (kappa, error_vector)
S = overlap_matrix               # AO overlap

for n = 1, 2, ...:
    # 1. Perform one SCF iteration (e.g., build Fock, diagonalize)
    C_new = scf_iteration(C_current)    # improved orbitals
    
    # 2. Compute total rotation angle from C_ref to C_new
    M = C_ref.T @ S @ C_new      # orthogonal matrix
    kappa = logm(M)              # matrix logarithm (antisymmetric)
    
    # 3. Build error vector (project out redundant OO and VV blocks)
    if n == 1:
        e = project_off_diagonal(kappa)
    else:
        delta = kappa - kappa_prev
        e = project_off_diagonal(delta)
    
    # 4. Store history (keep at most max_d items)
    history.append( (kappa, e) )
    if len(history) > max_d:
        history.pop(0)
    
    # 5. DIIS extrapolation if at least 2 history items exist
    if len(history) >= 2:
        try:
            d = len(history)
            # Build Gram matrix from error vectors
            B = [[np.dot(ei, ej) for ej in e_list] for ei in e_list]
            # Solve augmented system (Eq. (1))
            coeffs = solve_diis(B)       # returns list c_i
            kappa_extrap = sum(c_i * k_i for (k_i,_), c_i in zip(history, coeffs))
            U_new = expm(kappa_extrap)   # matrix exponential
            C_extrap = C_ref @ U_new
            C_next = C_extrap
        except np.linalg.LinAlgError:    # ill-conditioned
            # Fallback: reduce history or use unscaled update
            if len(history) > 2:
                history.pop(0)           # drop oldest and retry next iteration
            C_next = C_new
    else:
        C_next = C_new
    
    # 6. Prepare for next iteration
    kappa_prev = kappa
    C_current = C_next
    
    # 7. Convergence check
    if np.linalg.norm(kappa - kappa_prev) < tol:
        break