# FCLIB: LOGEM — Principal Matrix Logarithm

## 1. Algorithm Principles

### 1.1 Problem Statement

Given a real $N \times N$ matrix $X$, compute the principal matrix
logarithm $L = \ln(X)$ such that $\exp(L) = X$, with eigenvalues not
lying on the non-positive real axis (the branch cut).

### 1.2 Four-Step Spectral Method

**Step 1 — Eigendecomposition (DGEEV)**

$$
X \cdot V_R = V_R \cdot \Lambda
$$

$\Lambda$ contains eigenvalues (real or complex conjugate pairs);
$V_R$ is the real right eigenvector matrix.

**Step 2 — Invert eigenvectors (DGESV)**

Solve $V_R \cdot V_R^{-1} = I$ for $V_R^{-1}$.

**Step 3 — Build $D = \ln(\Lambda)$ as a real block-diagonal matrix**

- **Real eigenvalue** $\lambda_k > 0$: $D_{kk} = \ln(\lambda_k)$
  (if $\lambda_k \le 0$, error — branch cut violation)
- **Complex conjugate pair** $\lambda = a \pm ib = re^{\pm i\theta}$:
  $$
  D_{k:k+1,\,k:k+1} =
    \begin{bmatrix}
    \ln r & -\theta \\
    \theta & \ln r
    \end{bmatrix}
  $$

**Step 4 — Reconstruction (DGEMM × 2)**

$$
\ln X = V_R \cdot D \cdot V_R^{-1}
$$

The result is guaranteed real.

### Reference

Higham, N. J. *Functions of Matrices: Theory and Computation*,
SIAM, 2008, Ch. 11.

---

## 2. Code Implementation

### 2.1 Interface

```fortran
SUBROUTINE LOGEM(X, N, LNX, IERR)
IMPLICIT DOUBLE PRECISION (A-H, O-Z)
```

- **Input**: `X(N,*)` — real square matrix (column-major)
- **Output**: `LNX(N,*)` — $\ln(X)$
- **Output**: `IERR` — 0 = success, 1 = failure

### 2.2 Workspace Allocation

| Array | Size | Purpose |
|-------|------|---------|
| `A` | `N×N` | Copy of X (destroyed by DGEEV) |
| `WR, WI` | `N` | Eigenvalues (real/imag) |
| `VR` | `N×N` | Right eigenvectors |
| `VINV` | `N×N` | $V_R^{-1}$ |
| `D` | `N×N` | Block-diagonal $\ln(\Lambda)$ |
| `TEMP`, `RHS` | `N×N` | Workspace |
| `WORK` | $\max(1,8N)$ | LAPACK workspace |
| `IPIV` | `N` | Pivot indices |

### 2.3 Core Block — Constructing D

```fortran
K = 1
DO WHILE (K .LE. N)
  IF (ABS(WI(K)) .LT. TOL) THEN           ! real eigenvalue
    IF (WR(K) .LE. ZERO) THEN
      IERR = 1; RETURN
    END IF
    D(K,K) = LOG(WR(K))
    K = K + 1
  ELSE                                     ! complex pair
    A_VAL = WR(K)
    B_VAL = ABS(WI(K))
    R = SQRT(A_VAL*A_VAL + B_VAL*B_VAL)
    IF (R .LE. ZERO) THEN
      IERR = 1; RETURN
    END IF
    LR = LOG(R)
    THETA = ATAN2(B_VAL, A_VAL)
    D(K,   K  ) =  LR
    D(K,   K+1) = -THETA
    D(K+1, K  ) =  THETA
    D(K+1, K+1) =  LR
    K = K + 2
  END IF
END DO
```

Parameters: `TOL = 1.0D-14`, `ZERO = 0.0D+00`.

### 2.4 LAPACK Pipeline
