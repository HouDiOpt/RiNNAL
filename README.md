# RiNNAL

**RiNNAL version 1.0 — a MATLAB software for doubly nonnegative relaxations of mixed-binary quadratic programming problems**

### [Di Hou](https://houdiopt.github.io), [Tianyun Tang](https://ttymath.github.io/tianyuntang.github.io/), [Kim-Chuan Toh](https://blog.nus.edu.sg/mattohkc/)

This software is designed to solve doubly nonnegative relaxations of mixed-binary quadratic programming.

**Mixed-Binary Quadratic Programming (MBQP):**

```tex
minimize   x'Qx + 2c'x  
subject to Ax = b                                (equality)  
           x_i ∈ {0,1}, for i ∈ B                (binary)  
           x_i·x_j = 0, for (i, j) ∈ E           (complementarity)  
           x ≥ 0                                 (nonnegative)  
```

**Doubly Nonnegative (DNN) Relaxation:**

```tex
minimize   ⟨Q, X⟩ + 2c'x  
subject to Ax = b  
           A X = b x'  
           diag(X)_B = x_B  
           X_E = 0  
           X ≥ 0   (nonnegative)
           X ⪰ 0   (positive semidefinite)  
```

------

## Citation

- Di Hou, Tianyun Tang, Kim-Chuan Toh, **A low-rank augmented Lagrangian method for doubly nonnegative relaxations of mixed-binary quadratic programs,** [arXiv:2502.13849](https://arxiv.org/abs/2502.13849)

```bibtex
@article{hou2025low,
  title={A low-rank augmented Lagrangian method for doubly nonnegative relaxations of mixed-binary quadratic programs},
  author={Hou, Di and Tang, Tianyun and Toh, Kim-Chuan},
  journal={arXiv preprint arXiv:2502.13849},
  year={2025}
}
```

------

**Important note**:  

- The software is still under development. Thus it will invariably be buggy. We would appreciate your feedback and bugs’ report.
- This is a research software. It is not intended nor designed to be a general purpose software at the moment.

## Examples

1. run the script `run_BIQ_demo.m` in Matlab.
2. run the script `run_QKP_demo.m` in Matlab.
