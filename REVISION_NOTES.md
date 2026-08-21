# Revision and correction notes

This log records substantive changes made while converting the final Markdown
sources to the Palatino LaTeX book and Quarto Live companion site. Routine
spelling, punctuation, formatting, and notation cleanup is not itemized.

## Structure and scope

- The main matrix-algebra chapter now retains only the linear-algebra tools used
  later in the course. Row-operation details, loops, inverse formulas,
  eigendecomposition, singular-value decomposition, and matrix derivatives are
  preserved in Appendix A.
- The main probability chapter now focuses on distributions, conditional
  expectations, random sampling, laws of large numbers, and the central limit
  theorem. Measure-theoretic examples, longer simulations, joint-normal details,
  and the delta method are preserved in Appendix B.
- A new Chapter 11 covers canonical DID, staggered adoption, group-time effects,
  event studies, weighting problems in conventional TWFE, aggregation,
  pre-trend interpretation, clustering, and modern R implementations.

## Mathematical and conceptual corrections

### Chapters 1-4

- Corrected claims equating full rank with invertibility for rectangular
  matrices and stated the conditions under which Gram matrices are invertible.
- Removed the false implication that every nonsymmetric matrix is
  diagonalizable; clarified rank and eigenvalue statements.
- Corrected the sign in the least-squares matrix gradient.
- Distinguished a random variable's support from the underlying sample space
  and qualified the probability-as-length heuristic.
- Corrected convergence definitions, the covariance in the multivariate CLT,
  and the multivariate delta-method formula.
- Corrected conditional unbiasedness to
  `E(beta_hat | X) = beta`, distinguished sample orthogonality from conditional
  mean independence, and repaired the test-statistic, confidence-interval, and
  p-value explanations.

### Chapters 5-7

- Reworked regression adjustment around centered, fully interacted Lin-style
  specifications; distinguished finite-sample bias from consistency and fixed
  the superpopulation variance comparison.
- Added overlap to CATE identification and clarified what unconfoundedness does
  and does not make observable.
- Corrected the local-estimator standard-error rate to `(n h^p)^(-1/2)` and
  qualified tree/forest unbiasedness and convergence claims.
- Corrected the AIPW conditional mean-zero argument, the DML nuisance-product
  rate, cross-fitting exposition, strong-overlap requirement, and efficiency
  claims.

### Chapters 8-11

- In RD, corrected one-sided limit directions, separated lack of overlap from
  degenerate unconfoundedness, and added warnings about bias and high-order
  global polynomials.
- In IV/LATE, separated random assignment from exclusion, corrected the
  covariance/relevance algebra and treatment-selection notation, and qualified
  the common first-stage `F > 10` rule.
- In the neural-network chapter, corrected repeated input indices, coefficient
  subscripts, projection-proof notation, and the discussion of DML inference.
- In the new DID chapter, kept the estimand explicit as `ATT(g,t)`, clarified
  already-treated comparisons and TWFE weights, separated pointwise intervals
  from simultaneous bands, and emphasized cohort composition, anticipation,
  clustering, and the limits of pre-trend tests.

## Figures and code

- The supplied `illuTree.png` and `NeuralNets.png` are retained with revised
  captions and placement.
- Because the original `Ch8RD_Lee.png` was unavailable, the RD illustration was
  reconstructed from `data/house.csv`. The reproducible script is
  `scripts/recreate_rd_figure.py`; both PDF and PNG outputs are included.
- Browser examples use WebR-compatible package sets and reduced classroom
  settings. The RD and DID pages pin WebR 0.5.9 for package compatibility.
- A native R 4.5 Codespace configuration is included for full forest fits and
  final bootstrap inference.
