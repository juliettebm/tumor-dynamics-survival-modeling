# Individual Tumor Dynamics and Their Link to Survival (Simulated Data)

## Context

This project reproduces, on **fully simulated data**, the methodological approach used
in pharmacometrics to link individual tumor dynamics (longitudinal RECIST-type
measurements) to survival endpoints (PFS/OS) in oncology, a central question in
clinical development of anticancer therapies (e.g., Claret et al., Wang et al.,
*Clin Pharmacol Ther*, 2009).

No real patient data is used. The goal is to demonstrate methodological understanding,
simulation, individual parameter estimation, survival analysis, and joint modeling,
before applying the same approach to real-world data.

## Notebooks

- **`01_simulation_and_two_stage_model.Rmd`**: simulates 200 patients with individual
  tumor growth/regrowth parameters and survival outcomes, then estimates these
  parameters back from noisy longitudinal data using a classic "two-stage" approach
  (individual nonlinear fits followed by a Cox model). Includes a first attempt with a
  sparse observation design (weak, non-significant signal) and a second, improved design
  (significant signal), illustrating how observation quality affects signal detection.
- **`02_joint_model.Rmd`**: re-analyzes the same (improved-design) data using a joint
  model (`JM` package), estimating the longitudinal trajectory and survival
  simultaneously, and compares the results honestly against the two-stage approach.

## Model

Tumor growth inhibition (TGI) biexponential model (Wang et al., 2009):

```
SLD(t) = SLD0 * (exp(-g*t) + exp(d*t) - 1)
```

- **g**: tumor regression rate (treatment effect)
- **d**: tumor regrowth rate (resistance)
- **SLD0**: baseline tumor size (Sum of Longest Diameters)

Survival is simulated via a Cox proportional hazards structure in which `g` is
protective and `d` is detrimental.

## Key results

| Approach | Design | cor(g, g_hat) | cor(d, d_hat) | Signal on survival |
|---|---|---|---|---|
| Two-stage | Sparse (6-week scans, 15% noise) | 0.568 | 0.919 | Not significant (p = 0.60) |
| Two-stage | Improved (3-week scans, 8% noise) | 0.830 | 0.990 | Significant (p = 0.007) |
| Joint model | Improved design, simplified linear trend | - | - | Borderline (p = 0.056) |

**Main takeaway.** Detecting a true tumor-dynamics/survival association is highly
sensitive to the observation design (scan frequency, measurement precision),a real
but under-estimated design effect ("attenuation by measurement error"). The joint model
provides a more statistically rigorous framework than the two-stage approach, at the
cost of a simplified representation of tumor dynamics, illustrating a genuine
methodological trade-off between model realism and statistical tractability.

## Tools

R — packages: `minpack.lm`, `survival`, `survminer`, `nlme`, `JM`, `dplyr`, `ggplot2`

## Author

Juliette Bouli-Mengue 

