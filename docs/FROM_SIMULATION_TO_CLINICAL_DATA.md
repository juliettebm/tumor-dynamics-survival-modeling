# From simulated data to clinical data

This portfolio project demonstrates modeling concepts under a known data-generating
mechanism. Applying the same ideas to patient-level clinical-trial or real-world data
would require a substantially stronger data, estimand, validation, and governance
workflow. A small p-value in this simulation is evidence that the implementation can
recover a programmed signal; it is not external clinical validation.

## 1. Define the question and estimand first

Before fitting a model, I would specify the population, treatment conditions, endpoint,
intercurrent events, censoring rules, and intended use of the result. PFS and OS would
be derived from a version-controlled specification rather than inferred ad hoc. The
analysis would distinguish prognostic association, treatment-effect prediction, and
surrogate-endpoint validation: these are different claims requiring different designs.

## 2. Build an auditable analysis dataset

I would reconcile patient, lesion, scan, dose, exposure, progression, treatment-change,
and death records. Automated checks would cover identifiers, chronology, units,
duplicates, impossible values, baseline definitions, scan windows, and consistency
between lesion-level and patient-level tumor burden. Every exclusion and correction
would be traceable to source data and a documented rule.

Missing tumor assessments and irregular visit schedules are informative possibilities,
not merely formatting problems. I would describe their patterns, compare included and
excluded patients, and use methods appropriate to the likely missingness mechanism.
Sensitivity analyses would address delayed or missed scans and alternative progression
dates.

## 3. Implement endpoints prospectively

RECIST 1.1 requires more than the simplified target-lesion rule implemented here. A
clinical pipeline would account for target and non-target lesions, new lesions,
confirmation rules where applicable, unequivocal progression, scan timing, and
independent review. PFS censoring would explicitly cover no post-baseline assessment,
new anticancer therapy, missed assessments, loss to follow-up, and database cut-off.
OS would use verified vital status and a prespecified censoring date.

Treatment discontinuation, switching, dose modification, and therapies received after
progression can affect exposure, PFS, and especially OS. Their handling would follow the
estimand and trial protocol; alternatives such as treatment-policy, hypothetical, or
while-on-treatment strategies would be explored when scientifically relevant.

## 4. Fit models that reflect the data structure

I would examine between-patient, between-lesion, and between-center heterogeneity and
include clinically justified covariates without relying on automated significance
screening alone. PK/PD models would use actual dosing histories and measured exposure,
with diagnostics for structural and residual-error assumptions. Longitudinal and
survival components would be checked for functional form, proportional hazards,
influential observations, and sensitivity to association structure.

Uncertainty from every estimation stage must propagate to downstream predictions. A
two-stage analysis would therefore be accompanied by bootstrap or multiple-imputation
analyses, and compared with an appropriate joint or nonlinear mixed-effects model.
Numerical integration would use a validated ODE solver with convergence and step-size
checks rather than relying only on fixed-step Euler integration.

## 5. Validate before making clinical claims

Internal validation would use patient-level resampling or cross-validation that keeps
all records for one patient together. Performance would be reported with confidence
intervals using discrimination, calibration, prediction error, and clinically relevant
decision metrics—not p-values alone. Model development choices would be frozen before
evaluation on an external cohort from another study, site network, or time period.

Transportability would be assessed across treatment arms, centers, demographic and
clinical subgroups, assay or imaging workflows, and follow-up patterns. Failure modes,
applicability limits, and recalibration rules would be documented. Only then could the
analysis support a clinical-development claim, and any deployment would require
monitoring, versioning, privacy controls, and independent review.

## 6. Minimum sensitivity-analysis set

- alternative PFS censoring and progression-date conventions;
- missing-at-random and plausible missing-not-at-random scenarios;
- alternative scan-window and landmark definitions;
- treatment discontinuation and post-progression therapy strategies;
- alternative tumor-dynamics and survival association structures;
- influential-patient, center, and subgroup analyses;
- bootstrap uncertainty and external temporal/geographic validation;
- ODE solver tolerances and alternative residual-error models.

The practical deliverable would be a reproducible analysis package containing the data
specification, derivation code, quality-control report, statistical analysis plan,
model diagnostics, validation report, sensitivity analyses, and a limitations section
that separates supported conclusions from exploratory findings.

The endpoint rules follow the principles in RECIST 1.1 and the treatment of
intercurrent events and sensitivity analyses should follow ICH E9(R1); full citations
are listed in [`REFERENCES.md`](../REFERENCES.md).
