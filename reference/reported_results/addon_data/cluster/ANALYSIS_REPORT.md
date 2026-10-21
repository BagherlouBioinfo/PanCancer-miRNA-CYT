# Frozen-cluster incremental survival analysis

Discovery run: from_scratch_20261017T061900.082
Exploratory analysis added after peer review. Frozen labels; no reclustering, feature reselection, outcome cutpoints, or external prognostic validation.
Within each endpoint/cancer/tier, baseline and full models use identical complete cases. Age, CYT and ordinal stage are standardised within that set; constant nuisance variables are omitted.
Primary endpoint: OS. DSS and PFI are exploratory endpoint sensitivities, not universally recommended clinical endpoints across all cancers.
Primary comparisons: CYT + age + sex versus the same model plus fixed cluster; stage-adjusted comparisons use saved Stage 15 eligibility.
BH correction is separate by endpoint and tier across ALL attempted cancer comparisons, retaining non-estimable tests in the denominator.
Global and term-specific Schoenfeld PH tests are diagnostic, nominal and unadjusted. A non-significant test does not prove PH.
For every estimable OS comparison, a fixed time-varying sensitivity is attempted: all retained nuisance covariates have main and x*log((t+1)/366) terms in BOTH models; the full model adds cluster and tt(cluster), with a joint 2-df LRT.
The time-varying family has a separate BH correction per tier. EPP is applied to the doubled parameter count. This model is a sensitivity analysis, not a guaranteed correction of all non-PH forms. cox.zph is not applied to tt models.
HR(t) is reported at 1, 3 and 5 years only where at least 10 participants remain at risk. These are pointwise, not simultaneous confidence intervals.
Changes in concordance are in-sample apparent values, not validated prediction gains. Fixed cluster contrasts preserve original numeric label order.
Warning-generating, singular, non-finite or convergence-caution fits are excluded from inference and reported with reasons.

   endpoint              tier attempted estimated BH_significant
     <char>            <char>     <int>     <int>          <int>
1:       OS       CYT_age_sex        31        19              2
2:       OS CYT_age_sex_stage        31        15              0
3:      DSS       CYT_age_sex        31        18              2
4:      DSS CYT_age_sex_stage        31        13              0
5:      PFI       CYT_age_sex        31        23              1
6:      PFI CYT_age_sex_stage        31        17              0
   PH_flagged TV_estimated TV_BH_significant
        <int>        <int>             <int>
1:          4           15                 2
2:          6           10                 1
3:          3            0                 0
4:          4            0                 0
5:          3            0                 0
6:          4            0                 0
