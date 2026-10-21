# miR-146a/b complementary association analysis

Discovery run: from_scratch_20261017T061900.082
Fixed shared primary genes: CD3D, CD80, LCP2, SPI1
Exploratory analysis added after review; not part of original prespecification. No change to discovery set.
CPE is primary; LUMP is sensitivity. Eligible cancers are read from the final run, not selected by these outcomes.
Within each outcome/cancer/setting, both miRNAs, the outcome and purity must be finite. Identical complete cases are used for both conditional effects and reduced-model comparisons.
Variables are winsorised within that complete-case set at 0.1/99.9 percentiles (R quantile type 7), average-ranked, then standardised. This is a new rank-based conditional sensitivity, not an exact reproduction of the discovery residualisation.
Full model: standardised rank(outcome) ~ standardised rank(miR-146a-5p) + standardised rank(miR-146b-5p) + standardised rank(purity).
partial_rho is the correlation of outcome and focal-miRNA rank residuals after regressing both on the other miRNA and purity. It is distinct from beta_rank.
Inference: two-sided HC3 Wald z tests and pointwise normal-reference 95% intervals for beta_rank. These are approximate large-sample tests; no confidence interval for partial_rho is claimed.
Delta R2 compares the full model with the identical-sample reduced model omitting the focal miRNA. It is apparent in-sample fit improvement, not external predictive value.
BH and BY are separate by purity setting and outcome class (CYT vs shared genes), pooling BOTH miRNAs and ALL eligible cancer/outcome tests, with non-estimable planned tests retained in the denominator.
VIF >10 is flagged, not silently removed. Rank deficiency, constant purity or variables, leverage >=0.99, fewer than the configured minimum samples, or degenerate fits are non-estimable.
Conditional associations do not establish biological independence, causal regulation, direct targeting, or independence from immune-cell composition. Selected shared primary genes are not an unbiased transcriptome-wide validation sample.
Absence of significance does not prove redundancy; multicollinearity and power must be considered.

   setting outcome_type        mirna_id attempted estimated BH_supported
    <char>       <char>          <char>     <int>     <int>        <int>
1:     CPE          CYT hsa-miR-146a-5p        20        20           15
2:     CPE          CYT hsa-miR-146b-5p        20        20            8
3:     CPE  shared_gene hsa-miR-146a-5p        80        80           63
4:     CPE  shared_gene hsa-miR-146b-5p        80        80           40
5:    LUMP          CYT hsa-miR-146a-5p        19        19           12
6:    LUMP          CYT hsa-miR-146b-5p        19        19            8
7:    LUMP  shared_gene hsa-miR-146a-5p        76        76           63
8:    LUMP  shared_gene hsa-miR-146b-5p        76        76           42
   BY_supported VIF_over_10
          <int>       <int>
1:           13           0
2:            7           0
3:           54           0
4:           34           0
5:           12           0
6:            6           0
7:           51           0
8:           33           0
