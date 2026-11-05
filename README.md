================================================================================
TESTING IMPLICATIONS OF THE NETWORK THEORY OF MENTAL DISORDERS IN INTENSIVE
SHORT-TERM TREATMENT
================================================================================

DESCRIPTION
-----------
The preparatory and main analyses are described in the manuscript "Testing 
implications of the network theory of mental disorders in intensive short-term 
treatment" (Benz et al.) and in more detail in the supplemental material of 
the manuscript.

Large parts of the code are adapted from the preliminary study by Schumacher
et al. (2023): https://osf.io/fhqmk/overview


CONTENTS
--------
Two R scripts, intended to be run in the order listed.

1) 0_preparatory_analyses.R
   Fits and compares six multilevel logistic model variants (random intercept
   only; + treatment interaction; random slopes for lagged symptoms, with and
   without treatment interaction; random slopes for the symptom-by-week
   interaction, with and without treatment interaction) for each of six
   symptoms. Includes model comparison used to select the final model (Model 1).

2) 1_main_analyses.R
   Refits Model 1 with 3000 post-warmup posterior draws. Contains:
     - Model estimation and convergence checks
     - Results inspection
     - Model fit evaluation
     - Network analysis
     - Sensitivity analysis: repeats the above with an added fixed main
       effect of week, testing robustness to unmodeled linear trends in
       overall symptom prevalence


DATA AVAILABILITY
-----------------
Patient data are not included in this repository. 


DISCLAIMER
----------
Provided as-is for reproducibility purposes, without warranty of any kind.


CONTACT
-------
Sabrina Benz — sabrina.benz@psy.lmu.de

================================================================================
