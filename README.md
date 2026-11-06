TESTING IMPLICATIONS OF THE NETWORK THEORY OF MENTAL DISORDERS IN INTENSIVE
SHORT-TERM TREATMENT
================================================================================

DESCRIPTION
-----------
This repository contains the analysis code for the manuscipt "Testing 
implications of the network theory of mental disorders in intensive short-term 
treatment" (Benz et al.).

The analyses are described in the manuscript and in more detail in the
Supplement.

Large parts of the code are adapted from the preliminary study by Schumacher
et al. (2023): https://osf.io/fhqmk/overview


CONTENTS
--------
The repository condains 2 R scripts, intended to be run in the order listed.

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


PREREGISTRATION
-----------------
This study is a preregistered secondary analysis. The preregistration is available at:
https://osf.io/bkngs/overview


DATA AVAILABILITY
-----------------
Patient data are not included in this repository. 


CITATION AND ATTRIBUTION
------------------------
Large parts of this code were adapted from code accompanying Schumacher et al. (2023).
Original code: https://osf.io/fhqmk/overview

Please retain attribution to the original code when adapting or redistributing
these materials.


DISCLAIMER
----------
Provided as-is for reproducibility purposes, without warranty of any kind.


LICENSE
-------
This code is licensed under the Creative Commons Attribution-
NonCommercial-ShareAlike 4.0 International License (CC BY-NC-SA 4.0).

https://creativecommons.org/licenses/by-nc-sa/4.0/


CONTACT
-------
Sabrina Benz — sabrina.benz@psy.lmu.de
