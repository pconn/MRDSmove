# MRDSmove
Analysis of MRDS data subject to movement, measurement error, and (possible) observer dependence

This package contains files and code needed to conduct several analyses of double observer waterfowl survey data gathered by
Ray Alisauskas in northern Canada.  A simplistic approach using previously existing software tools for mark-recapture and 
mark-recapture distance sampling to examine factors affecting waterfowl detection is provided in the script /inst/analyze_waterfowl_mrds.R.
This file will recreate analyses in a manusript written for the Euring 2017 proceedings by Alisauskas and Conn.  However, the 
bulk of this packages consists of code to analyze MRDS models subject to movement and measurement error via maximum marginal likelihood, as
presented in a manuscript currently intended for JABES.  This includes a script to conduct such analyses on waterfowl data (see 
/inst/analyze_waterfowl_JABES.R), and scripts to conduct two simulation studies (/inst/simulation_study.R and /inst/simulation_study_pi.R).
The latter attempts to account for individual detection heterogeneity via an observer dependence specification that emulates point independence.
Functions to simulate data, to calculate the maximum marginal likelihood, to calculate Horvitz-Thompson-like abundance estimates, and to
conduct simulation-based goodness-of-fit testing are also availabe in the /R directory.  
