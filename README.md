# Laurie et al. A Reevaluation of Risk Factors for Overwintering of *Podosphaera macularis* on Hop  

The data set is an assessment of hop powdery mildew flag shoots in commercial hop yards assessed during 2014 to 2020. Variables used in the analysis reported in Laurie et al (2023) are given here, as explained below. NA indicates missing data for a given variable. 
	
These files are provided to enable full reproducibility of the results presented in the paper. We make no warranties regarding these programs.

Corresponding author: David H. Gent dave.gent@usda.gov 

---

## Description of Data Set
---

**Data.csv**

This is the primary data set used in the analyses reported.  

---

**Column descriptions** 
-	Unique_ID = ID number unique to each location-year
-	Planting year = year of planting
-	Year = year of assessment for flag shoots
-	Yard_Age = age of yard in year of assessment 
-	Variety = variety name. 
	Note that in a few instances, cultivar name was not disclosed or numerous experiments lines were presents. In these instances, a generic name was 	  assigned.
-	Variety_Susceptibility = susceptibility rating of variety categorized using a 0 to 5 scale as described in paper
-	FS_Binary = binary measure of flag shoots; Absent = 0 Present = 1
-	Pruning_Type = pruning method
	Chemical = only chemical pruning
	Mechanic = only mechanical pruning (specifically crowning)
	Mowed = only mowed
	Mow_chem = combination of mowing and chemical pruning
	Notprune = no spring pruning was conducted
-	Mech_Prune = binary measure for mechanical pruning (1 = yes; 0 = no)
-	Pruning_Quality = quality rating of location-year assessed using a 1 to 5 ordinal scale as described in paper
-	Poor_Prune = binary measure of pruning quality as described in paper
	Thorough = 0 (pruning score <2)
	Incomplete = 1 (pruning score ≥2)
-	Suscept_Binary = binary measure of variety susceptibility as described in paper
	Susceptible = 0 (susceptibility score ≥2
	Low susceptibility = 1 (susceptibility score <2)
-	Flagshoots = number of flag shoots per plant
-	Second_Year = binary measure of yard age in year of flag shoot assessment
	Established = 0 
	Second year = 1
---

### R Code
**RWLGIT.R**
-	R code for all modelling and plots. Brief description of code chunks are given as comments in the R script 

---

**Rstan files supplemental to R code**
These are R.stan files associated with each of the Bayesian logistic regression models fit and presented in the paper.
-	mildew1.stan
-	mildew2.stan
-	mildew3.stan
-	mildew4.stan
-	mildew5.stan

---

**Supplemental data for plotting**
-	Flagshoots.csv
-	Yardage.csv

