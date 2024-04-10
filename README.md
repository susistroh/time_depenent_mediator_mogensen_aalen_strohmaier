
Supplemental material to “Time-Dependent Mediators in Survival Analysis: Graphical Representation of Causal Assumptions” by Søren Wengel Mogensen, Odd Olai Aalen and Susanne Strohmaier. 

While data application in the main manuscript is based on the Systolic Blood Pressure Intervention Trial (SPRINT) data, we here want to illustrate the implementation in R using freely available data, the liver cirrhosis data attached to the timereg package. One potential downside of these data is that it might not comply so well with the causal assumption underlying the approach outlined in the theoretical sections of the paper. 

The liver cirrhosis data stem from a randomised trial conducted in the Copenhagen area between 1962 and 1969, including patients with histologically verified liver cirrhosis who were randomised to either receive prednisone or placebo at baseline. They were scheduled for regular visits after 3, 6 and 12 months and yearly thereafter and followed until death, end of study on October 1st, 1974 or loss to follow-up. For detailed description of the trial see for example Schlichting et al. 1983. 
Several versions of these data have previously been used to illustrate various methods for modelling time-to-event data with time-dependent covariates.  Andersen et al. (1993) use it at several occasions in their book, especially to illustrate multistate models. More recently, the data were used in Fosen et al. (2006) to exemplify dynamic path analysis using the additive hazards model. In the following, we address one of the questions also addressed in Fosen et al. namely, whether and to which extent the effect of prednisone on time to death is mediated through prothrombin.


Analysis of the effect of prednisone on death through prothrombin

The analysis set: 

For our illustration, we used information on 446 patients, 226 in the active treatment group and 220 in the placebo group, among whom 270 deaths (n=131 active group, n=139 placebo group) were observed. Following Fosen et al. (2006) we transformed the repeatedly measured prothrombin values by setting values of 70 or above to zero and subtracting 70 if the original values were below 70. As prothrombin values above 70 are considered ‘normal’, the transformed prothrombin values are zero for prothrombin values considered to be normal and negative for lower-than-normal values. 


The data set csl_CP is a modified version of the csl data included in the timereg package. 
It includes the following variables: 
Corresponding to the original variables as described here: https://rdrr.io/cran/timereg/man/csl.html

"id"     
"time"      
"prot"       
"dc"         
"eventT"     
"sex"        
"age"        
"treat"     
"prot.base"  
"prot.prev"  

The following variables are modified or newly created: 
“lt” and “rt” are the start and stop times of the relevant observation intervals 
"numbering"  - number of observation per individual 
"max_ob"    - maximum number of observation per individual 
"e_ind"  - event indicator (1=death)
"trans_prot" – transformed prothrombin value 

The analysis: 

We evaluated direct and indirect effects considering (i) the most recent value of prothrombin as well as (ii) an average of all previously measured values of prothrombin considering age and sex as confounding variables.

The R file code_04102023.R includes all relevant functions and function calls.

Required packages: 
library(survival)
library(timereg)
library(mvna)

Results: 

The estimated direct, indirect and total effects on a relative survival scale using the estimation procedure outlined in Section 6.1 of the main text are shown in Figure (file: overlay_04102023.pdf). The black lines display effect estimates when using the most recent value of prothrombin as mediating variable. The gray lines indicate the corresponding confidence intervals based on 1000 bootstrap samples. The blue lines display the effect estimates when considering the average of all previously measures of prothrombin as the mediating variable. 

While the total effect of prednisone seems fairly constant over the first 2.5 years, the slope indicates a beneficial effect thereafter. Based on the obtained direct and indirect effects it seems that this beneficial effect is mainly mediated through prothrombin. As indicated by the blue lines, it appears that incorporating information on the whole mediator history explains more of the total effect than simply using the most recent value (black lines). 


References: 

SCHLICHTING, POUL, CHRISTENSEN, ERIK, ANDERSEN, PER, FAUERHOLDT, LIS, JUHL, ERIK, POULSEN, HEMMING AND TYGSTRUP, NIELS. (1983). Prognostic factors in cirrhosis identified by Cox’s regression model. Hepatology 3(6), 889–895

ANDERSEN, PER KRAGH, BORGAN, ORNULF, GILL, RICHARD D AND KEIDING, NIELS. (1993). Statistical models based on counting processes. Springer. 

FOSEN, JOHAN, FERKINGSTAD, EGIL, BORGAN, ØRNULF AND AALEN, ODD O. (2006). Dynamic path analysis – a new approach to analyzing time-dependent covariates. Lifetime data analysis 12(2), 143–167






