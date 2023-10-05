Supplemental material to 
Time-dependent mediators in survival analysis: Graphical representation of causal assumptions
Morgensen, Aalen and Strohmaier 

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
lt and rt are the start and stop times of the relevant observation intervals 
"lt"         
"rt"        
"numbering"  - number of observation per individual 
"max_ob"    - maximum number of observation per individual 
"e_ind"  - event indicato (1=death)
"trans_prot" – transformed prothrombin value as described in the main text

code_04102023.R
As described in Section 7 we used the csl data to illustrate the estimation procedure in Section 6. The R file code_04102023.R includes all relevant functions and function calls.

Required packages: 
library(survival)
library(timereg)
library(mvna)
