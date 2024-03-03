# SNAP2: Critical Care Unit bed availability and postoperative outcome
R &amp; Stata code for conducting the analyses contained within the SNAP 2 study: Critical Care Unit bed availability and postoperative outcome – a multinational cohort study.   
The analyses were run using R Version 4.0.3 in RStudio Version 1.2.5033 & Stata Release 17.    

The order for running scripts is as follows: 
1. SNAP2_IV_Main_Analyses_Github.Rmd
2. SNAP2_IV_Sensitivity_Analyses_Github.Rmd
      
Stata scripts need to be run alongside these R markdown files - the relevant scripts to be run are indicated at the relevant points in the markdown files. This dependency is because some of the dataframes used in the Stata scripts are created during the course of the markdown files, and some of the results generated by the Stata scripts are used in the markdown files (e.g. for plotting of results).    

The required datasets to work with this code are subject to reasonable request from the corresponding author but will be subject to the approval of the data controller (The Royal College of Anaesthetists) and the sponsor (University College London) for the request. The data consists of an anonymised dataset with all variables necessary for the replication of the study, including occupancy data.    

If there are any issues with the code, please let us know by raising an issue and we shall get back to you as quickly as possible. 
