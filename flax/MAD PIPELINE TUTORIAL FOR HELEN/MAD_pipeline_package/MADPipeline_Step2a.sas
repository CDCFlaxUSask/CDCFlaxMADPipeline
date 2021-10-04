/*********************************************************************************************************
  Statistical analysis of the type 2 MAD
  Copyright 2013 by Dr. Frank You. All rights reserved

  Program: MDAPipeline_Step2a.sas
  Author:  Dr. Frank M. You
  Affiliation: Agriculture and Agri-Food Canada, Winnipeg, 
             Dafoe Road 195, R3T3N9, Manitoba, Canada
  
  Description: Perform joint ANOVA of MAD data if the same design with the same plot and subplot controls
               are used in all experiments

  Input: Adjusted values obtained from the program MAD_analysis_Step1b.pl
		MAD_adjusted_all_data_suggested.txt

  Output: If you have data from multiple years and locations, you can perform three types of Joint ANOVA. 
        Thus you can get output files for joint ANOVA for all years and location, by years, and by locations:
    	(1) multi_envirom_anova_stat_YL.txt: results for joint ANOVA for all years and locations
	    (2) multi_envirom_anova_stat_by_year.txt: results for joint ANOVA by years.  
	    (3) multi_envirom_anova_stat_by_loc.txt: results for joint ANOVA by locations.  
		These three files can be further summarized by the program "MAD_analysis_Step2b.pl to determine
        correct F tests if a mixed model used. Since the Proc GLM treats all effects fixed, you need to 
        correctly determine denominaters in F tests. 

  reference:
        Please see You et al. (2013) paper:

		Frank M. You1, Scott D. Duguid, Dinushika Thambugala and Sylvie Cloutier.(2013) Statistical analysis 
        and field evaluation of the type 2 modified augmented design (MAD) in phenotyping of flax (Linum usitatissimum) 
        germplasms in multiple environments. Australia J. Crop Science. 
**********************************************************************************************************/


/*********************************************************************************************************
   Create a library called MAD. All output data set will saved in this library 
**********************************************************************************************************/

libname MAD '.';

/* Step 1: input the adjusted phenotypic data */
data MAD.mad_phenotypic_data_adjusted;
	infile 'D:\loongf\US\helen\raja\temp\MAD_adjusted_all_data_suggested.txt' dlm='09'x LRECL=1000;

	length Genotype $30;
	input Record Plot Row Column Cp Csp Entry Year $ Location $ Genotype $ Trait $ Value;	

 	if _n_ > 1;
run;

/* Step 2a: Joint ANOVA for multiple location by year */
PROC SORT DATA=MAD.mad_phenotypic_data_adjusted;
	BY Year Trait;
RUN;

PROC GLM DATA=MAD.mad_phenotypic_data_adjusted outstat=MAD.multi_envirom_anova_stat_by_year;
	by Year Trait;
	class  Location Genotype;
	model Value = Location|Genotype /SS3;
run;

/* Step 2b: Joint ANOVA for multiple years by location */
PROC SORT DATA=MAD.mad_phenotypic_data_adjusted;
	BY Location Trait;
RUN;

PROC GLM DATA=MAD.mad_phenotypic_data_adjusted outstat=MAD.multi_envirom_anova_stat_by_loc;
	by Location Trait;
	class  Year Genotype;
	model Value = Year|Genotype /SS3;
run;

/* Step 2c: Joint analysis for  all years and locations*/
PROC SORT DATA=MAD.mad_phenotypic_data_adjusted;
	BY Trait;
RUN;

PROC GLM DATA=MAD.mad_phenotypic_data_adjusted outstat=MAD.multi_envirom_anova_stat_YL;
	by Trait;
	class  Year Location Genotype;
	model Value = Year|Location|Genotype /SS3;
	random Year; 
run;


/* output results tables */
/* ANOVA results of data from multiple years and locations  */
data _null_ ;          
    set MAD.multi_envirom_anova_stat_YL;
    FILE  'D:\loongf\US\helen\raja\temp\multi_envirom_anova_stat_YL.txt' DLM='09'x ;  
    PUT  Trait _SOURCE_ _TYPE_ DF SS F PROB;
run ; 

/* ANOVA results of data from multiple years  */
data _null_ ;          
    set MAD.multi_envirom_anova_stat_by_year;
    FILE  'D:\loongf\US\helen\raja\temp\multi_envirom_anova_stat_by_year.txt' DLM='09'x ;  
    PUT  Year Trait _SOURCE_ _TYPE_ DF SS F PROB;
run ; 

/* ANOVA results of data from multiple locations  */
data _null_ ;          
    set MAD.multi_envirom_anova_stat_by_loc;
    FILE  'D:\loongf\US\helen\raja\temp\multi_envirom_anova_stat_by_loc.txt' DLM='09'x ;  
    PUT  Location Trait _SOURCE_ _TYPE_ DF SS F PROB;
run ; 


/************************************ END OF ANALYSIS ****************************************/


