/*********************************************************************************************************
  Statistical analysis of the type 2 MAD
  Copyright 2013 by Dr. Frank You. All rights reserved
  
  Program: MDAPipeline_Step1a.sas
  Author:  Dr. Frank M. You
  Affiliation: Agriculture and Agri-Food Canada, Winnipeg, 
             Dafoe Road 195, R3T3N9, Manitoba, Canada
  
  Description: Perform ANOVA for plot control and subplot controls

  Input:

  Output: two files are exported for the downstream analysis
    	(1) unadjusted_subplot_anova_stats.txt
	    (2) unadjusted_plot_anova_stats.txt  

  Reference: 
         Please read the You et al. (2013) paper:

		Frank M. You1, Scott D. Duguid, Dinushika Thambugala and Sylvie Cloutier.(2013) Statistical analysis 
        and field evaluation of the type 2 modified augmented design in phenotyping of flax germplasms in 
        multiple environments. Australia J. Crop Science. 
**********************************************************************************************************/


/*********************************************************************************************************
   Create a library called MAD. All output data set will saved in this library 
**********************************************************************************************************/

libname MAD '.';



/*ODS TAGSETS.EXCELXP
file='C:\Users\rwm132\flax\madpipeline\sasoutput\BxN_Preston_Wilt_Data_2019_for_MAD_converted.xls'
STYLE=minimal
OPTIONS ( Orientation = 'landscape'
FitToPage = 'yes'
Pages_FitWidth = '1'
Pages_FitHeight = '100' );*/


/*ODS TAGSETS.EXCELXP
file='C:\Users\rwm132\flax\madpipeline\sasoutput\BxN_Preston_Wilt_Data_2019_for_MAD_converted.xls'
STYLE=Printer
OPTIONS ( Orientation = 'landscape'
FitToPage = 'yes'
Pages_FitWidth = '1'
Pages_FitHeight = '100' 
embedded_titles = 'yes');*/



/*ods excel file="C:\Users\rwm132\flax\madpipeline\sasoutput\BxN_Preston_Wilt_Data_2019_for_MAD_converted.xlsx" ;*/


ods pdf file="C:\Users\rwm132\flax\madpipeline\sasoutput\BxN_Preston_Wilt_Data_2019_for_MAD_converted..pdf";



/*********************************************************************************************************
   Input the phenotypic data of the type 2 MAD 
**********************************************************************************************************/

data MAD.mad_phenotypic_data;
	infile 'C:\Users\rwm132\flax\madpipeline\BxN_Preston_Wilt_Data_2019_for_MAD.txt_converted.txt' dlm='09'x LRECL=2400 firstobs=2;
	
	length Genotype $40  Whole_plot $10;
	input Record Plot Row Column Cp Csp Entry Year $ Location $ Genotype $ Trait $ Value;	

	Whole_plot = '.';
    if Csp > 0 or Cp = 1 then Whole_plot = cat (Row, '+', Column); 
	if _n_ > 1;
run;

proc print data=MAD.mad_phenotypic_data;
run;

/*********************************************************************************************************
 Preparing control plot and subplot control data for ANOVA of individual experiments
**********************************************************************************************************/

/**** plot and subplot control data ****/
proc sql;
	create table MAD.MDA_subplot_controls as
		select *
    	from MAD.mad_phenotypic_data
		where Csp > 0
		order by Year, Location, Trait;
quit;

/* put subplot controls and the plot control in the whole plot with subplot controls together 
  Only plot control data which have subplot controls. This data set is for ANOVA of subplot controls.
  In the paper of Lin and Voldeng (1989), it is suggested that the plot control data in the several 
  whole plots with subplot controls is also used for ANOVA of subplot controls to increase test power.
*/ 
proc sql;
	create table MAD.MDA_plot_subplot_controls as
		select distinct t1.* 
    	from MAD.mad_phenotypic_data as t1, MAD.MDA_subplot_controls as t2
		where t1.Whole_plot = t2.Whole_plot and t1.Year = t2.Year 
			and t1.Location = t2.Location and t1.Trait = t2.Trait and t1.Whole_plot ne '.'
		order by Year, Location, Trait, Whole_plot;
RUN;
quit;



/**** plot control data ****/
proc sql;
	create table MAD.MDA_plot_controls as
		select *
    	from MAD.mad_phenotypic_data
		where Cp = 1 and Csp = 0
		order by Year, Location, Trait;
quit;

/* ANOVA of plot controls and subplot controls for any number of traits
   The ANOVA statistics is exported to the data set MAD.plot_control_anova_stat
   and MAD.subplot_control_anova_stat
*/ 
PROC GLM DATA=MAD.MDA_plot_controls outstat=MAD.plot_control_anova_stat;
	BY Year Location Trait;
	CLASS  Row Column;
	MODEL Value = Row Column/SS3;
	RANDOM Row Column;
	LSMEANS Row Column; 
RUN;
quit;


PROC GLM DATA=MAD.MDA_plot_subplot_controls outstat=MAD.subplot_control_anova_stat;
	BY Year Location Trait;
	CLASS  Whole_plot Genotype;
	MODEL Value = Whole_plot Genotype/SS3;
	RANDOM Whole_plot; 
	LSMEANS Whole_plot;
RUN;
QUIT;


/*********************************************************************************************************
   Output results tables  to text files for the downstream analysis
**********************************************************************************************************/

/* Subplot control ANOVA results */
data _null_ ;          
    set MAD.subplot_control_anova_stat;
    FILE  'C:\Users\rwm132\flax\madpipeline\sasoutput\unadjusted_subplot_anova_stats.txt' DLM='09'x ;  
    PUT  Year Location Trait _SOURCE_ DF SS F PROB;
run ; 

/* Whole plot control ANOVA results */
data _null_ ;          
    set MAD.plot_control_anova_stat;
    FILE  'C:\Users\rwm132\flax\madpipeline\sasoutput\unadjusted_plot_anova_stats.txt' DLM='09'x ;  
    PUT  Year Location Trait _SOURCE_ DF SS F PROB;
run ; 

/************************************ END OF ANALYSIS ****************************************************/

/*ods excel close;*/
ods pdf close;
