********** 2005 **********************************************************  version 1.1  ****************;
*********************         SENSITIVITY ANALYSIS MACRO	              *******************************;
*********************************************************************************************************;													
* This macro will allow you to conduct a sensitivity analysis for a logistic regression analysis to     *;
* incorporate the effects of misclassification. By reconstructing your data set multiple times within   *;
* the range of parameters specified, the final analysis will produce a conventional logistic regression *;
* output (only random error), a sensitivity analysis output (only systematic error), and simulation     *;
* interval (both systematic and random error).                                                          *; 
*********************************************************************************************************;
* Because this involves multiple imputations of your data, the program may take long time to complete.	*;												 							*;
*********************************************************************************************************;
* To run the macro, change the code in this program. Before running, place the data set, and the       	*;
* macro program in the same directory. Parameters with an O next to them are optional.          		*;
*********************************************************************************************************;

/*%sensmac(
************************ MACRO SETTINGS *****************************************************************;
O		startover= 	Delete old repetitions? yes = delete , no = append to old repetitions
	(NOTE: IF MISSING WILL BE SET TO yes);
O		outset=		Name of the data set to keep multiple iterations in the libname specified below
	(NOTE: IF MISSING WILL BE SET TO "tempdata");
O		totalreps= 	Number of times to repeat the whole macro 
			0 = determine convergence before stopping or reach max iterations (3000)
	(NOTE: IF MISSING WILL BE SET TO 0);				
O		log= 		Regulates how much log to produce, 0 = full log, 1 = limited log, 2 = no log
	(NOTE: IF MISSING WILL BE SET TO 2);

************************ DATA SET AND VARIABLES ********************************************************;
		libname=	Name of library with the programs and where to store files - no quotes
		inset= 		Input data set containing all the variables
	(NOTE: Should be in the work library)
	(NOTE: ALL variables below must be either dummy variables or continuous variables)
		depend= 	Dependent variable of interest
	(NOTE: Should be a 1/0 dummy variable with 1 as the outcome of interest)
		exp= 		The exposure of interest
	(NOTE: Should be a 1/0 dummy variable with 1 as the outcome of interest)
O		indep=		String of independent predictors FOR FINAL MODEL

************************ MISCLASSIFICATION SETTINGS ***************************************************;
O		miscvar=	Variable that has been misclassified 
	(NOTE THIS MUST BE IN THE indep string OR IT MUST BE EITHER THE depend OR exp)
	(NOTE IF MISSING WILL BE SET TO &exp)
O		mis_ind=	Dichotomous variable in which prevalence of misclassified variable is expected to vary 
	(NOTE: IF MISSING WILL BE SET TO &depend)		
		sens_min=   Sensitivity of misclassified value minimum
O 		sens_mod=   Sensitivity of misclassified value mode1
O		sens_mod2=  Sensitivity of misclassified value mode2
		sens_max=	Sensitivity of misclassified value maximum
		spec_min=	Specificity of misclassified value minimum
O 		spec_mod=   Specificity of misclassified value mode1
O		spec_mod2=  Specificity of misclassified value mode2
		spec_max=	Specificity of misclassified value maximum
	 (NOTE: IF PARAMETERS BELOW ARE MISSING, MISCLASSIFICATION IS ASSUMED NON-DIFFERENTIAL)
O		sens_min_d= Sensitivity of misclassified value minimum (among mis_ind=0)
O 		sens_mod_d= Sensitivity of misclassified value mode1 (among mis_ind=0)
O		sens_mod_d2=Sensitivity of misclassified value mode2 (among mis_ind=0)
O		sens_max_d=	Sensitivity of misclassified value maximum (among mis_ind=0)
O		corrsens= 	Correlation between sensitivities if differential
	(NOTE: IF MISSING AND DIFFERENTIAL WILL BE SET TO 0.8)
0		spec_min_d=	Specificity of misclassified value minimum (among mis_ind=0)
O 		sens_mod_d= Specificity of misclassified value mode1 (among mis_ind=0)
O		sens_mod_d2=Specificity of misclassified value mode2 (among mis_ind=0)
O		spec_max_d=	Specificity of misclassified value maximum (among mis_ind=0)
0		corrspec= 	Correlation between specificities if differential 
	(NOTE: IF MISSING AND DIFFERENTIAL WILL BE SET TO 0.8)
);*/

* BELOW ARE EXAMPLES OF ANALYSES THAT COULD BE RUN ON A DATA SET CALLED 'EXAMPLE'
* TO RUN, SAVE THE DATASET AND THE FILES "Sensitivity.sas" AND Run "Sensitivity.sas"  
* IN A FOLDER C:\Sensitivity WHICH YOU WILL NEED TO CREATE ON YOUR C DRIVE;

libname sensdir 'C:\Sensitivity\';

 
* Create data set in work directory;                                                                             
data example; set sensdir.example;                      
run; 
title;

%include 'C:\Sensitivity\Sensitivity Analysis.sas';

*** NON-DIFFERENTIAL MISCLASSIFICATION;
%sensmac(                                                                               
    libname=    C:\Sensitivity\,                                                            
	inset=      example,
    depend=     case,
	outset=		exampleND,
	startover=	yes, 
	log=		2,
    exp=      	exp,
	totalreps=	5000,
    sens_min=   .75, 
	sens_mod=   .85,  
	sens_mod2=  .95,   
    sens_max=   1,                                                                 
    spec_min=   .75,                                                                
	spec_mod=   .85,  
	spec_mod2=  .95,
	spec_max=   1
);    

*** DIFFERENTIAL MISCLASSIFICATION;
%sensmac(                                                                               
    libname=    C:\Sensitivity\,                                                            
	inset=      example,
    depend=     case,
	outset=		exampleD,
	startover=	no, 
	log=		2,
    exp=      	exp,
	totalreps=	5000,
    sens_min=   .75,
	sens_mod=	.85,
	sens_mod2=	.95, 
    sens_max=   1,                                                                                                                               
    spec_min=   .75,
	spec_mod=	.85,
	spec_mod2=	.95,  
    spec_max=   1,                                                                
    sens_min_d= .70,
	sens_mod_d=	.80,
	sens_mod_d2=.90,  
    sens_max_d= .95,                                                                                                                                
    spec_min_d= .70,
	spec_mod_d=	.80,
	spec_mod_d2=.90, 
    spec_max_d= .95,
	corrsens=	.8,
	corrspec=	.8
);   


