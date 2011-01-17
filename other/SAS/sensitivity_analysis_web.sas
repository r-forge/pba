*set SAS page and line options; 
		options ps=40 ls=124;
*assign two formats;
		proc format;
			value 	t2x		1 = '1:Yes'
							0 = '0:No';
			value definf	1 = 'A:Less than Definitive'
							0 = 'B:Definitive';
*set the library name containing the data files;
*change the directory name in quotes to the directory containing the SAS data files;
		libname sensdir v6 'y:\sph\epidemiology\tlash\pubs\sens anal fink';
*set the SAS graphing options for the graphs that generate the cumulative distributions;
*all fonts are set to Swiss;
*legend option assumes three lines on the graph;
		legend1 across=3 frame fwidth=4 label=none mode=reserve
			position=(outside bottom center) value=(F=swiss H=9 pt);
*options for the y-axis;
		axis1 c=bl label=(A=90 F=swiss H=16 pt J=center 'Cumulative percentile') 
			major=(W=2) offset=(0,0) width=4
			order = (0 to 110 by 10) value=(F=swiss H=12 pt '0' '10%' '20%' '30%' '40%'
					'50%' '60%' '70%' '80%' '90%' '100%' ' ');
*options for the x-axis, on the arithmetic scale;
		axis2 c=bl label=(A=0 F=swiss H=16 pt J=center 'ln(Relative hazard)') 
			value=(F=swiss H=12 pt) major=(W=2) offset=(0,1) width=4;
*options for the x-axis, on the logarithmic scale;
		axis3 c=bl label=(A=0 F=swiss H=16 pt J=center 'Relative hazard') 
			value=(F=swiss H=12 pt) offset=(0,1) width=4 logbase=10 logstyle=expand;
*output of proc contents for SAS dataset senssamp;
/*
	#    Variable    Type    Len    Pos    Format     Label
	-----------------------------------------------------------------------------
	4    ADJCHEM     Num       8     24    T2X.       Adjuvant chemotherapy
	13   AGECAT1     Num       8     96    T2X.       65 to 74
	14   AGECAT2     Num       8    104    T2X.       75 to 90
	10   ALLCOMOR    Num       8     72               Sum of comorbidities
	7    AXIL        Num       8     48    T2X.       Axillary Dissection
	9    BCCAUSE     Num       8     64    T2X.       Death from breast cancer
	2    BCSSURG     Num       8      8    T2X.       Breast conserving surgery
	6    DEFNTHER    Num       8     40    DEFINF.    Definitive therapy
	12   EXCAT1      Num       8     88    T2X.       Regional
	5    GRPHOS      Num       8    112    T2X.       Hospital w/o tumor registry
	3    MASTSURG    Num       8     16    T2X.       Mastectomy
	11   MORTIME2    Num       8     80               Time to death (years)
	1    REID        Num       8      0    T2X.       Reidentified
	5    STOPPED     Num       8     32    T2X.       Chemo terminated, non-death
	8    TTORAD      Num       8     56               Time to radiation therapy
*/
*selection bias section;
*create a temporary data set with records for those who were reidentified (reid=1);
		data yes_reid1; set sensdir.senssamp;
			if reid eq 1;
		run;
*create a temporary data set with records for those who were not reidentified (reid=0);
		data no_reid1; set sensdir.senssamp;
			if reid eq 0;
		run;	
*use logistic regression to estimate risk of breast cancer mortality among those who were reidentified
 (data=yes_reid1) based on their characteristics;
*put the parameter estimates and the covariance matrix into a data set called logitout;
*do not print the results of the logistic regression to the output window;
        proc logistic data=yes_reid1 descending outest=logitout covout noprint;
                model bccause = excat1 agecat1 agecat2 allcomor defnther grphos;
		run;
*invoke SAS interactive matrix language to get the parameter estimates (central tendency plus a normal deviate
times the standard error, accounting also for covariance between parameters).  These parameter estimates will
be used to estimate the risk of breast cancer mortality for women who were not reidentified';
		proc iml;
*convert the parameter and covariance matrix (logitout) to a matrix called cova that contains
 only the covariance matrix (_type_=cov);
			use logitout;
			read all var {intercept excat1 agecat1 agecat2 allcomor defnther grphos} 
				into cova where(_type_='COV');
*convert the parameter and covariance matrix (logitout) to a row matrix called parama that contains
 only the parameters (_type_=parms);
			read all var {intercept excat1 agecat1 agecat2 allcomor defnther grphos} 
				into parama where(_type_='PARMS');
*transpose the parameter matrix to a column matrix;
			parama = t(parama);
*initialize a row matrix of standard deviates (stda);
			stda = {0,0,0,0,0,0,0};
*calculate the Cholesky root of the covariance matrix (covroot);
			covroot = root(cova);
*transpose the Cholesky root matrix (covrtt);
			covrtt = t(covroot);
*for the number of iterations (2000), calculate the parameter estimates as the parameter matrix plus the 
 product of the transposed covariance root and standard normal deviate;
*change the number of iterations (the number of iterations must be constant throughout);
*the optimal number of iterations depends on the size of the data set (# of records and # of variables);
*on a typical pc, start by dividing 1E6 by the number of records to estimate the number of iterations;
			do icount = 1 to 2000;
*fill the standard normal deviates matrix with standard normal deviates;
				stda = normal(stda);
*calculate the parameters column matrix (pouta);
				pouta = parama+covrtt*stda;
*transpose the column matrix to a row matrix (poutb);
				poutb = t(pouta);
*create a SAS data set (param) on the first iteration and write the row matrix to it;
				if icount = 1 then create param from poutb;
*if not the first iteration, then append the row matrix to the param data set;
				append from poutb;
			end;
*quit IML;
			quit;
*create a SAS data set (param2) from the output matrix, renaming the matrix columns with variable names;
*index the data set with a variable called repcnt, an integer counting from 1 to the number of iterations;
		data param2 (rename=(col1=interp col2=excat1p col3=agecat1p col4=agecat2p 
					col5=allcomp col6=deftherp col7=grphosp) index=(repcnt));
*for each iteration in the parameter data set, set the repcnt variable equal to the iteration;
				do obsnum=1 to last;
				set param point=obsnum nobs=last;
				repcnt=obsnum;
				output;
				end;
				stop;
		run;
*use the parameter data set, in combination with individual characteristics, to calculate the risk of breast cancer 
 mortality for each woman who was not reidentified.  Use the risk to calculate whether or not she would have died of
 breast cancer in this iteration;
*the data set no_reid2 will contain the breast cancer outcomes for the women who were not reidentified.  There will
 be the same number of groups of outcomes as there are iterations of the sensitivity analysis;
        data no_reid2;
*seed2 is the random number generator seed, so is set to -1 on the first iteration and retained 
 for subsequent iterations;
				retain seed2;
				if seed2 eq . then seed2 = -1;
*change the number of iterations (reptot);
				repcnt = 1; reptot = 2000;
*for each iteration, select all of the records in the data set of non-reidentified women, and
 for each record, append the corresponding parameter values;
*the first do loop counts iterations;
				do while (repcnt le reptot);
*the second do loop counts records;
					do reccnt = 1 to last;
*the first set statement reads the record from the data set of non-reidentified women;
					set no_reid1 point=reccnt nobs=last;
*the second set statement reads the parameter record for the iteration;
					set param2 key=repcnt;
*use the characteristics in the record and the parameter estimates for the iteration
 to calculate the logit of risk of breast cancer mortality;
                logite = interp+agecat1*agecat1p+agecat2*agecat2p 
                         +excat1*excat1p+allcomor*allcomp+grphos*grphosp;
*transform the logit of risk to risk;
                pprob = exp(logite) / (1+exp(logite));
*use the risk and random binary generator to assign breast cancer mortality;
*0 and 1 logic checks of the bounds are necessary to avoid errors when calling ranbin;
				if pprob eq 0 then bccause = 0;
					else if pprob eq 1 then bccause = 1;
					else call ranbin(seed2,1,pprob,bccause);
*uncomment the following statement to disable the selection bias component of the sensitivity analysis;
*				bccause = .;
*for non-reidentified women who were not assigned death due to breast cancer, assign a follow-up time of 5 years;
				if bccause eq 0 then mortime2 = 5;
*for non-reidentified women who were assigned death due to breast cancer, assign a follow-up time drawn
 from the distribution observed among reidentified women who died of breast cancer (bounded by 0 and 5 years);
					else if bccause eq 1 then do;
					mortime2 = max(0,min(5,2.63 + 1.37*rannor(seed2)));
					end;
					output;
				end;
				repcnt = repcnt+1;
				end;
				stop;
				run;
*combine the data set of reidentified women with the data set of non-reidentified women;
		data all_id1;
*change the number of iterations;
			repcnt=1; reptot=2000;
*step is the number of records in the non-reidentified database;
			step = last2/reptot;
*write a group of records from the reidentified database.  These records do not change with each iteration;
			do while (repcnt le reptot);
				do reccnt1 = 1 to last1;
					set yes_reid1 point=reccnt1 nobs=last1;
					output;
				end;
*calculate where to begin and end writing records from the database of non-reidentified women;
*these records do not change with each iteration;
			start=((repcnt-1)*step)+1;stop=repcnt*step;
*write a group of records from the non-reidentified database;
				do reccnt2 = start to stop;
					set no_reid2 point=reccnt2 nobs=last2;
					output;
				end;
			repcnt = repcnt+1;
			end;
			stop;
			run;
*the next section accounts for misclassification of stage;
*create a dataset of all reidentified women, read by the iteration number;
			data all_id2; set all_id1; by repcnt;
*retain seed variables for random number generators and misclassification variables;
			retain seeda1 seeda2 seeda3 ppv npv sens spec;
*initialize the seed variables on the first record;
			if seeda1 eq . then do;
            	seeda1=-1;
            	seeda2=-1;
				seeda3=-1;
			end;
*to disable the misclassification component of the sensitivity analysis add a 
 begin comment delimiter (/*) after this comment;

*initialize the values of the misclassification parameters each time a new iteration begins;
*values of the misclassification paramters derived from literature review;
*prevalence of regional disease is between 0.3 and 0.5, with median 139/342 (the prevalence among
 women staged pathologically in this study;
*the remaining fractions are abstracted from literature reports of the sensitivity and specificity of
 clinical staging versus the gold standard pathologic staging;
			if first.repcnt then do;
						prev = ((.5-.3)*rantri(seeda1,(139/342-.3)/(.5-.3))+.3);
                        sens = ((60/72-2/15)*rantri(seeda1,(995/2155-2/15)/(60/72-2/15))+2/15);
                        spec = ((152/155-24/58)*rantri(seeda2,(3001/3522-24/58)/(152/155-24/58))+24/58);
                        ppv = ((prev)*sens)/((prev)*sens + (1-prev)*(1-spec));
                        ppv = max(0,min(1,ppv));
                        npv = ((1-prev)*spec)/((1-prev)*spec + (prev)*(1-sens));
                        npv = max(0,min(1,npv));
            end;
*for women who did not receive an axillary dissection, implement the misclassification analysis;
        if axil eq 0 then do;
*for women staged clinically with local disease, examine whether to reclassify as regional;
        	if excat1 eq 0 then
                call ranbin(seeda1,1,(1-npv),excat1);
*for women staged clinically with regional disease, examine whether to reclassify as local;
        	else if excat1 eq 1 then
                call ranbin(seeda2,1,ppv,excat1);
*determine whether therapy was definitive, given potentially reclassified stage;
*definitive therapy for local disease;
        	if excat1 eq 0 then do;
                if mastsurg eq 1 then defther1 = 0;
                        else if bcssurg eq 1 and ttorad ne . and ttorad lt 152 then defther1 = 0;
                        else defther1 = 1;
        	end;
*definitive therapy for regional disease;
        	else if excat1 eq 1 then do;
                if mastsurg eq 1 and adjchem eq 1 and stopped eq 0 then defther2 = 0;
                else if bcssurg eq 1 and ttorad ne . and ttorad lt 152  and adjchem eq 1 then defther2 = 0;
                else defther2 = 1;
        	end;
        	if defther1 eq 0 or defther2 eq 0 then defnther = 0;
                else if defther1 eq 1 or defther2 eq 1 then defnther = 1;
		end;

*to end the disabled misclassification component of the sensitivity analysis add an 
 end comment delimiter (*/) before this comment;

*to disable the unknown confounder component of the sensitivity analysis add a 
 begin comment delimiter (/*) after this comment;

*calculate the prevalence of the unknown confounder among those with definitive therapy and who 
 did not die of breast cancer;
		base = 0.1*ranuni(seeda1)+0.3;
*calculate the additional prevalence among those with less than definitive therapy;
		adddef = 0.1*ranuni(seeda1)+0.1;
*calculate the additional prevalence among those who died of breast cancer;
		addbc = 0.1*ranuni(seeda1)+0.1;
*calculate the prevalence of the unknown confounder, given the values calculated above
 and the exposure and disease characteristics of the observation;
		punkwn = base + adddef*defnther + addbc*bccause;
*given the prevalence, call ranbin to determine unknown confounder status;
		call ranbin(seeda3,1,punkwn,unknown);

*to end the disabled unknown confounder component of the sensitivity analysis add an 
 end comment delimiter (*/) before this comment;
*keep only the variables necessary for the estimation of relative hazards;
        keep mortime2 bccause defnther excat1 agecat1 agecat2 
			  repcnt sens spec ppv npv unknown punkwn; 
		run;
*proportional hazards modeling for sensitivity analysis only;
*the by repcnt statement calculates an estimate of effect for each iteration;	
        proc phreg data=all_id2 outest=ther noprint covout;
                model mortime2*bccause(0) = defnther excat1 agecat1 agecat2 unknown;
				by repcnt;
*the data set pparam contains only the parameter estimates;
		data pparam (where=(_type_='PARMS') rename=(defnther=defnparm)); 
			set ther; run;
*the data set pstderr keeps only the variance estimates for the definitive therapy parameter;
		data pstderr (keep=_type_ _name_ repcnt defnther 
		where=(_type_='COV' and _name_='DEFNTHER')); 
			set ther; run;
*the data set paramstd keeps the paramter estimate and variance for the definitive therapy variable;
		data paramstd; 
			merge pparam pstderr (keep=repcnt defnther); 
			by repcnt; run; 
*phreg for true bootstrapping data sets to analyze sensitivity analysis and random error;
*create bootstrapping data set (boot1) by randomly selecting with replacement 449 records from each iteration group;
		data boot1;
			retain seed1;
*change the number of iterations;
			repcntb=1; reptotb=2000;
*initialize the random number generator seed on first record;
    		if seed1 eq . then seed1 = -1;
*the first do loop counts iterations of the sensitivity analysis;
			 	do while (repcntb le reptotb);
					ncountb = 1;
*the second do loop counts records selected with replacement until the same number of records as the original
 data set have been selected;
	        		do while (ncountb le 449);
*lookat is the record number to select;
			        	lookat = (repcntb-1)*449+round(449*ranuni(seed1),1)+1;
*the set statement reads the record at lookat, then it is output;
            			set all_id2 point=lookat;
                		output;
            			ncountb = ncountb+1;
            			end;
					repcntb=repcntb+1;
					end;
    			stop;
			run;
*proportional hazards modeling for sensitivity analysis and random error in the bootstrapped dataset;
*the by repcnt statement calculates an estimate of effect for each iteration;
        proc phreg data=boot1 outest=therb noprint;
                model mortime2*bccause(0) = defnther excat1 agecat1 agecat2 unknown;
				by repcntb;
*transfer the parameter estimate from the output dataset (therb) to a second dataset (trans1b),
 changing the names of the variables that are kept so they can be merged with the paramstd dataset;
		data trans1b (keep=paramadj repcnt); 
			set therb (rename=(defnther=paramadj repcntb=repcnt)); run;
*sort and merge the two parameter data sets;
		proc sort data=paramstd; by repcnt;
		proc sort data=trans1b; by repcnt;
        data combineb; 
			merge paramstd trans1b; by repcnt;
*initialize the random number generator seeds;
		retain seed1 seed2 seed3; if seed1 eq . then do; 
			seed1 = -1; seed2 = -1; seed3 = -1; 
		end;
*create the parameters of the bayesian prior distribution;
*bayeshr is the central tendency of the prior distribution and bayesse is the standard error
 of hte prior distribution;
*change the prior parameters;
		bayeshr = 1.5;
		bayesse = 0.2027;
		label
			expadj = 'Relative hazard, all error'
			expparm = 'Relative hazard, systematic error'
			exporg = 'Relative hazard, conventional analysis';
*create a distribution to represent the result without sensitivity analysis;
		paramorg = log(2) + rannor(seed1)*0.27755;
*create the bayesian prior distribution, assigning 50% probability to null effect and 
 50% of probability distribution to the prior specified above;
*change the upper limit (0.5) to 0, for example, to assign all of the probability to 
 the parameterized prior;
		priortst = ranuni(seed1);
			if 0 le priortst le 0.5 then pprior = 0;
				else pprior = max(0,log(bayeshr) + rannor(seed3)*bayesse);
		expadj = exp(paramadj);
		expparm = exp(defnparm);
		exporg = exp(paramorg);
		run;
*if necessary, combine multiple invocations of the sensitivity analysis by un-commenting the following code;
*the result will be to accumulate 2000 iterations (for example) per invocation;
*verify that no data set named sensdir.combined has been created before beginning a new analysis;
/*
*combining multiple sets;
		data sensdir.combined; set sensdir.combined combineb; 
		run;
*/
*use proc univariate to provide parameters necessary to combine prior distribution and 
 sensitivity analysis distribution.  Write results to dataset combinec;
*change combineb to sensdir.combined if multiple invocations;
		proc univariate data=combineb noprint;
			var defnparm paramadj paramorg pprior;
			output out=forbayes
				mean=parmmean padjmean porgmean pprimean
				var=parmvar padjvar porgvar pprivar;
		run;
*use the results of the proc univariate (combinec) concatenated with original results to
 generate the posterior distribution;
*change combineb to sensdir.combined if multiple invocations;
		data combinec;
			if _n_=1 then set forbayes;
			set combineb;
			retain seedb1; 
			if seedb1 eq . then seedb1 = -1;
			ppost = (paramadj/padjvar+ pprior/pprivar)/(1/padjvar+1/pprivar);
			exppost = exp(ppost);
			expprior = exp(pprior);
		run;
*use proc univariate to generate tabular results, written to data set pctls;
*note that the intervals about the percentiles can be used to ascertain convergence;
		proc univariate data=combinec cipctldf (type=asymmetric);
			var defnparm paramadj paramorg pprior ppost;
			output out=pctls
				min=parmmin padjmin porgmin pprimin ppstmin 
				mean=parmmean padjmean porgmean pprimean ppstmean
				median=parmmed padjmed porgmed pprimed ppstmed
				max=parmmax padjmax porgmax pprimax ppstmax
				var=parmvar padjvar porgvar pprivar ppstvar
				pctlpts = 2.5 97.5 pctlpre = parm padj porg ppri ppst
				pctlname = p025 p975;
		run;
*label the tabular results;
		data pctls; set pctls;
		label
		parmmean	=	mean, sensitivity analysis
		padjmean	=	mean, bootstrap analysis
		porgmean	=	mean, original analysis
		pprimean	=	mean, prior distribution
		ppstmean	=	mean, posterior distribution
		parmvar	=	variance, sensitivity analysis
		padjvar	=	variance, bootstrap analysis
		porgvar	=	variance, original analysis
		pprivar	=	variance, prior distribution
		ppstvar	=	variance, posterior distribution
		parmmax	=	maximum, sensitivity analysis
		padjmax	=	maximum, bootstrap analysis
		porgmax	=	maximum, original analysis
		pprimax	=	maximum, prior distribution
		ppstmax	=	maximum, posterior distribution
		parmmed	=	median, sensitivity analysis
		padjmed	=	median, bootstrap analysis
		porgmed	=	median, original analysis
		pprimed	=	median, prior distribution
		ppstmed	=	median, posterior distribution
		parmmin	=	minimum, sensitivity analysis
		padjmin	=	minimum, bootstrap analysis
		porgmin	=	minimum, original analysis
		pprimin	=	minimum, prior distribution
		ppstmin	=	minimum, posterior distribution
		parmp025	=	2.5%, sensitivity analysis
		parmp975	=	97.5%, sensitivity analysis
		padjp025	=	2.5%, bootstrap analysis
		padjp975	=	97.5%, bootstrap analysis
		porgp025	=	2.5%, original analysis
		porgp975	=	97.5%, original analysis
		pprip025		=	2.5%, prior distribution
		pprip975	=	97.5%, prior distribution
		ppstp025	=	2.5%, posterior distribution
		ppstp975	=	97.5%, posterior distribution;
		run;
*print the tabular results;
		proc means data=pctls median;
			var parmmean parmvar parmmin parmp025 parmmed parmp975 parmmax padjmean padjvar padjmin padjp025 padjmed
				padjp975 padjmax porgmean porgvar porgmin porgp025 porgmed porgp975 porgmax pprimean pprivar pprimin
				pprip025 pprimed pprip975 pprimax ppstmean ppstvar ppstmin ppstp025 ppstmed ppstp975 ppstmax; 	
		run;
*use proc rank to calculate the percentiles that will be plotted on the cumulative
 probability functions, outputting results to dataset combine2;
		proc rank data=combinec out=combine2 percent ties=mean;
			ranks rparam rparmadj rparmorg rparmpri rparmpst;
			var defnparm paramadj paramorg pprior ppost;
		run;
*label the output of proc ranks for legends on the plots;	
		data combine2; set combine2;
			label
				rparam = 'Sensitivity analysis (SA)'
				rparmadj = 'Bootstrapped SA'
				rparmorg = 'Conventional result'
				rparmpri = 'Bayesian prior'
				rparmpst = 'Bayesian posterior';
		run;				
*symbol definitions to create the lines for the plots;
		symbol1 ci=bl interpol=splineps line=20 value=none width=12;
		symbol2 ci=bl interpol=splineps line=1 value=none width=16;
		symbol3 ci=bl interpol=splineps line=41 value=none width=12;
		symbol4 ci=red interpol=splineps line=1 value=none width=12;
		symbol5 ci=red interpol=splineps line=20 value=none width=12;
*plot the sensitivity analysis, bootstrap, and original result on the arithmetic scale;
		proc gplot data=combine2 uniform;
			plot rparam*defnparm rparmadj*paramadj rparmorg*paramorg 
				/ legend=legend1 noframe
				overlay chref=gray cvref=gray vref=2.5 50 97.5 
				vaxis=axis1 haxis=axis2;
		run;
*plot the sensitivity analysis, bootstrap, and original result on the logarithmic scale;
		proc gplot data=combine2 uniform;
			plot rparam*expparm rparmadj*expadj rparmorg*exporg 
				/ legend=legend1 noframe
				overlay chref=gray cvref=gray vref=2.5 50 97.5 
				vaxis=axis1 haxis=axis3;
		run;				
*symbol definitions for the Bayesian plot;
		symbol1 ci=bl interpol=splines line=1 value=none width=16;
		symbol2 ci=bl interpol=splines line=41 value=none width=12;
		symbol3 ci=bl interpol=splines line=20 value=none width=12;
*plot the prior, bootstrap, and posterior on the arithmetic scale;
		proc gplot data=combine2 uniform;
			plot  rparmadj*paramadj 
				rparmpri*pprior rparmpst*ppost 
				/ legend=legend1 noframe
				overlay chref=gray cvref=gray vref=2.5 50 97.5 
				vaxis=axis1 haxis=axis2;
		run;
*plot the prior, bootstrap, and posterior on the logarithmic scale;
		proc gplot data=combine2 uniform;
			plot rparmadj*expadj 
				rparmpri*expprior rparmpst*exppost 
				/ legend=legend1 noframe
				overlay chref=gray cvref=gray vref=2.5 50 97.5 
				vaxis=axis1 haxis=axis3;
		run;
		quit;



		