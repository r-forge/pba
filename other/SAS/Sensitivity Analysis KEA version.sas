********** 2005 **********************************************************  version 1.1  *****;
***************        SENSITIVITY ANALYSIS MISCLASSIFICATION MACRO		    ******************;
**********************************************************************************************;

* this is the code that contains the macros that compute the sensitivity analysis;

* change to symbolgen to see the macro variables in the log, and clears the work library;
	options obs = max nosymbolgen spool ps=40 ls=124 nodate;

* KEA: logout is the filename where to store log output;
* dummy file to turn the SAS log on and off;
	*filename dumfl dummy;

* general macro to turn off log;
%macro logoff();
	proc printto log=&logfile; 
	run;

	* KEA: output is written to the file specified in outputfile;
	proc printto print=&outputfile;
	run;
%mend;

* general macro to turn on log;
%macro logon();
	proc printto; 
	run;
%mend;

* general macro to delete data sets;
%macro delset(libn,datset);			
	proc datasets library=&libn NOLIST NODETAILS;
		delete &datset;
	run;
%mend;

**********************************************************************************************;
*********************                  MACRO RUN                       ***********************;
**********************************************************************************************;
%macro sensmac(libname,startover,outset,totalreps,inset,log,logfile,outputfile,depend,exp,indep,miscvar,
			   mis_ind,sens_min,sens_max,sens_mod,sens_mod2,spec_min,spec_max,spec_mod,spec_mod2,
				sens_min_d,sens_max_d,sens_mod_d,sens_mod_d2,spec_min_d,spec_max_d,spec_mod_d,spec_mod_d2,
				corrsens,corrspec,strata_id);


	* run the macro parameter checks macro (at bottom);
		%checkerr;

	* set repetitions high if none specified otherwise add one for logistics;;
		%if &totalreps=0 %then %let totalreps=30000;
	
		%else %let totalreps=%eval(&totalreps+1);

	* if terminal errors found, terminate the macro;
		%if &quitmac %then %goto miscterm;

	* create libname;
		libname sensdir "&libname";

	* turn off the log if log option set to 1 or 2;
		%if &log = 2 or &log = 1 %then %logoff;	

	* delete old repetitions if startover set to yes;
		%if &startover = yes %then %do;
			%delset(sensdir,&outset drop_&outset converg_&outset);
		%end;

	* delete old temporary data sets if they exists;
		%delset(work,mistot combine2 combineb combinec combineplots conventional conventional1 converg converge1 converge2 );
		%delset(work,converge3 converge4 corrs count crude foostats misclas misclass new outpval pctls sstemp systematic);
		%delset(work,sys_parm sys_parms sys_var test2set testset univpcts);

	* turn on the log if log set to 1;
		%if &log = 1 %then %logon;

	* execute the macro;
	* create a loop counter (j) and set = 1;
		%global j;
		%let j=1;

		* loop until the number of loops is greater than the specified amount;
			%do %while (&j < &totalreps);

			* runs the sensitivity analysis;
				%miscsens(&inset,&depend,&exp,&indep,&log,&miscvar,&sens_min,
						&sens_max,&sens_mod,&sens_mod2,&spec_min,&spec_max,&spec_mod,&spec_mod2,
						&sens_min_d,&sens_max_d,&sens_mod_d,&sens_mod_d2,&spec_min_d,&spec_max_d,&spec_mod_d,
						&spec_mod_d2,&mis_ind,&corrsens,&corrspec);

			* check that iteration worked, then execute;
				%if &negat ne 1 %then %do;
				* creates the random error and bootstrap approximation datasets;	
					%randerr(&totalreps,&outset,&inset,&j,&depend,&exp,&indep,&log);
					%let j=%eval(&j+1);
				
				* checks convergence every repetition after the 20th or at last rep;
					%if &j > 20 or &j=&totalreps %then %conver(&totalreps);
				%end;	
				%if &negat=1 %then %do;
				* delete old temporary data sets if they exists;
					%delset(work,mistot combine2 combineb combinec combineplots conventional conventional1 converg converge1 converge2 );
					%delset(work,converge3 converge4 corrs count crude foostats misclas misclass new outpval pctls sstemp systematic);
					%delset(work,sys_parm sys_parms sys_var test2set testset univpcts);	
				%end;
		%end;

	* create graphs of sensitivity analysis;
		%graph(&exp);

	* turns back on the log if turned off;
		%if &log = 2 or &log = 1 %then %logon;

	* creates a marker for the line to goto if errors found;
		%miscterm:;

	* if macro was terminated, tell user in log;
		%if &quitmac=1 %then %put ERROR:  TERMINATING SENSITIVITY MACRO;
%mend;

*************************** MISCLASSIFICATION SECTION ***************************************;
%macro miscsens(inset,depend,exp,indep,log,miscvar,sens_min,
						sens_max,sens_mod,sens_mod2,spec_min,spec_max,spec_mod,spec_mod2,
						sens_min_d,sens_max_d,sens_mod_d,sens_mod_d2,spec_min_d,spec_max_d,
						spec_mod_d,spec_mod_d2,mis_ind,corrsens,corrspec);

	* create a marker for illogical repititions;
		%global negat;
		%let negat=0;

	* create a dataset and set the dummy counter (repcnt) for synchronization to 1;
		data misclas; set &inset;
			repcnt=1;
		run;

	* sort the dataset;
		proc sort data=misclas; by repcnt; run;

	* create dataset of observed collapsed data and output;
		proc freq data=&inset noprint;
			tables &mis_ind*&miscvar/out=crude;
		run;

	* collapse crude to 1 observation;
		data crude; set crude; where &miscvar ne .;
			retain ed_o eud_o ued_o ueud_o;
			if &miscvar=1 and &mis_ind=1 then ed_o=count;
			if &miscvar=0 and &mis_ind=1 then ued_o=count;
			if &miscvar=1 and &mis_ind=0 then eud_o=count;
			if &miscvar=0 and &mis_ind=0 then ueud_o=count;
			repcnt=1;
			if _N_=4 then output;
		run;

	* merge the collapsed data into each observation in the raw data;
		proc sort data=misclas; by repcnt; run;
		proc sort data=crude; by repcnt; run;

		data misclas; merge misclas crude; by repcnt;
		run;

	* set the misclassification data set by the repcnt variable;
		data misclas; set misclas; by repcnt;
			* retain seeds for random number generators and misclassification variables;
				retain seeda seedb seedc seedd seede seedf seedg seedh sens spec sens2 spec2;
			* initialize the seed variables at the first record;
				if seeda eq . then do; 
					seeda=-1; seedb=-1; seedc=-1; seedd=-1; seede=-1; seedf=-1;seedg=-1; seedh=-1; 
				end;

			* initialize values of the misclas parameters each time a new iteration begins;
				if first.repcnt then do; 
				* create a trapezoidal distribution based on input parameters and choose; 
				* sensitivity from it;
					p=ranuni(seeda);
					sens =  (p*(&sens_max+&sens_mod2-&sens_min-&sens_mod) + (&sens_min + &sens_mod))/2;
					if sens < &sens_mod then do; 
						sens  = &sens_min + sqrt((&sens_mod-&sens_min)*(2*sens-&sens_min-&sens_mod));
					end;
				   	else if sens > &sens_mod2 then do; 
						sens = &sens_max - sqrt(2*(&sens_max-&sens_mod2)*(sens-&sens_mod2));
					end;
				* create a trapezoidal distribution based on input parameters and choose; 
				* specificity from it;
					pt=ranuni(seedb);
					spec =  (pt*(&spec_max+&spec_mod2-&spec_min-&spec_mod) + (&spec_min + &spec_mod))/2;
					if spec < &spec_mod then do; 
						spec  = &spec_min + sqrt((&spec_mod-&spec_min)*(2*spec-&spec_min-&spec_mod));
					end;
				   	else if spec > &spec_mod2 then do; 
						spec = &spec_max - sqrt(2*(&spec_max-&spec_mod2)*(spec-&spec_mod2));
					end;

				* if nondifferential sensitivity, set equal to each other;
					if &ndsens = 1 then sens2=sens;
				* if not, create correlated sensitivities;
					%if &ndsens ne 1 %then %do;

					* creates correlated random variables;
						u=ranuni(seedc);
						e1=ranuni(seedd);
						e0=ranuni(seede);

						t = log(u/(1-u));
						
						f1 = log(e1/(1-e1));
						f0 = log(e0/(1-e0));

						p1 = exp(sqrt(&corrsens)*t+sqrt(1-&corrsens)*f1) / (1+(exp(sqrt(&corrsens)*t+sqrt(1-&corrsens)*f1)));
						p0 = exp(sqrt(&corrsens)*t+sqrt(1-&corrsens)*f0) / (1+(exp(sqrt(&corrsens)*t+sqrt(1-&corrsens)*f0)));

					* uses correlated variables p1 and p0 to choose from trapezoidals;
						sens =  (p1*(&sens_max+&sens_mod2-&sens_min-&sens_mod) + (&sens_min + &sens_mod))/2;
						if sens < &sens_mod then do; 
							sens  = &sens_min + sqrt((&sens_mod-&sens_min)*(2*sens-&sens_min-&sens_mod));
						end;
					   	else if sens > &sens_mod2 then do; 
							sens = &sens_max - sqrt(2*(&sens_max-&sens_mod2)*(sens-&sens_mod2));
						end;					
						sens2 =  (p0*(&sens_max_d+&sens_mod_d2-&sens_min_d-&sens_mod_d) + (&sens_min_d + &sens_mod_d))/2;
						if sens2 < &sens_mod_d then do; 
							sens2  = &sens_min_d + sqrt((&sens_mod_d-&sens_min_d)*(2*sens2-&sens_min_d-&sens_mod_d));
						end;
					   	else if sens2 > &sens_mod_d2 then do; 
							sens2 = &sens_max_d - sqrt(2*(&sens_max_d-&sens_mod_d2)*(sens2-&sens_mod_d2));
						end;
					%end;

				* if nondifferential specificity, set equal to each other;
					if &ndspec = 1 then spec2 = spec;
				* if not, create correlated specificites;
					%if &ndspec ne 1 %then %do;

					* creates correlated random variables;
						x=ranuni(seedf);
						e3=ranuni(seedg);
						e2=ranuni(seedh);

						v = log(x/(1-x));
						
						f3 = log(e3/(1-e3));
						f2 = log(e2/(1-e2));

						p3 = exp(sqrt(&corrspec)*v+sqrt(1-&corrspec)*f3) / (1+(exp(sqrt(&corrspec)*v+sqrt(1-&corrspec)*f3)));
						p2 = exp(sqrt(&corrspec)*v+sqrt(1-&corrspec)*f2) / (1+(exp(sqrt(&corrspec)*v+sqrt(1-&corrspec)*f2)));

						* uses correlated variables p3 and p2 to choose from trapezoidals;
							spec =  (p3*(&spec_max+&spec_mod2-&spec_min-&spec_mod) + (&spec_min + &spec_mod))/2;
							if spec < &spec_mod then do; 
								spec  = &spec_min + sqrt((&spec_mod-&spec_min)*(2*spec-&spec_min-&spec_mod));
							end;
						   	else if spec > &spec_mod2 then do; 
								spec = &spec_max - sqrt(2*(&spec_max-&spec_mod2)*(spec-&spec_mod2));
							end;					
							spec2 =  (p2*(&spec_max_d+&spec_mod_d2-&spec_min_d-&spec_mod_d) + (&spec_min_d + &spec_mod_d))/2;
							if spec2 < &spec_mod_d then do; 
								spec2  = &spec_min_d + sqrt((&spec_mod_d-&spec_min_d)*(2*spec2-&spec_min_d-&spec_mod_d));
							end;
						   	else if spec2 > &spec_mod_d2 then do; 
								spec2 = &spec_max_d - sqrt(2*(&spec_max_d-&spec_mod_d2)*(spec2-&spec_mod_d2));
							end;
					%end;
				end;

		* create titles for output;
			%global titles1 titles2 titles3 titles4;
			%let titles1=SENS 1 (min=&sens_min, mod1=&sens_mod, mod2=&sens_mod2, max=&sens_max);
			%let titles2=SPEC 1 (min=&spec_min, mod1=&spec_mod, mod2=&spec_mod2, max=&spec_max);
			%if &ndsens  ne 1 %then %let titles3=, SENS 2 (min=&sens_min_d, mod1=&sens_mod_d, mod2=&sens_mod_d2, max=&sens_max_d);
			%else %let titles3=, NON-DIFFERENTIAL;
			
			%if &ndspec ne 1 %then %let titles4=, SPEC 2 (min=&spec_min_d, mod1=&spec_mod_d, mod2=&spec_mod_d2, max=&spec_max_d);
			%else %let titles4=, NON-DIFFERENTIAL;

		* calculate expected truth based on chosen sens and spec;
				ed_t=(ed_o-(1-spec)*(ed_o+ued_o))/(sens-(1-spec));
				ued_t=(ed_o+ued_o)-ed_t;
				eud_t=(eud_o-(1-spec2)*(eud_o+ueud_o))/(sens2-(1-spec2));
				ueud_t=(eud_o+ueud_o)-eud_t;

		* from expected truth, calculate false positives and false negatives
		  and true positives and true negatives for group1;
				t_pos_d=(sens*ed_t);
				t_neg_d=(spec*ued_t);
				f_neg_d=((1-sens)*ed_t);
				f_pos_d=((1-spec)*ued_t);

		* from expected truth, calculate false positives and false negatives
		  and true positives and true negatives for group2;
				t_pos_ud=(sens2*eud_t);
				t_neg_ud=(spec2*ueud_t);
				f_neg_ud=((1-sens2)*eud_t);
				f_pos_ud=((1-spec2)*ueud_t);

		* from expected truth, calculate NPV and PPV for group1 and group2;
				PPV_d=t_pos_d/(t_pos_d+f_pos_d);
				NPV_d=t_neg_d/(t_neg_d+f_neg_d);
				PPV_ud=t_pos_ud/(t_pos_ud+f_pos_ud);
				NPV_ud=t_neg_ud/(t_neg_ud+f_neg_ud); 

		* if negative numbers then terminate;
				if 0 > PPV_d or PPV_d > 1 or 0 > NPV_d or NPV_d > 1 or 0 > PPV_ud or PPV_ud > 1 or 
				0 > NPV_ud or NPV_ud > 1 then call symput("negat",1);

				label sens='Sensitivity1' sens2='Sensitivity2'
					  spec='Specificity1' spec2='Specificity2'; 
			run;

	* check that parameters chosen are logical;
		%if &negat=1 %then %do;
				%put WARNING: Sens and Spec produced cells with negative numbers;
				%put Terminating iteration;
			* save dropped iterations in a data set;
				data misclass; set misclas; 
					contype="Negative cell counts";
					label contype="Type of error";
					if _N_=1 then output; 
					drop ed_o eud_o eud_o ueud_o u e1 e0 t f2 f1 f0 p1 p0 x e3 e2 v f3 f3 p3 p2 i case exp repcnt count percent seeda seedb seedc seedd seede seedf seedg seedh p pt;
				run;
				proc append base=sensdir.drop_&outset data=misclass force; run;
		%end; 
		%if &negat=1 %then %goto negterm;

	data misclas; set misclas; 
		retain seeda2;
			* initialize the seed variable at the first record;
				if seeda2 eq . then seeda2=-1;

	* reclassify those with &mis_ind = 0;
		if &mis_ind eq 0 then do;
		* reclassify those who have &miscvar = 0; 
			if &miscvar eq 0 then do;
			* check that NPV is not = 0 or 1;
				if npv_ud eq 0 then reclased = 1;
				else if npv_ud eq 1 then reclased = 0; 
			* determine who needs to be reclassified;
				else do;	
					mc=ranuni(seeda2)>npv_ud;
					reclased=mc=1;
				end;
			end;

			* reclassify those who have &miscvar = 1;
			else if &miscvar eq 1 then do;
				* check that PPV is not = 0 or 1;
				if ppv_ud eq 0 then reclased = 0;
				else if ppv_ud eq 1 then reclased = 1;
			* determine who needs to be reclassified;
				else do;	
					mc=ranuni(seeda2)>ppv_ud;
					reclased=mc=0;
				end;
			end;
		end;

		* reclassify those with &mis_ind = 1; 
		else if &mis_ind eq 1 then do;
			* reclassify those who have &miscvar = 0; 
			if &miscvar eq 0 then do;
				* check that CNPV is not = 0 or 1;
				if npv_d eq 0 then reclased = 1;
				else if npv_d eq 1 then reclased = 0; 
				* determine who needs to be reclassified;
				else do;	
					mc=ranuni(seeda2)>npv_d;
					reclased=mc=1;
				end;
			end;

			* reclassify those who have &miscvar = 1;
			else if &miscvar eq 1 then do;
				* check that CPPV is not = 0 or 1;
				if ppv_d eq 0 then reclased = 0;
				else if ppv_d eq 1 then reclased = 1;
				* determine who needs to be reclassified;
				else do;	
					mc=ranuni(seeda2)>ppv_d;
					reclased=mc=0;
				end;
			end;
		end;
	* correct the variable;
		&miscvar=reclased;

	* clean up data set;
		drop count percent seeda seedb seedc seedd seede seedf seedg seedh;
	run;

	proc means data=misclas noprint; class &depend &exp; output out=zero n=n; run;
	data zero; set zero; where _TYPE_=3; run;
	proc means data=zero noprint; var n; output out=zeros n=n; run;
	data zeros; set zeros; call symput ("zero",n);run;

		%if %eval(&zero)<4 %then %do;
				%let negat=1;
				%put WARNING: zero cells found;
				%put Terminating iteration;
			* save dropped iterations in a data set;
				data misclass; set misclas; 
					contype="Zero cells found";
					label contype="Type of error";
					if _N_=1 then output; 
					drop ed_o eud_o eud_o ueud_o u e1 e0 t f2 f1 f0 p1 p0 x e3 e2 v f3 f3 p3 p2 i case exp repcnt count percent seeda seedb seedc seedd seede seedf seedg seedh p pt;
				run;
				proc append base=sensdir.drop_&outset data=misclass force; run;
		%end; 
		%if %eval(&zero)<4 %then %goto negterm;	

* create dataset for saving the chosen sensitivity and specificity;
	proc means data=misclas noprint; var sens spec sens2 spec2;
		output out=sstemp mean=sens spec sens2 spec2;
	run;

* make chosen values macro variables;
	data sstemp; set sstemp; 
		call symput("SE_a",sens);
		call symput("SE_b",sens2);
		call symput("SP_a",spec);
		call symput("SP_b",spec2);
		drop _TYPE_ _FREQ_; 
	run;

	%put Iteration &j;
	%put Chosen values -- Sensitivity1: &SE_a , Sensitivity2: &SE_b;
	%put Chosen values -- Specificity1: &SP_a , Specificity2: &SP_b;

	* KEA: also make dataset containing the current iteration number;
	data iteration_number;
		current_iteration = &j;
	run;
	libname lib "&libname";
	data lib.iteration_number;
		set iteration_number;
	run;
	* KEA: end;

********************************* MODELED OUTCOMES ********************************************;
* runs the logistic regression with the current data set;

	* run a logistic regression model for the corrected misclas data set;
	/* KEA: use phreg to make conditional logistic regression - note that an additional macro 
	  variable &strata_id must be specified

	  KEA: because phreg doesn't have a descending option	*/
	data misclas;
		set misclas;
		if &depend=0 then case=1; else case=0;
	run;
	proc phreg data=misclas nosummary;
		model case = &exp &indep / covb ties=discrete;
		strata &strata_id;
		ods output ParameterEstimates=ParaEst CovB=Cov ConvergenceStatus=ConvStatus;
	run;
	* to make a dataset like 'systematic' from PROC LOGISTIC;
	proc transpose data=ParaEst out=ParaEst2;
		id variable;
	run;
	data ParaEst2;
		set ParaEst2;
		if _name_ = 'Estimate';
		drop _name_ _label_;
	run;
	data ParaEst2;
		format _NAME_ $26.;
		set ParaEst2;
		_TYPE_ = 'PARMS';
		_NAME_ = 'CASE';
	run;
	data cov2;
		set cov;
		_NAME_ = RowName;
		_TYPE_ = 'COV';
		drop RowName;
	run;
	data cov_paraest;
		set ParaEst2 cov2;
	run;
	data ConvStatus2;
		set ConvStatus;
		_STATUS_ = status;
		drop status reason;
	run;
	* systematic now mimics the output from PROC LOGISTIC;
	data systematic;
		if _n_ = 1 then set ConvStatus2;
		set cov_paraest;
	run;
	* KEA: end;

	* remove observations if logistic model does not converge, and tell user;
		data systematic; set systematic;
	* KEA: Convergence output from phreg is different;
			*if _STATUS_='1 Warning' then do;
			if _STATUS_^=0 then do;
					delete;
					call symput("negat",1);
			end;
		run;

		%if &negat=1 %then %do;
				%put WARNING: REMOVING OBSERVATIONS THAT DO NOT CONVERGE (DATA SEPARATION) - RESULTS MAY BE INACCURATE;
				%put Terminating iteration;
			* save dropped iterations in a data set;
				data misclass; set misclas; 
					contype="Non-convergence";
					label contype="Type of error";
					if _N_=1 then output; 
					drop ed_o eud_o eud_o ueud_o i case exp repcnt count percent seeda seedb seedc seedd seede seedf seedg seedh p pt;
				run;
				proc append base=sensdir.drop_&outset data=misclass force; run;
		%end; 
		%if &negat=1 %then %goto negterm;

	* the data set sys_parm contains only the parameter estimates;
		data sys_parm (where=(_type_='PARMS') rename=(&exp=sys_parm)); 
			set systematic; 
		run;

	* the data set sys_var keeps only the variance estimates for the outcome parameter;
		data sys_var; set systematic (keep=_type_ _name_ &exp); 
			where (_type_='COV' and upcase(_name_)=upcase(resolve('&exp')));
			sys_var=&exp;
		run;

	* the data set sys_parms keeps the parameter estimate and variance for the outcome variable;
		data sys_parms; 
			merge sys_parm sys_var (keep = sys_var); 
		run;

	%negterm:;
%mend;

%macro randerr(totalreps,outset,inset,j,depend,exp,indep,log,clearout);

	/* KEA: use phreg to make conditional logistic regression - note that an additional macro 
	  variable &strata_id must be specified
	  KEA: because phreg doesn't have a descending option	*/

	data inset;
		set &inset;
		if &depend=0 then case=1; else case=0;
	run;
	proc phreg data=inset nosummary;
		model case = &exp &indep / covb ties=discrete;
		strata &strata_id;
		ods output ParameterEstimates=ParaEst_conven CovB=Cov_conven ConvergenceStatus=ConvStatus_conven;
	run;
	* to make a dataset like 'conventional1' from PROC LOGISTIC;
	proc transpose data=ParaEst_conven out=ParaEst_conven2;
		id variable;
	run;
	data ParaEst_conven2;
		set ParaEst_conven2;
		if _name_ = 'Estimate';
		drop _name_ _label_;
	run;
	data ParaEst_conven2;
		format _NAME_ $26.;
		set ParaEst_conven2;
		_TYPE_ = 'PARMS';
		_NAME_ = 'CASE';
	run;
	data cov_conven2;
		set cov_conven;
		_NAME_ = RowName;
		_TYPE_ = 'COV';
		drop RowName;
	run;
	data cov_paraest_conven;
		set ParaEst_conven2 cov_conven2;
	run;
	data ConvStatus_conven2;
		set ConvStatus_conven;
		_STATUS_ = status;
		drop status reason;
	run;
	* conventional1 now mimics the output from PROC LOGISTIC;
	data conventional1;
		if _n_ = 1 then set ConvStatus_conven2;
		set cov_paraest_conven;
	run;
	* KEA: end;

	* create data set of parameter and variance for the variable of interest;
		data conventional; set conventional1;
			if _TYPE_="PARMS" or upcase(_name_)=upcase(resolve('&exp'));
		run;

	* place parameter and variance into macro vars to create distribution for 
	  the original analysis;
		data conventional; set conventional;
			if _TYPE_="PARMS" then call symput("convBETA",&exp);
			if _TYPE_="COV" then call symput("convVAR",&exp);
		run;	

	* merge the two parameter datasets;
	    data combineb; 
				merge sys_parms sstemp;

		* initialize the random number generator seeds;
				retain zns1 zns2; 
				if zns1 eq . then do; zns1 = -1; zns2 = -1; end;
			* create a distribution to represent result without sensitivity analysis;
				znor1=rannor(zns1);
				znor2=rannor(zns2);
				conv_var=&convVAR;
			* create a random error estimate using the conventional variance;	
				conv_parm = &convBETA - znor1*sqrt(conv_var);
			* create a distribution to represent approximation to bootstrapping analysis;
				boot_parm = sys_parm - znor2*sqrt(conv_var);
			* create odds ratio, systematic error;
				exp_sys = exp(sys_parm); 
			* create odds ratio, conventional analysis;
				exp_conv = exp(conv_parm);
			* create odds ratio, boot approx;
				exp_boot = exp(boot_parm);

				label 	conv_parm='Beta RANDOM ERROR'
						boot_parm='Beta SYSTEMATIC and RANDOM ERROR'
					  	sys_parm='Beta SYSTEMATIC' 
					  	znor1="Standard Normal Deviate for RANDOM ERROR"
					  	znor2="Standard Normal Deviate for BOOT APPROX"
					  	sys_var='Variance of SYSTEMATIC DATASET'
						exp_sys = 'Odds Ratio, Systematic Error'
						exp_conv = 'Odds Ratio, Conventional Analysis'
						exp_boot = 'Odds Ratio, Total Error Analysis'
						conv_var = 'Conventional Variance';
			* KEA: no intercept in phreg and no _link_ and _lnlike_ in systematic;
			* drop _TYPE_ _NAME_ _LINK_ _STATUS_ intercept _LNLIKE_;
			drop _TYPE_ _NAME_ _STATUS_;
		run;

	* combining data from multiple loops;
		proc append base=sensdir.&outset data=combineb force; 
		run;

	* read full data set into combinec;
		data combinec;	
			set sensdir.&outset;
		run;
%mend;

************************ DETERMINE CONVERGENCE ***********************************;
%macro conver(totalreps);
	* assume non-convergence;
		%global converge totrp; 
		%let converge=0;%let temp=0;

	* generate distribution information on the total error dataset;
		proc univariate data=sensdir.&outset CIPCTLDF noprint;
			var boot_parm;
			output out=univpcts pctlpre=p pctlpts=2.5,97.5 n=n;
		run;

	* create formats for convergence;
		proc format; value conv 1="Converged" 0="Did not converge"; run;

	* make the total number of iterations into a macro variable;
		data count; set univpcts; 
			keep n; call symput('n',n);
		run;

	* set the iterations data set by the exposure variable;
		proc sort data= sensdir.&outset; by boot_parm; run; 

	* create convergence parameters;
		data converge1; set sensdir.&outset;
		* add in the count data set;
			if _N_ = 1 then set count;
		* create an id variable;
			retain id; id=id+1; if id=. then id=1;
		* set the integer value for convergence for 97.5th and 2.5th percentiles;
			i1=int(n*.025)+1; i2=int(n*.975)+1;
		* determine the cumulative binomial probability for 97.5th and 2.5th percentiles;
			cumbin1=cdf('Binomial',id-1,.025,n);
			cumbin2=cdf('Binomial',id-1,.975,n);
		* calculate absolute distance from the 2.5th and 97.5th percentiles;
			dist1t=abs(i1-id); dist2t=abs(i2-id);
		* label variables;
			label n='Total number of iterations' dist1t='Distance from the 2.5th percentile' 
				  dist2t='Distance from the 97.5th percentile';
			keep cumbin1 cumbin2 n i1 i2 id dist1t dist2t n boot_parm;
		run;

	* sort data by distance from the 2.5th percentile;
		proc sort data=converge1; by dist1t cumbin1; run;
			
		data converge2; set converge1; by dist1t cumbin1; 
			retain tempdist cilow1 cihigh1;
		* keep only records with 2 distances (one above and one below);
			if first.dist1t = 1 and last.dist1t = 1 then delete;
		* calculate the difference in the cumulative binomials for equidistant points;
			if first.dist1t then do;
				tempdist = cumbin1;
				cilow1 = boot_parm;
			end;
			if last.dist1t then do;
				dist1 = cumbin1 - tempdist;
				cihigh1 = boot_parm;
				output;
			end;
		run;

	* sort by difference in binomials;
		proc sort data=converge2; by i1 dist1; run;

	* keep only records that have more than 95% coverage;
		data converge2; set converge2; by i1; where dist1 > .95;
			* keep only the first place where probability is above 95%;
				if first.i1 then output;
				keep cihigh1 cilow1;	
		run;

	* sort data by distance from the 97.5th percentile;
		proc sort data=converge1; by dist2t cumbin2; run;
			
		data converge3; set converge1; by dist2t cumbin2; 
			retain tempdist cilow2 cihigh2;
		* use only records with 2 distances;
			if first.dist2t = 1 and last.dist2t = 1 then delete;
			if first.dist2t then do;
				tempdist = cumbin2;
				cilow2 = boot_parm;
			end;
			if last.dist2t then do;
				dist2 = cumbin2 - tempdist;
				cihigh2 = boot_parm;
				output;
			end;
		run;

	* sort by difference in binomials;
		proc sort data=converge3; by i2 dist2; run;

	* keep only records that have more than 95% coverage;
		data converge3; set converge3; by i2; where dist2 > .95;
		* use only records with 2 distances;
			if first.i2 then output;
			keep cihigh2 cilow2;	
		run;

	* merge the two convergence datasets;
		data converge4; merge converge2 converge3;  
		* calculate the width of the intervals about the percentiles;
			widthlow  = cihigh1 - cilow1;
			widthhigh = cihigh2 - cilow2;
		* calculate the maximum of the two;
			width=max(widthhigh,widthlow);
		* label variables;
			label cihigh1='Upper 95% CI About 2.5th percentile'  	cihigh2='Upper 95% CI About 97.5th percentile' 
				  cilow1='Lower 95% CI About 2.5th percentile' 		cilow2='Lower 95% CI About 97.5th percentile' 
				  widthlow='Width of CI about 2.5th percentile' 	widthhigh='Width of CI about 97.5th percentile'
				  width='Max of two intervals';
		run;

	* merge the convergence dataset with the 2.5th and 97.5th percentiles to determine convergence;
		data converg; merge converge4 univpcts;
		* calculate width from the 2.5th to 97.5th percentiles;
			cipctl=p97_5-p2_5;
		* calculate 1/20th of the width from the 2.5th to 97.5th percentiles;
			cipctl20=cipctl/20;
		* create a marker of the total repetitions;
			numconv=n;
		* determine convergence if the max of intervals about percentiles is less than 1/20th the width from 
		  the 2.5th to 97.5th percentiles;
			if cipctl ne . and width ne . then converge =(cipctl20 > width);	
			else converge=0;	
		* format variables;
			format converge conv.;
		* label variables;
			label 	cipctl="Width from 2.5th to 97.5th% percentile"
					cipctl20="1/20th of width of interval about percentiles"
					width="Width of maximum interval about 2.5th and 97.5th percentiles"
					converge="Did the algorithm converge?" 
					numconv="Number of iterations used";
		* drop unnecessary variables;
			*keep cipctl cipctl20 converge numconv width ;
		* create a macro variable to indicate convergence;
			call symput("converge",converge);
		* create a macro variable to indicate total number of iterations;
			call symput("totrp",n);
		run;

	* save convergence info;
		/*proc append base=sensdir.converg_&outset data=converg force; 
		run;*/

	* if convergence asked for and reached, set j = 30,000 so looping stops;
		%if &totalreps=30000 and &converge=1 %then %let j=&totalreps;
%mend;

********************************  GRAPH  *****************************************;
%macro graph(exp);
	* use proc univariate to generate tabular results, written to data set pctls; 
		proc univariate data=combinec cipctldf (type=asymmetric) noprint; 
			var sys_parm conv_parm boot_parm;
		* create output dataset with min, max, mean, median and mode;
			output out=pctls
				min = sys_min conv_min boot_min     mean = sys_mean conv_mean boot_mean
				median = sys_med conv_med boot_med	max = sys_max conv_max boot_max
				var = sys_var conv_var boot_var     pctlpts = 2.5 97.5 
				pctlpre = sys_ conv_ boot_		 	pctlname = p025 p975;
		run;

	* determine dropped number;
		%global dropnum;%let dropnum=0;
		%let droplen=%sysfunc(open(sensdir.drop_&outset,i));
		%if &droplen > 0 %then %do;
			proc means data=sensdir.drop_&outset noprint; 
				var sens;
				output out=dropstat n=dropnum;
			run;

			data dropstat; set dropstat;
				call symput("dropnum",dropnum);
			run;
		%end;
		%if &droplen = 0 %then %let dropnum=0;
	* label the tabular results;
		data pctls; set pctls;
			label
			sys_mean	=	mean, sensitivity analysis
			conv_mean	=	mean, conventional analysis
			boot_mean   =   mean, total error analysis

			sys_var		=	variance, sensitivity analysis
			conv_var	=	variance, conventional analysis
			boot_var	=	variance, total error analysis

			sys_max		=	maximum, sensitivity analysis
			conv_max	=	maximum, conventional analysis
			boot_max	=	maximum, total error analysis

			sys_med		=	median, sensitivity analysis
			conv_med	=	median, conventional analysis
			boot_med	=	median, total error analysis

			sys_min		=	minimum, sensitivity analysis
			conv_min	=	minimum, conventional analysis
			boot_min	=	minimum, total error analysis

			sys_p025	=	2.5%, sensitivity analysis
			sys_p975	=	97.5%, sensitivity analysis

			conv_p025	=	2.5%, conventional analysis
			conv_p975	=	97.5%, conventional analysis

			boot_p025	=	2.5%, total error analysis
			boot_p975	=	97.5%, total error analysis;
		run;

	* get tabular results;
		proc means data=pctls median noprint;
			var sys_mean sys_var sys_min sys_p025 sys_med sys_p975 sys_max  
				conv_mean conv_var conv_min conv_p025 conv_med conv_p975 conv_max
				boot_mean boot_var boot_min boot_p025 boot_med boot_p975 boot_max; 
				title1 "Tabular Results";	
				output out=foostats;
		run;

		data foostats; set foostats;
			where _STAT_ ne 'N';
			if _N_ = 1 then output;
		run;

	* create output dataset;
		data new; set foostats;
			labe='Conventional Analysis';
			param=conv_med;
			param5=conv_p025;
			param95=conv_p975;
			output;

			labe='Sensitivity Analysis'; 
			param=sys_med;  
			param5=sys_p025; 
			param95=sys_p975;
			output;

			labe='Total Error Analysis';  
			param=boot_med; 
			param5=boot_p025; 
			param95=boot_p975;
			output;

			keep labe param param5 param95;
		run;

	* create data for final output table;
		data new; set new;
			exparam=exp(param);
			exparam5=exp(param5);
			exparam95=exp(param95);
			width=exparam95/exparam5;
			label  labe="Analysis" exparam = "OR median estimate" exparam5 = "(a) OR  2.5 percentile " exparam95 = "(b) OR 97.5 percentile"
					width='Width of Interval (b/a)';
		run;

		data new; set new;
			retain base;
			if _N_=1 then base=width;
			chn=(width-base)/base;
			chng=round(chn,.001)*100;
			change='         ';
			change= chng;
			label change='% change from conventional';
		run;

			%if &converge=1 %then %let ti="Converged";
			%else %let ti="Did Not Converge";
	* determine correlation of sens and spec;
		proc corr data=combinec noprint outs=corrs;
			var sens sens2 spec spec2;
		run;

		data corrs; set corrs; by _TYPE_; where _type_='CORR';
			retain senscorr speccorr;
			if _NAME_= 'sens' then senscorr=round(sens2,.0001);
			if _NAME_= 'spec' then speccorr=round(spec2,.0001);
			if last._TYPE_ then output;
			keep senscorr speccorr;
		run;

		data corrs; set corrs;
			call symput('senscorr',senscorr);
			call symput('speccorr',speccorr);
		run;

	* Calculate p-values for each;
		data combinec; set combinec;
			if sys_parm > 0 then sys_pval=1;else sys_pval=0;
			if conv_parm > 0 then conv_pval=1;else conv_pval=0;
			if boot_parm > 0 then boot_pval=1;else boot_pval=0; 
			label sys_pval='Is RANDOM ERROR estimate above 0?'
			      conv_pval='Is SYSTEMATIC ERROR estimate above 0?'
			      boot_pval='Is SYSTEMATIC and RANDOM ERROR estimate above 0?';
		run;

		proc means data=combinec noprint;
			var sys_pval conv_pval boot_pval;
			output out=outpval mean=sys_pval conv_pval boot_pval;
		run;
		
		data outpval; set outpval;
			sys_pval=round(sys_pval,.0001);
			conv_pval=round(conv_pval,.0001); 
			boot_pval=round(boot_pval,.0001);
		run;

		data outpval; set outpval;
			call symput('sys_pval',sys_pval); 
			call symput('conv_pval',conv_pval);  
			call symput('boot_pval',boot_pval); 
		run;

		data new; set new;
			if labe='Conventional Analysis' then pval=1-&sys_pval;
			if labe='Sensitivity Analysis' then pval=1-&conv_pval;
			if labe='Total Error Analysis' then pval=1-&boot_pval;
			label pval='% of observations < 1';
		run;

		proc print data=new label noobs;
			var labe exparam5 exparam exparam95 width change pval
				/STYLE=[FONT_FACE="Times New Roman" font_size=4];;
			title1 "SENSITIVITY ANALYSIS RESULTS"; title2;
			title3 "Algorithm:  " &ti "," %trim(&totrp) " repetitions in dataset: " &outset;
			title4; 
			title5 "Exposure: " &exp  ", Outcome  : " &depend " , Deleted iterations: " &dropnum; 
			title6 "Misclassified variable: "  &miscvar ", within levels of: " &mis_ind;
			title7 ;
			title8 &titles1 &titles3 color=blue "Corr: " %trim(&senscorr);
			title9 ;
			title10 &titles2 &titles4 color=blue "Corr: " %trim(&speccorr);
		run;

	quit;		
	* clear the title;
		title1; title2; title3; title4; title5; title6; title7; title8; title9; title10; 

	* create a data set for plots that do not exceed 300 to make plots work;
		data combineplots; set combinec;
			where exp_sys < 300 and exp_boot < 300 and exp_conv < 300;
		run;
	* use proc rank to calculate percentiles that will be plotted on the cumulative
 	  probability functions, outputting results to dataset combine2;
		proc rank data=combineplots out=combine2 percent ties=mean;
			ranks conv_pctl sys_pctl boot_pctl;
			var conv_parm sys_parm boot_parm;
		run;

	* label the output of proc ranks for legends on the plots;	
		data combine2; set combine2;
			label	conv_pctl 		= 'Conventional Result'
					sys_pctl 		= 'Sensitivity Analysis'
					boot_pctl 		= 'Total Error';
		run;				

	* set SAS graphing options for graphs that generate cumulative distributions;
	* all fonts set to Swiss, legend option assumes three lines on the graph;
		legend1 across=3 frame fwidth=4 label=none mode=reserve
			position=(outside bottom center) value=(F=swiss H=9 pt);
	* options for the y-axis;
		axis1 c=bl label=(A=90 F=swiss H=16 pt J=center 'Cumulative percentile') 
			major=(W=2) offset=(0,0) width=4
			order = (0 to 110 by 10) value=(F=swiss H=12 pt '0' '10%' '20%' '30%' '40%'
					'50%' '60%' '70%' '80%' '90%' '100%' ' ');
	* options for the x-axis, on the arithmetic scale;
		axis2 c=bl label=(A=0 F=swiss H=16 pt J=center 'ln(Odds Ratio)') 
			value=(F=swiss H=12 pt) major=(W=2) offset=(0,1) width=4 ;
	* options for the x-axis, on the logarithmic scale;
		axis3 c=bl label=(A=0 F=swiss H=16 pt J=center 'Odds Ratio') 
			value=(F=swiss H=12 pt) offset=(0,1) width=4 logbase=10 logstyle=expand;
	* symbol definitions to create the lines for the plots;
		symbol1 ci=bl interpol=splineps line=20 value=none width=12;
		symbol2 ci=bl interpol=splineps line=1 value=none width=12;
		symbol3 ci=bl interpol=splineps line=41 value=none width=12;
		symbol4 ci=red interpol=splineps line=1 value=none width=12;
		symbol5 ci=red interpol=splineps line=20 value=none width=12;

	* plot conventional result, systematic error and bootstrap approximation and on arithmetic scale;
		proc gplot data=combine2 uniform;
			plot sys_pctl*sys_parm boot_pctl*boot_parm conv_pctl*conv_parm/ legend=legend1 
				noframe overlay chref=gray cvref=gray vref=2.5 50 97.5 
				vaxis=axis1 haxis=axis2;
				title1 "Sensitivity Plots";
		run;

	* plot conventional result, systematic error and bootstrap approximation on logarithmic scale;
		proc gplot data=combine2 uniform;
			plot sys_pctl*exp_sys boot_pctl*exp_boot conv_pctl*exp_conv/ legend=legend1 
				noframe overlay chref=gray cvref=gray vref=2.5 50 97.5 
				vaxis=axis1 haxis=axis3;
				title1 "Log Scale Sensitivity Plots";
		run;	
		quit;
	
	* options for the y-axis;
		axis1 c=bl label=(A=90 F=swiss H=16 pt J=center 'Frequency') 
			major=(W=2) offset=(0,0) width=4 value=(F=swiss H=12 pt);
	* options for the x-axis;
		axis2 c=bl label=(A=0 F=swiss H=16 pt J=center) value=(F=swiss H=12 pt) 
			offset=(1 cm,1 cm) width=4 value=(F=swiss H=12 pt);
	* plot sens and spec actual;
		proc gchart data=combinec; 
			vbar sens sens2 spec spec2/raxis=axis1 maxis=axis2; 
			title "Actual Sensitivity and Specificity Plots";
		run; 
		quit;
%mend;

*********************** ERROR CHECKING MACRO ***********************;
%macro checkerr();
	* create global macro variable to decide termination and set to 0;
		%global quitmac ndspec ndsens;%let quitmac=0;%let ndspec=0;%let ndsens=0;

	* determine if data set and exposure and outcome variables exist;
		%let datlen=%sysfunc(open(&inset,i));
		%if &datlen=0 %then %do;
		    %put ERROR:  DATASET &inset DOES NOT EXIST;
			%let quitmac=1;%goto stoper;
		%end; 
		%let explen=%sysfunc(varnum(&datlen,&exp));
		%if &explen=0 %then %do;
		    %put ERROR:  VARIABLE &exp DOES NOT EXIST;
			%let quitmac=1;%goto stoper;
		%end;  
		%let deplen=%sysfunc(varnum(&datlen,&depend));
		%if &deplen=0 %then %do;
		    %put ERROR:  VARIABLE &depend DOES NOT EXIST;
			%let quitmac=1;%goto stoper;
		%end;  

	* if log not specified, set equal to 2;
		%if %length(&log)=0 %then %let log = 2;

	* if startover not specified, set equal to yes;
		%if %length(&startover)=0 %then %do;
		      %let startover = yes;
			  %put WARNING:  &startover NOT SPECIFIED, SETTING = yes;
		%end;    	
	      	
	* if outset not specified, set equal to tempdata;
		%if %length(&outset)=0 %then %do;
		     %let outset = tempdata;
			 %put WARNING:  NO OUTPUT DATASET SPECIFIED, SETTING = tempdata;
		%end; 

	* if no libname specified, terminate;
		%if %length(&libname)=0 %then %do;
			%put ERROR:  NO LIBNAME SPECIFIED;
			%let quitmac=1;%goto stoper;
		%end; 

	* if input dataset not given, terminate;
		%if %length(&inset)=0 %then %do;
		    %put ERROR:  NO INPUT DATASET SPECIFIED;
			%let quitmac=1;%goto stoper;
		%end; 

	* if dependent variable not given, terminate;
		%if %length(&depend)=0 %then %do;
		    %put ERROR:  NO DEPENDENT VARIABLE SPECIFIED;
			%let quitmac=1;%goto stoper;
		%end; 

	* if predictor variable not given, terminate;
		%if %length(&exp)=0 %then %do;
		    %put ERROR:  NO PREDICTOR VARIABLE SPECIFIED;
			%let quitmac=1;%goto stoper;
		%end; 

	* if no misclassification variable is specified, use exp;
		%if %length(&miscvar)=0 %then %do;
			%put WARNING:  NO MISCLASSIFICATION VARIABLE SPECIFIED, SETTING = &exp;
			%let miscvar=&exp;
		%end;
		%let mvarlen=%sysfunc(varnum(&datlen,&miscvar));
		%if &mvarlen=0 %then %do;
		    %put ERROR:  VARIABLE &miscvar DOES NOT EXIST;
			%let quitmac=1;%goto stoper;
		%end;  

	* make sure all misclassification parameters are specified;	
		* is min sensitivity specified?;
			%if %length(&sens_min)=0 %then %do;
			    %put ERROR:  NO MINIMUM SENSITIVITY FOR MISCLASSIFICATION SPECIFIED;
				%let quitmac=1;%goto stoper;
			%end;
		* is max sensitivity specified?;
			%if %length(&sens_max)=0 %then %do;
			    %put ERROR:  NO MAXIMUM SENSITIVITY FOR MISCLASSIFICATION SPECIFIED;
				%let quitmac=1;%goto stoper;
			%end;
		* is min specificity specified?;
			%if %length(&spec_min)=0 %then %do;
			    %put ERROR:  NO MINIMUM SPECIFICITY FOR MISCLASSIFICATION SPECIFIED;
				%let quitmac=1;%goto stoper;
			%end;
		* is max specificity specified?;
			%if %length(&spec_max)=0 %then %do;
			    %put ERROR:  NO MAXIMUM SPECIFICITY FOR MISCLASSIFICATION SPECIFIED;
				%let quitmac=1;%goto stoper;
			%end;
		* if indicator is not specified, use the dependent variable;
			%if %length(&mis_ind)=0 %then %do;
				%let mis_ind=&depend;
			    %put WARNING:  NO MISCLASSIFICATION INDICATOR GIVEN, SETTING = &depend;	
			%end;
			%let mivarlen=%sysfunc(varnum(&datlen,&mis_ind));
			%if &mivarlen=0 %then %do;
			    %put ERROR:  VARIABLE &mis_ind DOES NOT EXIST;
				%let quitmac=1;%goto stoper;
			%end;  
	* make sure exp is not in the indep;
		%let cov_pred=%index(%upcase(&indep),%upcase(&exp));
		%if &cov_pred > 0 %then %do;
			%put ERROR: YOU CANNOT INCLUDE THE exp VARIABLE IN THE indep STRING;
			%let quitmac=1;%goto stoper;
		%end;

	* make sure that the misc var is in predictors or is one of the main variables;
		%let mis1=%index(%upcase(&indep),%upcase(&miscvar));
		%let mis2=%index(%upcase(&depend),%upcase(&miscvar));
		%let mis3=%index(%upcase(&exp),%upcase(&miscvar));
		%let mistake=&mis1+&mis2+&mis3;
		%if &mistake=0 %then %do;
			%put ERROR: MISCLASSIFIED VARIABLE NOT INCLUDED IN indep STRING AND NOT THE depend OR exp VARIABLE;
			%let quitmac=1;%goto stoper;
		%end;

	* if no totalreps specified, find own convergence;
		%if %length(&totalreps)=0 %then %do;
			%let totalreps=0;
			%put WARNING: &totalreps NOT SPECIFIED, WILL LOOP TO CONVERGENCE;
		%end;

	* make sure the miscvar is not the miscind;
		%if &mis_ind=&miscvar %then %do;
			%put ERROR: THE miscvar = THE mis_ind VARIABLE;
			%let totalreps=0;
			%let quitmac=1;%goto stoper;
		%end;

	* if non-differential sens, set sens to be same;
		%if %length(&sens_min_d)=0 and %length(&sens_max_d)=0 %then %do;
			%let ndsens=1;
			%let sens_min_d=0;
			%let sens_max_d=0;
			%let sens_mod_d=0;
			%let sens_mod_d2=0;
		%end;

	* if non-differential spec, set sens to be same;
		%if %length(&spec_min_d)=0 and %length(&spec_max_d)=0 %then %do;
			%let ndspec=1;
			%let spec_min_d=0;
			%let spec_max_d=0;
			%let spec_mod_d=0;
			%let spec_mod_d2=0;
		%end;

	* If differential misclassification but no correlation specified, set to .8 and warn user;
		%if &ndsens ne 1 %then %do;
			%if %length(&corrsens)=0 %then %do;
				%let corrsens=.8;
				%put WARNING: DIFFERENTIAL sensitivity SPECIFIED, BUT WITHOUT A CORRELATION;
				%put WARNING: CORRELATION IS SET TO 0.8 BUT USER SHOULD TRY ALTERNATE SPECIFICATIONS;
			%end;
		%end;
		%if &ndspec ne 1 %then %do;
			%if %length(&corrspec)=0 %then %do;
				%let corrspec=.8;
				%put WARNING: DIFFERENTIAL specificity SPECIFIED, BUT WITHOUT A CORRELATION;
				%put WARNING: CORRELATION IS SET TO 0.8 BUT USER SHOULD TRY ALTERNATE SPECIFICATIONS;			%end;
		%end;

	* make sure all misclassification parameters are specified;	
		%if %length(&sens_min_d)>0 %then %do;
		* is min sensitivity specified?;
			%if %length(&sens_min_d)=0 %then %do;
			    %put ERROR:  NO MINIMUM SENSITIVITY DIFFERENTIAL FOR MISCLASSIFICATION SPECIFIED;
				%let quitmac=1;%goto stoper;
			%end;
		* is max sensitivity specified?;
			%if %length(&sens_max_d)=0 %then %do;
			    %put ERROR:  NO MAXIMUM SENSITIVITY DIFFERENTIAL FOR MISCLASSIFICATION SPECIFIED;
				%let quitmac=1;%goto stoper;
			%end;
		* make sure that the min is not greater than the max, etc;
			%if &sens_min_d>&sens_max_d %then %do;
				%put ERROR: Sensitivity DIFFERENTIAL PARAMETERS ORDER WRONG;
				%let quitmac=1;%goto stoper;
			%end; 
		%end;

		%if %length(&spec_min_d)>0 %then %do;
		* is min specificity specified?;
			%if %length(&spec_min_d)=0 %then %do;
			    %put ERROR:  NO MINIMUM SPECIFICITY DIFFERENTIAL FOR MISCLASSIFICATION SPECIFIED;
				%let quitmac=1;%goto stoper;
			%end;
		* is max specificity specified?;
			%if %length(&spec_max_d)=0 %then %do;
			    %put ERROR:  NO MAXIMUM SPECIFICITY DIFFERENTIAL FOR MISCLASSIFICATION SPECIFIED;
				%let quitmac=1;%goto stoper;
			%end;

		* make sure that the min is not greater than the max, etc;
			%if &spec_min_d>&spec_max_d %then %do;
				%put ERROR: Specificity DIFFERENTIAL PARAMETERS ORDER WRONG;
				%let quitmac=1;%goto stoper;
			%end; 
		%end;  

		%if %length(&sens_mod)=0 and %length(&sens_mod2)=0 %then %do;
			%let sens_mod=&sens_min;
			%let sens_mod2=&sens_max;
			%put WARNING: Sensitivity Uniform Distribution;
		%end;
		%else %if %length(&sens_mod2)=0 %then %do;
			%let sens_mod2=&sens_mod;
			%put WARNING: Sensitivity Triangular Distribution;
		%end;
		%if %length(&spec_mod)=0 and %length(&spec_mod2)=0 %then %do;
			%let spec_mod=&spec_min;
			%let spec_mod2=&spec_max;
			%put WARNING: Specificity Uniform Distribution;
		%end;
		%else %if %length(&spec_mod2)=0 %then %do;
			%let spec_mod2=&spec_mod;
			%put WARNING: Specificity Triangular Distribution;
		%end;

		%if %length(&sens_min_d) and %length(&sens_mod_d)=0 and %length(&sens_mod_d2)=0 %then %do;
			%let sens_mod_d=&sens_min_d;
			%let sens_mod_d2=&sens_max_d;
			%put WARNING: Sensitivity2 Uniform Distribution;
		%end;
		%else %if %length(&sens_min_d) and %length(&sens_mod_d2)=0 %then %do;
			%let sens_mod_d2=&sens_mod_d;
			%put WARNING: Sensitivity2 Triangular Distribution;
		%end;
		%if %length(&spec_min_d) and %length(&spec_mod_d)=0 and %length(&spec_mod_d2)=0 %then %do;
			%let spec_mod_d=&spec_min_d;
			%let spec_mod_d2=&spec_max_d;
			%put WARNING: Specificity2 Uniform Distribution;
		%end;
		%else %if %length(&spec_min_d)  and %length(&spec_mod_d2)=0 %then %do;
			%let spec_mod_d2=&spec_mod_d;
			%put WARNING: Specificity2 Triangular Distribution;
		%end;
		* make sure that the min is not greater than the max, etc;
			%if &sens_min > &sens_max or &sens_min > &sens_mod or &sens_min > &sens_mod2
			or &sens_mod > &sens_mod2 or &sens_mod > &sens_max or &sens_mod2 > &sens_max %then %do;
				%put ERROR: Sensitivity PARAMETERS ORDER WRONG;
				%let quitmac=1;%goto stoper;
			%end; 
		* make sure that the min is not greater than the max, etc;
			%if &spec_min > &spec_max or &spec_min > &spec_mod or &spec_min > &spec_mod2
			or &spec_mod > &spec_mod2 or &spec_mod > &spec_max or &spec_mod2 > &spec_max %then %do;
				%put ERROR: Specificity PARAMETERS ORDER WRONG;
				%let quitmac=1;%goto stoper;
			%end; 
		* make sure that the min is not greater than the max, etc;
			%if &sens_min_d > &sens_max_d or &sens_min_d > &sens_mod_d or &sens_min_d > &sens_mod_d2
			or &sens_mod_d > &sens_mod_d2 or &sens_mod_d > &sens_max_d or &sens_mod_d2 > &sens_max_d %then %do;
				%put ERROR: Sensitivity2 PARAMETERS ORDER WRONG;
				%let quitmac=1;%goto stoper;
			%end; 
		* make sure that the min is not greater than the max, etc;
			%if &spec_min_d > &spec_max_d or &spec_min_d > &spec_mod_d or &spec_min_d > &spec_mod_d2
			or &spec_mod_d > &spec_mod_d2 or &spec_mod_d > &spec_max_d or &spec_mod_d2 > &spec_max_d %then %do;
				%put ERROR: Specificity2 PARAMETERS ORDER WRONG;
				%let quitmac=1;%goto stoper;
			%end; 
		proc printto log=dumfl; 
		run;

	* Check that &miscvar is a 0/1 dummy variable;
		proc means noprint data=&inset noprint;
			var &miscvar;
			output out=test2set min=min max=max;
		run;

		data test2set; set test2set;
			if min or max ne 1 then error=1;
			else error=0;
			call symput("err",error);
		run;			

	* Check that &mis_ind is a 0/1 dummy variable;
		proc means noprint data=&inset noprint;
			var &mis_ind;
			output out=testset min=min max=max;
		run;

		data testset; set testset;
			if min or max ne 1 then error=1;
			else error=0;
			call symput("er",error);
		run;			

		proc printto; run;

		%if &err=1 %then %do;
			%put &miscvar not a 0/1 dummy variable;
			%let quitmac=1;%goto stoper;
		%end;			

		%if &er=1 %then %do;
			%put &mis_ind not a 0/1 dummy variable;
			%let quitmac=1;%goto stoper;
		%end;
		%stoper:;
%mend;


