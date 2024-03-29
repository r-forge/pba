*set the output library name;
libname out 'c:\';
*create the basic data set;
data multiple;
	input a b c d;
	cards;
	118 832 103 884
	;
data multiple2; set multiple;
*calculate conventional standard error;
	se = sqrt(1/a+1/b+1/c+1/d);
*simulation p-value calculator;
	with=0; without=0;
*set random number seed by computer clock;
	seed = -1;
*begin iterative section, set to 100 to proof check code or 100,000 to run;
	do count1 = 1 to 100000;
*set parameters for misclassification bias analysis;
*set parameters for case sensitivity trapezoidal;
	casens_min=0.45;
	casens_lmode=0.5;
	casens_umode=0.6;
	casens_max=0.65;
*set parameters for case specificity trapezoidal;
	caspec_min=0.95;
	caspec_lmode=0.97;
	caspec_umode=0.99;
	caspec_max=1.0;
*set parameters for control sensitivity trapezoidal (a little lower than case);
	cosens_min=0.4;
	cosens_lmode=0.48;
	cosens_umode=0.58;
	cosens_max=0.63;
*set parameters for control specificity trapezoidal (a little higher than case);
	cospec_min=0.96;
	cospec_lmode=0.98;
	cospec_umode=0.99;
	cospec_max=1.0;
*set parameters for selection bias analysis;
*set AD use prevalence in controls, assume normally distributed;
	coprev=0.2; coprevsd=0.05;
*set parameters for AD users to non-users participation in cases, trapezoidal;
	capart_min=0.75; 
	capart_lmode=0.85; 
	capart_umode=0.95; 
	capart_max=1;
*set AD use prevalence in cases, assume normally distributed;
	caprev=0.25; caprevsd=0.075;
*set parameters for AD users to non-users participation in controls, trapezoidal (slightly lower than cases);
	copart_min=0.7; 
	copart_lmode=0.8; 
	copart_umode=0.9; 
	copart_max=1;
*set parameters for unmeasured confounder bias analysis;
*set parameters associating incident breast cancer with exercise, trapezoidal;
	exbc_min=0.2; 
	exbc_lmode=0.58; 
	exbc_umode=1.01; 
	exbc_max=1.24;
*set exercise prevalence in AD users;
	exprevAD=0.30; exprevADsd=0.05;
*set exercise prevalence in AD nonusers;
	exprevnonAD=0.44; exprevnonADsd=0.05;
*set the strength of correlation between random numbers;
	r=0.8;
*misclassification bias analysis;
*Find two sets of correlated random standard uniform;
*set one;
	u1=ranuni(seed); u2=ranuni(seed); u3=ranuni(seed);
	u1=log(u1/(1-u1)); u2=log(u2/(1-u2)); u3=log(u3/(1-u3));
		g1=exp(sqrt(r)*u1+sqrt(1-r)*u2)/(1+exp(sqrt(r)*u1+sqrt(1-r)*u2));
		g2=exp(sqrt(r)*u1+sqrt(1-r)*u3)/(1+exp(sqrt(r)*u1+sqrt(1-r)*u3));
*set two;
	u4=ranuni(seed); u5=ranuni(seed); u6=ranuni(seed);
	u4=log(u4/(1-u4)); u5=log(u5/(1-u5)); u6=log(u6/(1-u6));
		g3=exp(sqrt(r)*u4+sqrt(1-r)*u5)/(1+exp(sqrt(r)*u4+sqrt(1-r)*u5));
		g4=exp(sqrt(r)*u4+sqrt(1-r)*u6)/(1+exp(sqrt(r)*u4+sqrt(1-r)*u6));
*Find the case sensitivities and specificities;
		casens =  (g1*(casens_umode+casens_max-casens_min-casens_lmode) + (casens_min + casens_lmode))/2;
		if casens < casens_lmode then do; 
			casens  = casens_min + sqrt((casens_lmode-casens_min)*(2*casens - casens_min - casens_lmode));
		end;
	   	else if casens > casens_umode then do; 
			casens = casens_max - sqrt(2*(casens_max-casens_umode)*(casens-casens_umode));
		end;
		caspec =  (g3*(caspec_umode+caspec_max-caspec_min-caspec_lmode) + (caspec_min + caspec_lmode))/2;
		if caspec < caspec_lmode then do; 
			caspec  = caspec_min + sqrt((caspec_lmode-caspec_min)*(2*caspec - caspec_min - caspec_lmode));
		end;
	   	else if caspec > caspec_umode then do; 
			caspec = caspec_max - sqrt(2*(caspec_max-caspec_umode)*(caspec-caspec_umode));
		end;
*Find the control sensitivities and specificities;
		cosens =  (g2*(cosens_umode+cosens_max-cosens_min-cosens_lmode) + (cosens_min + cosens_lmode))/2;
		if cosens < cosens_lmode then do; 
			cosens  = cosens_min + sqrt((cosens_lmode-cosens_min)*(2*cosens - cosens_min - cosens_lmode));
		end;
	   	else if cosens > cosens_umode then do; 
			cosens = cosens_max - sqrt(2*(cosens_max-cosens_umode)*(cosens-cosens_umode));
		end;
		cospec =  (g4*(cospec_umode+cospec_max-cospec_min-cospec_lmode) + (cospec_min + cospec_lmode))/2;
		if cospec < cospec_lmode then do; 
			cospec  = cospec_min + sqrt((cospec_lmode-cospec_min)*(2*cospec - cospec_min - cospec_lmode));
		end;
	   	else if cospec > cospec_umode then do; 
			cospec = cospec_max - sqrt(2*(cospec_max-cospec_umode)*(cospec-cospec_umode));
		end;
*find corrected cells (under assumed sensitivities and specificities);
		afixm = (a-(1-caspec)*(a+b))/(casens+caspec-1);
		bfixm = a+b-afixm;
		cfixm = (c-(1-cospec)*(c+d))/(cosens+cospec-1);
		dfixm = c+d-cfixm;
		ORm= (afixm/bfixm) / (cfixm/dfixm);
*Selection bias analysis;
*Find two sets of correlated random standard uniform;
*set one;
	u7=ranuni(seed); u8=ranuni(seed); u9=ranuni(seed);
	u7=log(u7/(1-u7)); u8=log(u8/(1-u8)); u9=log(u9/(1-u9));
		g5=exp(sqrt(r)*u7+sqrt(1-r)*u8)/(1+exp(sqrt(r)*u7+sqrt(1-r)*u8));
		g6=exp(sqrt(r)*u7+sqrt(1-r)*u9)/(1+exp(sqrt(r)*u7+sqrt(1-r)*u9));
*set two;
	u10=ranuni(seed); u11=ranuni(seed); u12=ranuni(seed);
	u10=log(u10/(1-u10)); u11=log(u11/(1-u11)); u12=log(u12/(1-u12));
		g7=exp(sqrt(r)*u10+sqrt(1-r)*u11)/(1+exp(sqrt(r)*u10+sqrt(1-r)*u11));
		g8=exp(sqrt(r)*u10+sqrt(1-r)*u12)/(1+exp(sqrt(r)*u10+sqrt(1-r)*u12));
*find case participation;
		caprev=caprev-caprevsd*probit(g5);
		capart =  (g7*(capart_umode+capart_max-capart_min-capart_lmode) + (capart_min + capart_lmode))/2;
		if capart < capart_lmode then do; 
			capart  = capart_min + sqrt((capart_lmode-capart_min)*(2*capart - capart_min - capart_lmode));
		end;
	   	else if capart > capart_umode then do; 
			capart = capart_max - sqrt(2*(capart_max-capart_umode)*(capart-capart_umode));
		end;
		capartnonAD = 0.806/(1-caprev+caprev*capart);
		capartAD = capartnonAD*capart;
*find control participation;
		coprev=coprev-coprevsd*probit(g6);
		copart =  (g8*(copart_umode+copart_max-copart_min-copart_lmode) + (copart_min + copart_lmode))/2;
		if copart < copart_lmode then do; 
			copart  = copart_min + sqrt((copart_lmode-copart_min)*(2*copart - copart_min - copart_lmode));
		end;
	   	else if copart > copart_umode then do; 
			copart = copart_max - sqrt(2*(copart_max-copart_umode)*(copart-copart_umode));
		end;
		copartnonAD = 0.806/(1-coprev+coprev*copart);
		copartAD = copartnonAD*copart;
*find corrected cells (under assumed participation);
		afixs = a/capartAD;
		bfixs = b/capartnonAD;
		cfixs = c/copartAD;
		dfixs = d/copartnonAD;
		ORs = (afixs/bfixs) / (cfixs/dfixs);
*find corrected cells (under assumed misclassification and participation);
		afixms = afixm/capartAD;
		bfixms = bfixm/capartnonAD;
		cfixms = cfixm/copartAD;
		dfixms = dfixm/copartnonAD;
		ORms = (afixms/bfixms) / (cfixms/dfixms);
*Unmeasured confounder bias analyis;
*Find one set of correlated random standard uniform;
*set one;
	u13=ranuni(seed); u14=ranuni(seed); u15=ranuni(seed);
	u13=log(u13/(1-u13)); u14=log(u14/(1-u14)); u15=log(u15/(1-u15));
		g9=exp(sqrt(r)*u13+sqrt(1-r)*u14)/(1+exp(sqrt(r)*u13+sqrt(1-r)*u14));
		g10=exp(sqrt(r)*u13+sqrt(1-r)*u15)/(1+exp(sqrt(r)*u13+sqrt(1-r)*u15));
*find two exercise prevalence;
		exprevAD=exprevAD-exprevADsd*probit(g9);
		exprevnonAD=exprevnonAD-exprevnonADsd*probit(g10);
*find exercise-BC association;
		exbc =  (ranuni(seed)*(exbc_umode+exbc_max-exbc_min-exbc_lmode) + (exbc_min + exbc_lmode))/2;
		if exbc < exbc_lmode then do; 
			exbc  = exbc_min + sqrt((exbc_lmode-exbc_min)*(2*exbc - exbc_min - exbc_lmode));
		end;
	   	else if exbc > exbc_umode then do; 
			exbc = exbc_max - sqrt(2*(exbc_max-exbc_umode)*(exbc-exbc_umode));
		end;								
*find strata cell frequencies, confounding only;
		aex1 = (exbc*(a+c)*exprevAD*a)/(exbc*(a+c)*exprevAD+(a+c)-(a+c)*exprevAD);
		aex0 = a-aex1;
		cex1 = (a+c)*exprevAD-aex1;
		cex0 = c-cex1;
		bex1 = (exbc*(b+d)*exprevnonAD*b)/(exbc*(b+d)*exprevnonAD+(b+d)-(b+d)*exprevnonAD);
		bex0 = b-bex1;
		dex1 = (b+d)*exprevnonAD-bex1;
		dex0 = d-dex1;
		ORc = (aex1*dex1/(aex1+bex1+cex1+dex1)+aex0*dex0/(aex0+bex0+cex0+dex0)) /
			  (bex1*cex1/(aex1+bex1+cex1+dex1)+bex0*cex0/(aex0+bex0+cex0+dex0));									
*find strata cell frequencies, confounding added to misclassification and selection;
		aaex1 = (exbc*(afixms+cfixms)*exprevAD*afixms)/(exbc*(afixms+cfixms)*exprevAD+(afixms+cfixms)-(afixms+cfixms)*exprevAD);
		aaex0 = afixms-aaex1;
		ccex1 = (afixms+cfixms)*exprevAD-aaex1;
		ccex0 = cfixms-ccex1;
		bbex1 = (exbc*(bfixms+dfixms)*exprevnonAD*bfixms)/(exbc*(bfixms+dfixms)*exprevnonAD+(bfixms+dfixms)-(bfixms+dfixms)*exprevnonAD);
		bbex0 = bfixms-bbex1;
		ddex1 = (bfixms+dfixms)*exprevnonAD-bbex1;
		ddex0 = dfixms-ddex1;
		ORmsc = (aaex1*ddex1/(aaex1+bbex1+ccex1+ddex1)+aaex0*ddex0/(aaex0+bbex0+ccex0+ddex0)) /
			  (bbex1*ccex1/(aaex1+bbex1+ccex1+ddex1)+bbex0*ccex0/(aaex0+bbex0+ccex0+ddex0));
*incorporate sampling error;
		s = se*rannor(seed);
		ORr = exp(log((a/c)/(b/d))-s);
		ORmr = exp(log(ORm)-s);
		ORsr = exp(log(ORs)-s);
		ORcr = exp(log(ORc)-s);
		ORmsr = exp(log(ORms)-s);
		ORmscr = exp(log(ORmsc)-s);
*simulation p-value calculator;
		if ORmsc gt 1 then without=without+1;
		if ORmscr gt 1 then with=with+1;
		output;
	end;
put with without;
run;
*store output data set for future work and plots;
	data out.multiplefin; set multiple2; run;
*diagnostics for the bias analysis;
*confirm expected correlations;
	proc corr data=multiple2;
			var  g1 g3 g5 g7 g9  casens caspec caprev capart exprevAD;
			with g2 g4 g6 g8 g10 cosens cospec coprev copart exprevnonAD;
		run;
*examine distributions;
	proc univariate data=multiple2;
		var casens cosens caspec cospec 
		capart copart caprev coprev capartAD capartnonAD copartAD copartnonAD
		exprevAD exprevnonAD
		s;
*use proc univariate to generate median and interval for each OR;
	proc univariate data=multiple2 cipctldf (type=asymmetric) noprint;
		 var ORm ORs ORc ORms ORmsc
		 ORr ORmr ORsr ORcr ORmsr ORmscr;
		* create output dataset ;
			output out=results
				median=mdORm mdORs mdORc mdORms mdORmsc
				 mdORr mdORmr mdORsr mdORcr mdORmsr mdORmscr
				
				pctlpts = 2.5 97.5 

				pctlpre = ORm ORs ORc ORms ORmsc
		 		ORr ORmr ORsr ORcr ORmsr ORmscr

				pctlname = p025 p975;
	run;
proc print data=results; run;
