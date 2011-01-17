*set the library with the iterations;
	libname out 'c:\';
*set the SAS graphing options for the graphs that generate the plots;
*all fonts are set to Century, global search and replace centx with swiss for a font without sariffs;
*options for the axes 
	(prefix 1 in the axis number denotes y-axis, 
	 prefix 5 in the axis number denotes x-axis);
*classification densities plot;
	axis11 c=bl label=(A=90 F=centx H=12 pt J=center 'Probability Density') major=(H=5 pt N=4) 
		minor=none value=(A=0 F=centx H=10 pt J=right)
		order=(0 to 40 by 10);
	axis51 c=bl label=(A=0 F=centx H=12 pt J=center 'Classification proportion') 
		major=(H=5 pt) minor=none value=(A=-300 R=360 F=centx H=10 pt J=left);
*classification results plot;
	axis12 c=bl label=(A=90 F=centx H=12 pt J=center 'Cases') 
		major=(H=5 pt) minor=none value=(A=-300 R=360 F=centx H=10 pt J=left)
		order=(0 to 1 by 0.1);
	axis52 c=bl label=(A=0 F=centx H=12 pt J=center 'Controls') 
		major=(H=5 pt) minor=none value=(A=-300 R=360 F=centx H=10 pt J=left)
		order=(0 to 1 by 0.1);
*participation proportions plot;
	axis13 c=bl label=(A=90 F=centx H=12 pt J=center 'Cumulative percent of iterations') 
		major=(H=5 pt) minor=none value=(A=-300 R=360 F=centx H=10 pt J=left)
		order=(0 to 100 by 10);
	axis53 c=bl label=(A=0 F=centx H=12 pt J=center 'Selection proportion') 
		major=(H=5 pt) minor=none value=(A=-300 R=360 F=centx H=10 pt J=left)
		order=(0 to 1 by 0.1);
*confounder associations plot;
	axis14 c=bl label=(A=90 F=centx H=12 pt J=center 'Cumulative percent of iterations') 
		major=(H=5 pt) minor=none value=(A=-300 R=360 F=centx H=10 pt J=left)
		order=(0 to 100 by 10);
	axis54 c=bl label=(A=0 F=centx H=12 pt J=center 'Odds ratio') 
		major=(H=5 pt) minor=none value=(A=-300 R=360 F=centx H=10 pt J=left)
		order=(0.2 to 1.4 by 0.2);
*results histograms;
	axis15 c=bl label=(A=90 F=centx H=12 pt J=center 'Percent') major=(H=5 pt N=4) 
		minor=none value=(A=0 F=centx H=10 pt J=right)
		order=(0 to 25 by 5);
	axis55 c=bl label=(A=0 F=centx H=12 pt J=center) 
		minor=none value=(A=-300 R=360 F=centx H=10 pt J=left '0.55'
		'0.58' '0.61' '0.64' '0.67' '0.70' '0.74' '0.78' '0.82' '0.86'
		'0.90' '0.95' '1.00' '1.05' '1.11' '1.16' '1.22' '1.28' '1.35'
		'1.42' '1.49' '1.57' '1.65' '1.73' '1.82');	
*histogram pattern option;
		pattern1 c=bl v=s;
*options for the plot lines;
	*classification plot lines;
		symbol1 ci=gr interpol=join line=24 value=none width=12;
		symbol2 ci=bl interpol=join line=3 value=none width=12;
		symbol3 ci=gr interpol=join line=15 value=none width=12;
		symbol4 ci=bl interpol=join line=1 value=none width=12;
	*classification results lines;
		symbol5 cv=bl interpol=none value=circle h=.1 cm;
		symbol6 cv=bl interpol=none value=dot h=.1 cm;
		symbol7 ci=bl interpol=join line=1 value=none width=5;
	*selection plot lines;
		symbol8 ci=gr interpol=splineps line=24 value=none width=8;
		symbol9 ci=gr interpol=splineps line=3 value=none width=8;
		symbol10 ci=bl interpol=splineps line=15 value=none width=8;
		symbol11 ci=bl interpol=splineps line=1 value=none width=8;
	*confounder association plot lines;
		symbol12 ci=bl interpol=splineps line=15 value=none width=12;
		symbol13 ci=bl interpol=splineps line=3 value=none width=12;
		symbol14 ci=gr interpol=splineps line=1 value=none width=12;
*options for the legends;
	*classification legend;
		legend1 across=2 down=2 frame fwidth=4 mode=reserve label=none
			position=(outside bottom center) value=(F=centx H=8 pt);
	*classification results legend;
		legend2 across=2 down=2 frame fwidth=4 mode=reserve label=none
			position=(outside bottom center) value=(F=centx H=8 pt);
	*selection proportion legend;
		legend3 across=2 down=2 frame fwidth=4 mode=reserve label=none
			position=(outside bottom center) value=(F=centx H=8 pt);
	*confounder association legend;
		legend4 across=2 down=2 frame fwidth=4 mode=reserve label=none
			position=(outside bottom center) value=(F=centx H=8 pt);
*classification distributions;
data classifications;
*set parameters for case sensitivity trapezoidal;
	casens_min=0.45;
	casens_lmode=0.5;
	casens_umode=0.6;
	casens_max=0.65;
	casensp=2/(casens_max+casens_umode-casens_lmode-casens_min);
*set parameters for case specificity trapezoidal;
	caspec_min=0.95;
	caspec_lmode=0.97;
	caspec_umode=0.99;
	caspec_max=1.0;
	caspecp=2/(caspec_max+caspec_umode-caspec_lmode-caspec_min);
*set parameters for control sensitivity trapezoidal (a little lower than case);
	cosens_min=0.4;
	cosens_lmode=0.48;
	cosens_umode=0.58;
	cosens_max=0.63;
	cosensp=2/(cosens_max+cosens_umode-cosens_lmode-cosens_min);
*set parameters for control specificity trapezoidal (a little higher than case);
	cospec_min=0.96;
	cospec_lmode=0.98;
	cospec_umode=0.99;
	cospec_max=1.0;
	cospecp=2/(cospec_max+cospec_umode-cospec_lmode-cospec_min);
*label the lines for the legend;
	label
		trapezoidal_1 = 'Case sensitivity'
		trapezoidal_2 = 'Control sensitivity'
		trapezoidal_3 = 'Case specificity'
		trapezoidal_4 = 'Control specificity'
		;
*step through the range 0 to 1 and calculate densities;
do i=0 to 1 by 0.00001;
		if i lt casens_min or i gt casens_max then trapezoidal_1=0;
			else if casens_min le i lt casens_lmode then trapezoidal_1=casensp*(i-casens_min)/(casens_lmode-casens_min);
			else if casens_lmode le i lt casens_umode then trapezoidal_1=casensp;
			else if casens_umode le i lt casens_max then trapezoidal_1=casensp*(casens_max-i)/(casens_max-casens_umode);
		if i lt cosens_min or i gt cosens_max then trapezoidal_2=0;
			else if cosens_min le i lt cosens_lmode then trapezoidal_2=cosensp*(i-cosens_min)/(cosens_lmode-cosens_min);
			else if cosens_lmode le i lt cosens_umode then trapezoidal_2=cosensp;
			else if cosens_umode le i lt cosens_max then trapezoidal_2=cosensp*(cosens_max-i)/(cosens_max-cosens_umode);
		if i lt caspec_min or i gt caspec_max then trapezoidal_3=0;
			else if caspec_min le i lt caspec_lmode then trapezoidal_3=caspecp*(i-caspec_min)/(caspec_lmode-caspec_min);
			else if caspec_lmode le i lt caspec_umode then trapezoidal_3=caspecp;
			else if caspec_umode le i lt caspec_max then trapezoidal_3=caspecp*(caspec_max-i)/(caspec_max-caspec_umode);
		if i lt cospec_min or i gt cospec_max then trapezoidal_4=0;
			else if cospec_min le i lt cospec_lmode then trapezoidal_4=cospecp*(i-cospec_min)/(cospec_lmode-cospec_min);
			else if cospec_lmode le i lt cospec_umode then trapezoidal_4=cospecp;
			else if cospec_umode le i lt cospec_max then trapezoidal_4=cospecp*(cospec_max-i)/(cospec_max-cospec_umode);
		output;
	end;
	run;
*plot the density graphs;
*classification densities graph;
	title h=2 c=bl f=centx j=left 'A';
	proc gplot data=classifications uniform;
			plot trapezoidal_1*i=1 trapezoidal_2*i=2
			trapezoidal_3*i=3 trapezoidal_4*i=4 / vaxis=axis11 haxis=axis51 noframe overlay legend=legend1;
		run;
*classification results graph;
*select only the first two hundred iterations because 100,000 points obscure the results;
	title h=2 c=bl f=centx j=left 'B';
	data classplot; set out.multiplefin (where=(count1 le 200));
		label casens = 'Sensitivity'
			  caspec = 'Specificity'
			  id1 = 'Non-differential'
			  id2 = 'Non-differential'
			;
		id1=count1/200;
		id2=id1;
	run;

	proc gplot data=classplot uniform;
			plot casens*cosens=5 caspec*cospec=6 id1*id2=7
			/ vaxis=axis12 haxis=axis52 noframe overlay legend=legend2;
		run;
*selection proportions plot;
*use proc rank to calculate the percentiles that will be plotted on the cumulative
 probability functions, outputting results to dataset selection;
	proc rank data=out.multiplefin out=selection percent ties=mean;
		ranks rcapartAD rcapartnonAD rcopartAD rcopartnonAD;
		var capartAD capartnonAD copartAD copartnonAD;
	run;
*add labels;
	data selection; set selection;
		label
			rcapartAD = 'Cases, AD+'
			rcapartnonAD = 'Cases, AD-'
			rcopartAD = 'Controls, AD+'
			rcopartnonAD = 'Controls, AD-'
		;
	run;
*selection proportions plot;
	title h=2 c=bl f=centx j=left 'C';
	proc gplot data=selection uniform;
		plot  rcapartAD*capartAD=8 rcapartnonAD*capartnonAD=9 rcopartAD*copartAd=10 rcopartnonAD*copartnonAD=11
			/ noframe vaxis=axis13 haxis=axis53 overlay legend=legend3;
	run;		
*unmeasured confounder plot;
	data confounder; set out.multiplefin;
*calculate association between exposure and confounder;
		ORec = (exprevAD/(1-exprevAD))/(exprevnonAD/(1-exprevnonAD));
/*calculate bound on confounding;
		ORconf = 1/exprevnonAD;
			if abs(log(ORec)) lt abs(log(ORconf)) then
				ORconf = ORec;
			if abs(log(exbc)) lt abs(log(ORconf)) then
				ORconf = exbc;
			if abs(log((ORec/((1-exprevnonAD)+ORec*exprevnonAD)))) lt abs(log(ORconf)) then
				ORconf = (ORec/((1-exprevnonAD)+ORec*exprevnonAD));
			if abs(log((exbc/((1-exprevnonAD)+exbc*exprevnonAD)))) lt abs(log(ORconf)) then
				ORconf = (exbc/((1-exprevnonAD)+exbc*exprevnonAD));
*/
		ORconf = 1.22/ORc;
	run;
*use proc rank to calculate the percentiles that will be plotted on the cumulative
 probability functions, outputting results to dataset selection;
	proc rank data=confounder out=rconfounder percent ties=mean;
		ranks rORec rexbc rORconf;
		var ORec exbc ORconf;
	run;
*add labels;
	data rconfounder; set rconfounder;
		label
			rORec = 'AD-exercise OR'
			rexbc = 'BC-exercise OR'
			rORconf =  'Bound on confounding'
		;
	run;
*confounder OR plot;
	title h=2 c=bl f=centx j=left 'D';
	proc gplot data=rconfounder uniform;
		plot  rORec*ORec=12 rexbc*exbc=13 rORconf*ORconf=14
			/ noframe vaxis=axis14 haxis=axis54 overlay legend=legend4;
	run;		
*plot the results histograms;
	*take ln of all plotted OR;
			data histoplots; set out.multiplefin;
				lnORr = log(ORr);
				lnORm = log(ORm);
				lnORs = log(ORs);
				lnORc = log(ORc);
				lnORmsc = log(ORmsc);
				lnORmscr = log(ORmscr);
				label
					lnORr = 'Conventional OR'
					lnORm = 'Misclassification analysis OR'
					lnORs = 'Selection bias analysis OR'
					lnORc = 'Unmeasured confounder analysis OR'
					lnORmsc = 'All three bias analyses OR'
					lnORmscr = 'All three bias analysis and random error OR'
				;
			run;
	*plot the sampling error only histogram;
	title h=2 c=bl f=centx j=left 'A';
	proc gchart data=histoplots;
		vbar lnORr / noframe space=0 type=percent raxis=axis15 maxis=axis55 
			coutline=w woutline=3 midpoints=(-0.6 to 0.6 by 0.05);
		run;
	*plot the misclassification analysis histogram;
	title h=2 c=bl f=centx j=left 'B';
	proc gchart data=histoplots;
		vbar lnORm / noframe space=0 type=percent raxis=axis15 maxis=axis55 
			coutline=w woutline=3 midpoints=(-0.6 to 0.6 by 0.05);
		run;
	*plot the selection bias analysis histogram;
	title h=2 c=bl f=centx j=left 'C';
	proc gchart data=histoplots;
		vbar lnORs / noframe space=0 type=percent raxis=axis15 maxis=axis55 
			coutline=w woutline=3 midpoints=(-0.6 to 0.6 by 0.05);
		run;
	*plot the unmeasured confounder bias analysis histogram;
	title h=2 c=bl f=centx j=left 'D';
	proc gchart data=histoplots;
		vbar lnORc / noframe space=0 type=percent raxis=axis15 maxis=axis55 
			coutline=w woutline=3 midpoints=(-0.6 to 0.6 by 0.05);
		run;
	*plot the three bias combined analysis histogram;
	title h=2 c=bl f=centx j=left 'E';
	proc gchart data=histoplots;
		vbar lnORmsc / noframe space=0 type=percent raxis=axis15 maxis=axis55 
			coutline=w woutline=3 midpoints=(-0.6 to 0.6 by 0.05);
		run;
	*plot the three bias combined plus sampling error histogram;
	title h=2 c=bl f=centx j=left 'F';
	proc gchart data=histoplots;
		vbar lnORmscr / noframe space=0 type=percent raxis=axis15 maxis=axis55 
			coutline=w woutline=3 midpoints=(-0.6 to 0.6 by 0.05);
		run;
	quit;
