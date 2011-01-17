*! v 1.0.0 N.Orsini 28sep2006

capture program drop episens
program  episens, rclass
version 8.2
syntax varlist(min=2 max=3) [if] [in] [fw] [ , UNConf(string) MIExp(string) Format(string) STudy(string) COrder(string) * ]

	tokenize `varlist'
	local cas `1'
	local ex `2'
	local tim `3'

	marksample use

// check display format

	if "`corder'" != "" {

	if "`unconf'" == "" | "`miexp'" == "" {
		di as err "specify both the option unconf and miexp"
		exit 198
	}
		tokenize "`corder'"
		confirm   number `1'
		tempname first 
		scalar `first' = `1'	

		if inlist(`first', 1,2) == 0 {
			di as err "the number of option corder() must be 1 or 2"
			exit 198
		}
 
	}   

// check display format

	if "`format'" == "" {
		local fmt = "%3.2f"
	}   
	else {
		local fmt = "`format'"
	} 

// check type of study

	if "`study'" == "" {
		local type = "cc"
	}   
	else {
		local type = "`study'"
	} 

	if "`tim'" != "" {
		local type = "ir"
	}   

	if inlist("`type'", "cc", "ir", "cs") != 1 {
		di in red "choose cc, ir, or cs for the option study()"
		exit 198
	} 

// case-control and cumulative incidence data

	if inlist("`type'", "cc", "cs")==1 {

		quietly { 
				tempvar WGT one
				quietly gen byte `one'=1

 			if `"`weight'"'!="" { 
				qui gen double `WGT' `exp' if `use'
				local w "[fweight=`WGT']"
			}

			safesum `one' `w' if `cas' & `ex' & `use'
			local a=r(sum)
			safesum `one' `w' if `cas' & `ex'==0 & `use'
			local b=r(sum)
			safesum `one' `w' if `cas'==0 & `ex' & `use'
			local c=r(sum)
			safesum `one' `w' if `cas'==0 & `ex'==0 & `use'
			local d=r(sum)
		}
	}

// incidence rate data

	if "`type'" == "ir" {
		quietly { 
			sum `cas' `weight' if `ex' & `use'
			local a=int(r(N)*r(mean)+.5)
			sum `tim' `weight' if `ex' & `use'
			local c = r(N)*r(mean)
			sum `cas' `weight' if `ex'==0 & `use'
			local b=int(r(N)*r(mean)+.5)
			sum `tim' `weight' if `ex'==0 & `use'
			local d= r(N)*r(mean)
		}
	}

	*di " `a'  `b'  `c'  `d'"

	confirm integer number `a'
	confirm integer number `b'
	confirm integer number `c'
	confirm integer number `d'

	if `a'<0 | `b'<0 | `c'<0 | `d'<0 { 
		di in red "negative numbers invalid"
		exit 498
	}

// Unmeasured confounding 

if "`unconf'" != "" {
	tokenize "`unconf'"
	confirm   number `1'
	confirm   number `2'
	confirm   number `3'
	tempname prz1 prz0 rrdz
	scalar `prz1' = `1'
	scalar `prz0' = `2'
	scalar `rrdz' = `3'
}   

// misclassification of the exposure

if "`miexp'" != "" {

	* get sensitivities and specifities of cases and non-cases

      tokenize "`miexp'"
	confirm   number `1'
	confirm   number `2'
	confirm   number `3'
	confirm   number `4'

	tempname seec spec seeo speo 

	scalar `seec' = `1'
	scalar `spec' = `2'
	scalar `seeo' = `3'
	scalar `speo' = `4'
}   

// get the observed relative risk
 
if "`type'" == "cc" {	
	qui cci `a'  `b'  `c'  `d'
	tempname arrdx lbarr ubarr 
	scalar `arrdx' = r(or)
	scalar `lbarr' = r(lb_or)
	scalar `ubarr' = r(ub_or)
	local effect = "Odds Ratio"
}

if "`type'" == "ir" {	
	qui iri `a'  `b'  `c'  `d'
	tempname arrdx lbarr ubarr 
	scalar `arrdx' = r(irr)
	scalar `lbarr' = r(lb_irr)
	scalar `ubarr' = r(ub_irr)
	local effect = "Rate Ratio"

} 

if "`type'" == "cs" {	
	qui csi `a'  `b'  `c'  `d'
	tempname arrdx lbarr ubarr n1 n0
	scalar `arrdx' = r(rr)
	scalar `lbarr' = r(lb_rr)
	scalar `ubarr' = r(ub_rr)
	scalar `n1' = `a' + `c' 
	scalar `n0' = `b' + `d' 
	local effect = "Risk Ratio"
}

// Display the observed relative risk and 95% CI

  	di _col(1) _n as text "Observed `effect' [95% Conf. Interval]= " `fmt' as res `arrdx' ///
          in g " [" in y `fmt' `lbarr' in gr ", " in y `fmt' `ubarr' in gr "]" 

// Sensitivity analysis of a binary unmeasured confounder

if "`unconf'" != "" & "`corder'" == "" {
		
	episens_unc `a'  `b'  `c'  `d' , prz1(`prz1') prz0(`prz0') rrdz(`rrdz') arrdx(`arrdx') type(`type')
	
	di _col(1) as text _n "Unmeasured confounding"
  	di _col(4) as text "External adjusted `effect' = " `fmt' as res r(rrdx_unc)
	di _col(4) as text "Percent bias = " %3.2f as res r(bias_unc) "%"
	
		// Saved results

		return scalar rrdx_unc = r(rrdx_unc)
		return scalar bias_unc = r(bias_unc)
		return scalar arrdx = `arrdx'
}

// Sensitivity analysis for misclassification of the exposure

if "`miexp'" != "" & "`corder'" == "" {

	episens_mie `a'  `b'  `c'  `d' , seec(`seec') spec(`spec') seeo(`seeo') speo(`speo') arrdx(`arrdx') type(`type')
	
	di _col(1) as text _n "Misclassification of the exposure"
  	di _col(4) as text "External adjusted `effect' = " `fmt' as res r(rrdx_mie)
	di _col(4) as text "Percent bias = " %3.2f as res r(bias_mie) "%"

		// Saved results

		return scalar rrdx_mie = r(rrdx_mie)
		return scalar bias_mie = r(bias_mie)
		if "`unconf'" == "" return scalar arrdx = `arrdx'
}

// Combined corrections - unmeasured confounder and misclassification of the exposure 

if "`corder'" != "" {

	if `first' == 1 { 
		episens_unc `a'  `b'  `c'  `d' , prz1(`prz1') prz0(`prz0') rrdz(`rrdz') arrdx(`arrdx') type(`type')
	
		// get the counts, round the counts, and put them into locals

		tempname a1 b1 c1 d1 

		scalar `a1' = r(a1)
		scalar `b1' = r(b1)
		scalar `c1' = r(c1)
		scalar `d1' = r(d1)

		local a1 = `a1' 
		local b1 = `b1'
		local c1 = `c1'
		local d1 = `d1'

		di as text _n "Order of the corrections:"
		di _col(1) as text  "1. Unmeasured confounding"
  		di _col(4) as text "External adjusted `effect' = " `fmt' as res r(rrdx_unc)
		di _col(4) as text "Percent bias = " %3.2f as res r(bias_unc) "%"
		return scalar rrdx_unc = r(rrdx_unc)
		return scalar bias_unc = r(bias_unc)

		episens_mie `a1' `b1'  `c1' `d1' , seec(`seec') spec(`spec') seeo(`seeo') speo(`speo') arrdx(`arrdx') type(`type')
		di _col(1) as text _n "2. Misclassification of the exposure"
  		di _col(4) as text "External adjusted `effect' = " `fmt' as res r(rrdx_mie)
		di _col(4) as text "Percent bias = " %3.2f as res r(bias_mie) "%"

		// Saved results

		return scalar rrdx_mie = r(rrdx_mie)
		return scalar bias_mie = r(bias_mie)
		return scalar arrdx = `arrdx'
	}
	else {
		episens_mie `a'  `b'  `c'  `d' , seec(`seec') spec(`spec') seeo(`seeo') speo(`speo') arrdx(`arrdx') type(`type')

		// get the counts, round the counts, and put them into locals

		tempname a1 b1 c1 d1

		scalar `a1' = r(a1)
		scalar `b1' = r(b1)
		scalar `c1' = r(c1)
		scalar `d1' = r(d1)

		local a1 = `a1' 
		local b1 = `b1'
		local c1 = `c1'
		local d1 = `d1'

*	 di  `a1' " " `b1' " "  `c1' "  "  `d1'

		di as text _n "Order of the corrections:"
		di _col(1) as text   "1. Misclassification of the exposure"
  		di _col(4) as text "External adjusted `effect' = " `fmt' as res r(rrdx_mie)
		di _col(4) as text "Percent bias = " %3.2f as res r(bias_mie) "%"
		return scalar rrdx_mie = r(rrdx_mie)
		return scalar bias_mie = r(bias_mie)

		episens_unc `a1' `b1'  `c1' `d1' , prz1(`prz1') prz0(`prz0') rrdz(`rrdz') arrdx(`arrdx') type(`type')
		di _col(1) as text _n "2. Unmeasured confounding"
  		di _col(4) as text "External adjusted `effect' = " `fmt' as res r(rrdx_unc)
		di _col(4) as text "Percent bias = " %3.2f as res r(bias_unc) "%"

		// Saved results

		return scalar rrdx_unc = r(rrdx_unc)
		return scalar bias_unc = r(bias_unc)
		return scalar arrdx = `arrdx'

	}
	
}


end
 

capture program drop episens_unc
program episens_unc, rclass
syntax anything   ,  prz1(string) prz0(string) rrdz(string) arrdx(string)  type(string)

	gettoken a 0 : 0, parse(" ,")
	gettoken b 0 : 0, parse(" ,")
	gettoken c 0 : 0, parse(" ,")
	gettoken d 0 : 0, parse(" ,")
 
	tempname  rrdx rrxz  bias id b11 b01 a11 a01 rr 
 
	 scalar `rrxz' = [(`prz1')*(1-`prz0')]/[(1-`prz1')*(`prz0')] 
	 scalar `b11' = `prz1' * `c'  
	 scalar `b01' = `prz0' * `d'  
	 scalar `a11' = (`rrdz'*`a'*`b11')/(`rrdz'*`b11' +`c'-`b11')  
	 scalar `a01' = (`rrdz'*`b'*`b01')/(`rrdz'*`b01'+`d'-`b01') 

	if "`type'" == "cc"  scalar `rrdx' = (`a11'*`b01')/(`b11'*`a01')  
	if "`type'" == "ir"  scalar `rrdx' = (`a11'/`b11')/(`a01'/`b01')  
	if "`type'" == "cs"  scalar `rrdx' = (`a11'/(`a11'+`b11'))/(`a01'/(`a01'+`b01'))
	
	scalar `bias' = (`arrdx'-`rrdx')/(`rrdx')*100  

// Saved results

 	return scalar a1 = `a11'
	return scalar b1 = `a01'
	return scalar c1 = `b11'
	return scalar d1 = `b01'
	return scalar rrdx_unc = `rrdx'
	return scalar bias_unc = `bias'
end


capture program drop episens_mie
program episens_mie, rclass
syntax anything   ,  seec(string) spec(string)  seeo(string) speo(string) arrdx(string)  type(string)

	gettoken a 0 : 0, parse(" ,")
	gettoken b 0 : 0, parse(" ,")
	gettoken c 0 : 0, parse(" ,")
	gettoken d 0 : 0, parse(" ,")

	tempname fnec fpec fneo fpeo

	scalar `fnec' = 1-`seec'
	scalar `fpec' = 1-`spec'
	scalar `fneo' = 1-`seeo'
	scalar `fpeo' = 1-`speo'

// Sensitivity analysis for missclassification of the exposure
 
	tempname b1s b0s m0 b0 b1 a1 a0 rrdx bias
 
	scalar `b1' = (`speo'*`c'-`fpeo'*`d')/(`seeo'*`speo'-`fneo'*`fpeo')
	scalar `b0' = (`c'+`d'-`b1')
	scalar `a1' = (`spec'*`a'-`fpec'*`b')/(`seec'*`spec'-`fnec'*`fpec')
	scalar `a0' = (`a'+`b')-`a1'

	if "`type'" == "cc"  scalar `rrdx' = (`a1'*`b0')/(`b1'*`a0')  
	if "`type'" == "ir"  scalar `rrdx' = (`a1'/`b1')/(`a0'/`b0')  
	if "`type'" == "cs"  scalar `rrdx' = (`a1'/(`a1'+`b1'))/(`a0'/(`a0'+`b0'))
	scalar `bias' = (`arrdx'-`rrdx')/(`rrdx')*100  

// Saved results

 	return scalar a1 = `a1'
	return scalar b1 = `a0'
	return scalar c1 = `b1'
	return scalar d1 = `b0'
	return scalar rrdx_mie = `rrdx'
	return scalar bias_mie = `bias'
end
