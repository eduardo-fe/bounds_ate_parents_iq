/*
Youth intelligence

Solutions for the raven test are:

a11 4; b12 5;c4	8;c12 2;d7  5;d12 6;e1. 7;e7. 1;

2112 kids (84%) have full ravens; another 137 have 8 ravens (5.74%) and 
another 70 have 7 ravens; the latter have worse overall scores than the main
group; we will keep these. 

*/

program drop standardize_var
program define standardize_var
    // Declare the syntax for the program
    syntax varlist(min=1 max=1) [if] [in], BYVAR(varname) NEWVAR(string)

    // Parse the input
    
	sum `varlist'
	sum `byvar'
	tempvar y
	
	gen `y'= .
	sum `y'
	
	levelsof `byvar', local(levels)
	foreach level of local levels {
		summarize `varlist' if `byvar' == `level'
		local mu = r(mean)
        local sd = r(sd)

        quietly replace `y' = (`varlist' - `mu') / `sd' if `byvar' == `level'
        
	} 
	
	gen `newvar' = `y'
	
	
end

* Create a dataset with parental education levels
local datadir "/Users/user/Documents/datasets/UKDA-6614-stata/stata/stata13_se/ukhls/"

clear
tempfile master
save `master', emptyok

* Define the list of waves
local waves a b c d e f g h i j k l

* Loop through each wave
foreach wave in `waves' {
    use pidp `wave'_qfhigh using `datadir'`wave'_indresp.dta, clear
    * Rename the variable to keep track of the wave
    rename `wave'_qfhigh qfhigh
    * Append to the master dataset
    append using `master'
    save `master', replace
}

* Load the final dataset
use `master', clear
bysort pidp: egen maxed = max(qfhigh)
codebook qfhigh
label values maxed  l_qfhigh
duplicates drop pidp , force

save "`datadir'parents_education.dta", replace




local datadir "/Users/user/Documents/datasets/UKDA-6614-stata/stata/stata13_se/ukhls/"
use "`datadir'j_youth.dta", clear

* Keep relevant variables and merge with parent-child relationship data
keep pidp pid j_hidp j_pno 
merge 1:m pid pidp using "`datadir'j_egoalt.dta"
keep if _merge==3
keep if j_rel_dv==13

* Rename variables for clarity
rename pidp pidp_child
rename apidp pidp 
rename j_asex parent_sex
keep pidp_child j_hidp j_pno pidp parent_sex 

* Drop duplicate records and save parent data
duplicates drop pidp, force
merge 1:1 pidp using "`datadir'parents_education.dta", keepusing(maxed)
keep if _merge ==3 
drop _merge

save "`datadir'parents_wave10.dta", replace

* Merge with parent cognitive ability data
merge 1:1 pidp using "`datadir'c_indresp.dta", keepusing(c_cgsrmem_dv ///
	c_cgwri_dv c_cgwrd_dv c_cgs7cs_dv c_cgs7ca_dv ///
	c_cgns1sc6_dv c_cgns1sc10_dv c_cgns2sc6_dv c_cgns2sc10_dv ///
	c_cgvfc_dv c_cgvfw_dv c_big5a_dv c_big5c_dv c_big5e_dv ///
	c_big5n_dv c_big5o_dv)
keep if _merge ==3
drop _merge

* Create parent identifiers
rename c_* *
gen dad = "_dad" if parent_sex==1
replace dad = "_mum" if parent_sex==2
drop parent_sex j_pno

* Merge with parent income and job status data
merge 1:1 pidp using "`datadir'j_indresp.dta", ///
	keepusing(j_fimngrs_dv j_fimnnet_dv j_jbstat j_age_dv)

keep if _merge==3
drop _merge

* Handle missing values and variable labels
local vars big5a_dv big5c_dv big5e_dv big5n_dv big5o_dv ///
            cgsrmem_dv cgwri_dv cgs7cs_dv cgwrd_dv cgs7ca_dv ///
            cgns1sc6_dv cgns1sc10_dv cgns2sc6_dv cgns2sc10_dv ///
            cgvfc_dv cgvfw_dv j_fimngrs_dv j_fimnnet_dv ///
			j_age_dv j_jbstat j_hidp maxed

foreach var of local vars {
    local lbl`var' : variable label `var'
	replace `var' = . if `var' <0
}

* Reshape data for mother and father separately
reshape wide   big5a_dv big5c_dv big5e_dv big5n_dv big5o_dv ///
	cgsrmem_dv cgwri_dv cgs7cs_dv cgwrd_dv cgs7ca_dv ///
	cgns1sc6_dv cgns1sc10_dv cgns2sc6_dv cgns2sc10_dv ///
	cgvfc_dv cgvfw_dv pidp j_fimngrs_dv j_fimnnet_dv j_jbstat ///
	j_hidp j_age_dv maxed, i(pidp_child) j(dad) string

* Label variables
foreach var of local vars {
    local newvar `var'_dad
    label variable `newvar' "`lbl`var''"
	local newvar `var'_mum
    label variable `newvar' "`lbl`var''"
}


* Rename child ID variable and save the cleaned dataset
rename pidp_child pidp
label variable pidp "child's id"
save "`datadir'parents_wave10.dta", replace

* Merge with youth dataset for Raven test scores
merge m:1 pidp using "`datadir'j_youth.dta"

* Define Raven test variables and correct answers
local rav_vars j_ypravena11 j_ypravenb12 j_ypravenc4 j_ypravenc12 j_ypravend7 ///
               j_ypravend12 j_ypravene1 j_ypravene5 j_ypravene7
local rav_answers 4 5 8 2 5 6 7 1 1

* Generate indicators for correct answers and calculate total Raven score
forvalues i = 1/9 {
    local var : word `i' of `rav_vars'
    local answer : word `i' of `rav_answers'
    gen rav`i'ok = `var' == `answer'
    replace rav`i'ok = . if `var' < 0
}

egen ravScore = rowtotal(rav1ok rav2ok rav3ok rav4ok rav5ok rav6ok rav7ok rav8ok rav9ok)
egen check = rownonmiss(rav1ok rav2ok rav3ok rav4ok rav5ok rav6ok rav7ok rav8ok rav9ok)

* Standardize Raven score by age
standardize_var ravScore, byvar(j_dvage) newvar(ravScore_sd)
qui sum ravScore_sd, det
gen highIQ_kid = ravScore_sd >= r(p50) & ravScore_sd !=. 
replace highIQ_kid =. if ravScore_sd==.

* Calculate average scores for each variable across parents
local vars big5a_dv big5c_dv big5e_dv big5n_dv big5o_dv ///
            cgsrmem_dv cgwri_dv cgs7cs_dv cgwrd_dv cgs7ca_dv ///
            cgns1sc6_dv cgns1sc10_dv cgns2sc6_dv cgns2sc10_dv ///
            cgvfc_dv cgvfw_dv

foreach var of local vars {

	egen `var'_mean = rowmean(`var'_mum `var'_dad)
	label variable `var'_mean "`lbl`var''"
}

* Drop merge indicators and merge with additional youth data
drop _merge
merge m:1 pidp using "`datadir'k_youth.dta", keepusing(k_ypsave)
drop _merge

* Handle negative values and calculate parental net income
replace j_fimnnet_dv_dad= . if j_fimnnet_dv_dad <0
egen parental_net_inc = rowtotal(j_fimnnet_dv_dad j_fimnnet_dv_mum)

merge m:1 j_hidp using "`datadir'j_hhresp.dta", ///
	keepusing(j_fihhmngrs_dv j_fihhmnnet1_dv j_xpmg ///
	j_xpfood1_g3 j_xpfdout_g3 j_xpaltob_g3)

keep if _merge ==3 
drop _merge 


 

* Generate cumulative distribution functions and plot
cumul  j_fihhmnnet1_dv if  cgvfc_dv_mean >=24, gen(cdf_price_24p)
cumul  j_fihhmnnet1_dv if  cgvfc_dv_mean <24, gen(cdf_price_24m)
twoway(line cdf_price_24p  j_fihhmnnet1_dv if  cgvfc_dv_mean >=24 & j_fihhmnnet1_dv<10000, sort)(line cdf_price_24m  j_fihhmnnet1_dv if  cgvfc_dv_mean <24 &j_fihhmnnet1_dv<10000, sort)


* Standarise parental congitive scores.
local parent dad mum
foreach v of local parent{
	egen mem_score_`v' = rowtotal(cgwri_dv_`v' cgwrd_dv_`v')
	standardize_var mem_score_`v' , byvar(j_age_dv_`v') newvar(memScore_sd_`v')
	standardize_var cgvfc_dv_`v' , byvar(j_age_dv_`v') newvar(words_sd_`v')
	
	standardize_var cgs7ca_dv_`v', byvar(j_age_dv_`v') newvar(sub7_sd_`v')
	
	egen try_`v' = rowtotal(cgns1sc6_dv_`v' cgns2sc6_dv_`v'), missing
	standardize_var try_`v', byvar(j_age_dv_`v') newvar(series_sd_`v')
	
	egen cogindex_`v' =rowmean(  memScore_sd_`v' words_sd_`v' sub7_sd_`v'  series_sd_`v')
	
	qui sum cogindex_`v', det
	gen highIQ_`v' = cogindex_`v' > r(p50) & cogindex_`v'!=.
	replace highIQ_`v' =. if cogindex_`v'==.

} 



* Child's outcomes.

gen saves_longterm = k_ypsave == 2
replace saves_longterm=. if k_ypsave ==4

local parentEngage ypdisbuk ypdistv ypgetbuk ypfadmus ypfadspt ypfadttr
foreach v of local parentEngage {
	
	gen aux_`v' = j_`v' == 1
	replace  aux_`v' = . if j_`v' <0 | j_`v' ==.
	
}
egen parentEngage = rowtotal(aux_ypdisbuk aux_ypdistv aux_ypgetbuk ///
	aux_ypfadmus aux_ypfadspt aux_ypfadttr)

gen does_musicArt = (j_yposclas1 ==1 |j_yposclas2 ==1 |j_yposclas8_code==1 |j_yposclas8_code==2)
gen does_sport = (j_yposclas3 ==1 |j_yposclas4 ==1 |j_yposclas8_code==3 ///
	|j_yposclas8_code==4|j_yposclas8_code==7|j_yposclas8_code==8|j_yposclas8_code==14)

//gen does_tuition = (j_yposclas1 ==5 |j_yposclas8_code==5  ) // Just 50
 
gen doesntRead =  j_ypnbuks==0

gen smokeOccasionallyRisk = j_ypsmrsk1== 4
gen smokeDailyRisk = j_ypsmrsk2 == 4
gen drinkDailyRisk = j_ypalcrsk1 ==4

gen noIntentionALevel = j_yplvsc2do != 2 
replace noIntentionALevel = . if j_yplvsc2do <0 
gen does_volunteer = j_ypfvolunt !=6 
replace does_volunteer = . if j_ypfvolunt<0

gen doesSomething = (does_musicArt ==1 | ///
	does_sport==1 | /// 
	does_volunteer==1 )
	
* Household Chars
sum j_fihhmngrs_dv, det
gen wealthy =j_fihhmngrs_dv >= r(p50)

gen parentDisengaged = parentEngage ==0
gen intWxD = wealthy * parentDise

* Analysis

* Define the outcome variables and covariates

local outcomedir "/Users/user/Dropbox/Econometrics/youthIntelligence/outcomes/"
local outcomes does_musicArt does_sport doesntRead smokeDaily ///
    drinkDaily noIntention does_volunteer saves
local x1 j_dvage
local x2 j_dvage wealthy
local x3 j_dvage parentDi
local x4 j_dvage wealthy parentDi intWxD
local x5 j_dvage wealthy parentDi intWxD cgvfc_dv_mean

 
* Create subdatasets

preserve
keep if parentDisengage ==1
save "/Users/user/Dropbox/Econometrics/youthIntelligence/code/smaller_datasets/data_disengage_1.dta", replace
restore
preserve
keep if parentDisengaged==0
save "/Users/user/Dropbox/Econometrics/youthIntelligence/code/smaller_datasets/data_disengage_0.dta", replace
restore

preserve
keep if wealthy==0
save "/Users/user/Dropbox/Econometrics/youthIntelligence/code/smaller_datasets/data_wealthy_0.dta", replace
restore

preserve
keep if wealthy==1
save "/Users/user/Dropbox/Econometrics/youthIntelligence/code/smaller_datasets/data_wealthy_1.dta", replace
restore

tab j_ypsex, gen(youth_sex)
tab j_sex_dv, gen(prnt_sex)


local outputdir "/Users/user/Dropbox/Econometrics/youthIntelligence/outcomes/"




local outputdir "/Users/user/Dropbox/Econometrics/youthIntelligence/outcomes/"

estpost summarize mem_score_dad cgvfc_dv_dad cgs7ca_dv_dad try_dad cogindex_dad highIQ_dad ///
mem_score_mum cgvfc_dv_mum cgs7ca_dv_mum try_mum cogindex_mum highIQ_mum ///
highIQ_kid ///
parentEngage ///
j_fihhmngrs_dv ravScore j_age_dv_mum j_age_dv_dad j_dvage youth_sex2 prnt_sex2

esttab using "`outputdir'summary_stats2.tex" , cells("mean sd min max count") replace
 

	 
local outcomedir "/Users/user/Dropbox/Econometrics/youthIntelligence/outcomes/"
local outcomes does_musicArt does_sport doesntRead smokeDaily ///
    drinkDaily noIntention does_volunteer saves
local x1 j_dvage
local x2 j_dvage wealthy
local x3 j_dvage parentDi
local x4 j_dvage wealthy parentDi intWxD
local x5 j_dvage wealthy parentDi intWxD cgvfc_dv_mean

* Loop over each outcome variable
foreach y of local outcomes {
    * Initialize a counter for storing estimates
    local counter = 1
    
    * Loop over the different sets of covariates
    forvalues i = 1/5 {
        * Run the probit regression with the specified covariates
        probit `y' ravScore_sd `x`i'', vce(robust)
        
        * Compute the marginal effects
        mfx
        
        * Store the marginal effects for each model
        eststo ME_`y'_`counter'
        
        * Increment the counter
        local counter = `counter' + 1
    }

    * Export the stored estimates for each outcome variable
	estout  ME_`y'_* using "`outcomedir'out`y'.tex", margin style(tex) ///
		cells(b(star fmt(3)) se) stats(N) replace 
 
}

eststo clear
 
local outcomedir "/Users/user/Dropbox/Econometrics/youthIntelligence/outcomes/"
local outcomes  j_yposclas1
local x1 wealthy
local x2 wealthy parentDisengaged
local x3 wealthy parentDisengaged highIQ_kid

forvalues j = 1/8 {
	
	forvalues k = 1/3 {
	
		probit j_yposclas`j' `x`k'' , cluster(j_hidp)
		mfx
        eststo ME_j_yposclas`j'_reg`k'
	}
	estout  ME_j_yposclas`j'_reg* using "`outcomedir'yposclas`j'.tex", margin style(tex) ///
		cells(b(star fmt(3)) se) stats(N) replace 
	
}
 * Export the stored estimates for each outcome variable
	
 
probit j_yposclas1  wealthy parentEngage highIQ_dad highIQ_mum j_age_dv_dad j_age_dv_mum  big5a_dv_dad big5c_dv_dad big5e_dv_dad big5n_dv_dad big5o_dv_dad big5a_dv_mum big5c_dv_mum big5e_dv_mum big5n_dv_mum big5o_dv_mum if j_yposclas4!=. , vce(robust)
