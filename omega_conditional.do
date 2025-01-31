mata:  
    real matrix eigenvalue_adjust(string scalar __v) {
        real matrix v
        v = st_matrix(__v)
        
        // Set small values in eigenvalues to zero
        v = mm_cond(v :< 1e-4, 0, v)
        
        st_matrix("__vsqrt", sqrt(v))
		 
        return(sqrt(v))
    }
end

mata:
void phat(real scalar n, real scalar __ub_pc,  real scalar __lb_pc, real scalar __ub75_pc,  real scalar __lb75_pc, real scalar __ub25_pc,  real scalar __lb25_pc)
{
	scalar gammaHat, gammaHatPlus, rho, tau, p
	
	gammaHat = __ub_pc - __lb_pc
	
	gammaHatPlus = rowmax((0, gammaHat))
	
	rho = rowmax(((__ub75_pc - __ub25_pc) , (__lb75_pc - __lb25_pc)))
	
	tau = 1/(rho*log(n))
	p = 1 - normal(tau*gammaHatPlus)*0.05
	st_numscalar("phat", p)
	
	
}
end 

mata:

void rho_adjust(real scalar N, string scalar bound, string scalar sigma, string scalar sigmaHalf, real scalar lb, real scalar c)
{
	matrix omega, omegaHalf, vShat, Z, gnu, Znu, cap,phi, test, keep, theta
	scalar k, i, normGnu, cn, kappa, R, sol, bnd
	R=10000
		
	omega = st_matrix(sigma)
	k = rows(omega)
	
	omega
	if(k==1 & abs(J(1,k,1)*omega*J(k,1,1))<=1e-4){
		st_numscalar("__bnd", st_matrix(bound))
		st_numscalar("__ci", 0)
	} 
	else{
	
	omegaHalf 	= st_matrix(sigmaHalf)
	omegaHalf
	theta = st_matrix(bound)
	
	vShat = (sqrt(diagonal(omega)/N))
	
	
	
	Z = rnormal(k, R,0,1)
	
	cap = J(1,R,.)
	
	for(i= 1; i<= k; i++) {
	
		gnu = omegaHalf[i,]; 

		normGnu = sqrt(gnu*gnu'); 

		Znu = (gnu*Z)/normGnu;
		cap = cap \ Znu
		
	}
	
	cap = cap[2..k+1,]
	
	phi = colmax(cap)

	
	cn = 1-(0.1/log(N));
	
	kappa = mm_quantile(phi',1, (cn))
	
	keep = J(k,1,.)
	
	for(i= 1; i<= k; i++) {
		if(lb !=1) {
			 
			test = theta  + kappa*vShat + 2*kappa *(vShat[i,]*J(k,1,1))
			test = colmin(test)
			test = mm_cond(theta[i,] :<= test, 1, 0)
			
			keep[i,] = test
		}
		else{
			test = theta  - kappa*vShat - 2*kappa *(vShat[i,]*J(k,1,1))
			test = colmax(test)
			test = mm_cond(theta[i,] :>= test, 1, 0)
			
			keep[i,] = test
		}
	}
	
	 
	cap = select(cap, keep)
	
	phi= colmax(cap)
	
	kappa = mm_quantile(phi',1, c)
	
	
	if(lb !=1) {
		sol = colmin(theta+ kappa*vShat)
		bnd =colmin(theta)
		 
	}
	else{
		sol = colmax(theta- kappa*vShat)
		bnd =colmax(theta)
		 
	}
	st_numscalar("__bnd", bnd)
	st_numscalar("__ci", sol)
	}
	
}
end


program drop cnpbound
program cnpbound, rclass

	syntax varname(numeric max=3) [if] [, Positive]	
	
	tempvar  Y Df Ds
	
	tokenize `varlist'
			
	gen `Y' = `1'
	gen `Df' = `2'
	gen `Ds' = `3'
	//gen `CondVar' = `4'
	
		
	marksample touse

	
	sum `Y' // Note here we compute the max in the full sample, not touse
	local ymin = r(min)
	local ymax = r(max)
	
	sum `Ds' if `Df'==0
	local pi_s1f0 = r(mean)
	local pi_s0f0 = 1-`pi_s1f0'

	sum `Df' if `Ds'==0
	local pi_f1s0 = r(mean)
	local pi_f0s0 = 1-`pi_f1s0'
	
	sum `Y' if `Df'==0 & `Ds' ==1  & `touse'
	local y_01 = r(mean)
	sum `Y' if `Df'==0 & `Ds' ==0  & `touse'
	local y_00 = r(mean)
	sum `Y' if `Df'==1  &`Ds' ==0  & `touse'
	local y_10 = r(mean)
	

	
	local lbY01_f0_nab = `ymin' *`pi_s0f0'  + `y_01'*`pi_s1f0'
	local ubY01_f0_nab = `ymax'*`pi_s0f0' + `y_01'*`pi_s1f0'
		
	local lbY01_f0_mtr = `y_00' *`pi_s0f0'  + `y_01'*`pi_s1f0'
	local ubY01_f0_mtr =  `ymax'*`pi_s0f0' + `y_01'*`pi_s1f0'
			
	local lbY01_f0_mts = `ymin' *`pi_s0f0'  + `y_01'*`pi_s1f0'
	local ubY01_f0_mts = `y_01'*`pi_s0f0' + `y_01'*`pi_s1f0'
		
	local lbY01_f0_mtrs= `y_00' * `pi_s0f0' + `y_01'*`pi_s1f0'
	local ubY01_f0_mtrs=  `y_01' * `pi_s0f0' + `y_01'*`pi_s1f0'
		
	local lbY00_f0_nab = `y_00' *`pi_s0f0'  + `ymin'*`pi_s1f0'
	local ubY00_f0_nab = `y_00' *`pi_s0f0' + `ymax'*`pi_s1f0'
		
	local lbY00_f0_mtr =`y_00' * `pi_s0f0' + `ymin'*`pi_s1f0'
	local ubY00_f0_mtr =`y_00' * `pi_s0f0' + `y_01'*`pi_s1f0'	
	
	local lbY00_f0_mts = `y_00' *`pi_s0f0'  + `y_00'*`pi_s1f0'
	local ubY00_f0_mts = `y_00' *`pi_s0f0' + `ymax'*`pi_s1f0'
	
	local lbY00_f0_mtrs=  `y_00' *`pi_s0f0'  + `y_00' * `pi_s1f0'
	local ubY00_f0_mtrs=  `y_00' *`pi_s0f0' +  `y_01' * `pi_s1f0'
	

	local assumption nab mtr mts mtrs
	foreach j of local assumption {
		
		
		local ate_f0_l_`j' =  `lbY01_f0_`j'' - ( `ubY00_f0_`j'')
		local ate_f0_u_`j' =  `ubY01_f0_`j'' - ( `lbY00_f0_`j'')
		
				
		display "---------------------------------------------"
		display ""
		display "ATE: `j' E(Y(0,1)-Y(0,0)|Z_f=0): " _column(16) "[" %9.4f `ate_f0_l_`j'' "," %9.4f `ate_f0_u_`j''  "]"	
		display ""
		matrix mate_f0_l_`j' = `ate_f0_l_`j''
		matrix mate_f0_u_`j' = `ate_f0_u_`j''
		
		return matrix mate_f0_`j'_l = mate_f0_l_`j'
		return matrix mate_f0_`j'_u = mate_f0_u_`j'
			
	}
	
		   
	local lbY10_s0_nab = `ymin' * `pi_f0s0' + `y_10'*`pi_f1s0'
	local ubY10_s0_nab = `ymax' * `pi_f0s0' + `y_10'*`pi_f1s0'
		  
	local lbY10_s0_mtr = `y_00' * `pi_f0s0' + `y_10'*`pi_f1s0'
	local ubY10_s0_mtr = `ymax' * `pi_f0s0' + `y_10'*`pi_f1s0'
	
	local lbY10_s0_mts = `ymin' * `pi_f0s0' + `y_10'*`pi_f1s0'
	local ubY10_s0_mts = `y_10' * `pi_f0s0' + `y_10'*`pi_f1s0'
	
	local lbY10_s0_mtrs = `y_00' * `pi_f0s0' + `y_10'*`pi_f1s0'
	local ubY10_s0_mtrs = `y_10' * `pi_f0s0' + `y_10'*`pi_f1s0'

	
	local lbY00_s0_nab =`y_00' *`pi_f0s0' + `ymin'*`pi_f1s0'
	local ubY00_s0_nab = `y_00' *`pi_f0s0' + `ymax'*`pi_f1s0'
		
	local lbY00_s0_mtr = `y_00' * `pi_f0s0' + `ymin'*`pi_f1s0'
	local ubY00_s0_mtr = `y_00' * `pi_f0s0' + `y_10'*`pi_f1s0'	
	
	local lbY00_s0_mts =  `y_00' *`pi_f0s0'  + `y_00'*`pi_f1s0'
	local ubY00_s0_mts =  `y_00' *`pi_f0s0' + `ymax'*`pi_f1s0'
	
	local lbY00_s0_mtrs =  `y_00' *`pi_f0s0'  + `y_00'*`pi_f1s0'
	local ubY00_s0_mtrs =  `y_00' *`pi_f0s0' + `y_10'*`pi_f1s0'
	
	local assumption nab mtr mts mtrs
	foreach j of local assumption {
		
		local ate_s0_l_`j' =  `lbY10_s0_`j'' - (`ubY00_s0_`j'')
		local ate_s0_u_`j' =  `ubY10_s0_`j'' - (`lbY00_s0_`j'')
		
		
		display "---------------------------------------------"
		display ""
		display "ATE: `j' E(Y(1,0)-Y(0,0)|Z_s=0): " _column(16) "[" %9.4f `ate_s0_l_`j'' "," %9.4f `ate_s0_u_`j''  "]"	
		display ""
			
		matrix mate_s0_l_`j' = `ate_s0_l_`j''
		matrix mate_s0_u_`j' = `ate_s0_u_`j''
		
		return matrix mate_s0_`j'_l = mate_s0_l_`j'
		return matrix mate_s0_`j'_u = mate_s0_u_`j'
	}
	
	
	
	 
	
end


program m_sqrt, rclass
// bound are like mate_s0_`j'_l
    // Parse the input arguments
    syntax varname(numeric max=3) [if] [in]  [,BOUND(string)  Replications(real 9) ]
	tempvar  Y Df Dm 
	
	tokenize `varlist'
			
	gen `Y' = `1'
	gen `Df' = `2'
	gen `Dm' = `3'
	
	qui sum `Y' 
	local N = r(N)
	
    // Number of bootstrap samples
    local B = `replications'
	 
	qui cnpbound `Y' `Df' `Dm'
	
	display "we want `bound'_u"
     
	matrix __mRES_u =  r(`bound'_u) 
	matrix __mRES_l =  r(`bound'_l)
	
    forvalues b = 1/`B' {
        preserve
        bsample
        qui cnpbound `Y' `Df' `Dm'
		
        matrix temp =  r(`bound'_u)  
        matrix __mRES_u = __mRES_u , temp
		matrix temp =  r(`bound'_l)  
        matrix __mRES_l = __mRES_l , temp
        restore
    }

    matrix __mRES_u =sqrt(`N')*__mRES_u[.,2...]  // <------------ Multiply root-N
	matrix __mRES_l =sqrt(`N')*__mRES_l[.,2...]  // <------------ Multiply root-N
	

    // Calculate the mean of the results
    matrix __mu_u = (__mRES_u * J(`B', 1, 1) / `B')
	
    // Calculate the outer product of mRES
    matrix __Outter_u = __mRES_u * __mRES_u'

    // Calculate the product of mu and its transpose, scaled by `B'
    matrix __outMu_u = `B' * __mu_u * __mu_u'

    // Calculate the covariance matrix
    matrix __sigma_u = (1 / (`B' - 1)) * (__Outter_u - __outMu_u)
	 
	matrix list __mRES_l
	matrix list __mRES_u
     
    // Eigenvalue decomposition to find square root of matrix
    matrix symeigen __X_u __v_u = __sigma_u
    mata: eigenvalue_adjust("__v_u")
	

    // Calculate the square root of the covariance matrix
    matrix __sigmaHalf_u = __X_u * diag(__vsqrt) * __X_u'
    
    // Return the matrices
    return matrix sigma_u = __sigma_u
    return matrix sigmaHalf_u = __sigmaHalf_u
	
	//********
	// Calculate the mean of the results
    matrix __mu_l = (__mRES_l * J(`B', 1, 1) / `B')
	
    // Calculate the outer product of mRES
    matrix __Outter_l = __mRES_l * __mRES_l'

    // Calculate the product of mu and its transpose, scaled by `B'
    matrix __outMu_l = `B' * __mu_l * __mu_l'

    // Calculate the covariance matrix
    matrix __sigma_l = (1 / (`B' - 1)) * (__Outter_l - __outMu_l)
    
    // Eigenvalue decomposition to find square root of matrix
    matrix symeigen __X_l __v_l = __sigma_l
    mata: eigenvalue_adjust("__v_l")
	

    // Calculate the square root of the covariance matrix
    matrix __sigmaHalf_l = __X_l * diag(__vsqrt) * __X_l'
	
	
    // Return the matrices
    return matrix sigma_l = __sigma_l
    return matrix sigmaHalf_l = __sigmaHalf_l

end


program engine_cbound, rclass

    // Parse the input arguments
    syntax varname(numeric max=3) [if] [in]  [,BOUND(string)  Replications(real 9) ]
	tempvar  Y Df Dm 
	
	tokenize `varlist'
			
	gen `Y' = `1'
	gen `Df' = `2'
	gen `Dm' = `3'
	
	local what "`bound'"
	 

	

	m_sqrt `Y' `Df' `Dm' , bound(`what') replications(`replications')  

	matrix sigma_u = r(sigma_u)
	matrix sigma_l = r(sigma_l)
	matrix sigmaHalf_u = r(sigmaHalf_u)
	matrix sigmaHalf_l = r(sigmaHalf_l)

	cnpbound  highIQ_kid highIQ_dad highIQ_mum 
	matrix bound_u= r(`what'_u) 
	matrix bound_l= r(`what'_l) 

	mata: rho_adjust(9000, "bound_l", "sigma_l", "sigmaHalf_l",1, 0.5)
	scalar __lb = __bnd
	scalar __lb_pc = __ci
	mata: rho_adjust(9000, "bound_u", "sigma_u", "sigmaHalf_u",0, 0.5)
	scalar __ub = __bnd
	scalar __ub_pc = __ci
	
	
	mata: rho_adjust(9000, "bound_l", "sigma_l", "sigmaHalf_l",1, 0.25)
	scalar __lb25 = __bnd
	scalar __lb25_pc = __ci
	
	mata: rho_adjust(9000, "bound_u", "sigma_u", "sigmaHalf_u",0, 0.25)
	scalar __ub25 = __bnd
	scalar __ub25_pc = __ci
	
	mata: rho_adjust(9000, "bound_l", "sigma_l", "sigmaHalf_l",1, 0.75)
	scalar __lb75 = __bnd
	scalar __lb75_pc = __ci
	
	mata: rho_adjust(9000, "bound_u", "sigma_u", "sigmaHalf_u",0, 0.75)
	scalar __ub75 = __bnd
	scalar __ub75_pc = __ci
	
	local lbpc = __lb_pc
	local ubpc = __ub_pc
	local lbpc75 = __lb75_pc
	local ubpc75 = __ub75_pc
	local lbpc25 = __lb25_pc
	local ubpc25 = __ub25_pc
	
	
	mata:phat(9000, `ubpc',  `lbpc', `ubpc75',  `lbpc75',`ubpc25',  `lbpc25' )
	local pHat = phat
	
	mata: rho_adjust(9000, "bound_l", "sigma_l", "sigmaHalf_l",1, `pHat')
	scalar __lbci = __bnd
	scalar __lbci_pc = __ci
	
	mata: rho_adjust(9000, "bound_u", "sigma_u", "sigmaHalf_u",0, `pHat')
	scalar __ubci = __bnd
	scalar __ubci_pc = __ci
	
	display __ub ", " __ub_pc ", " __ub75_pc
	display __lb_pc ", " __lb75_pc
	display __ubci_pc 
	display __lbci_pc
	
end




matrix mRES = J(8,4,.)
local subindex f s
local assumptions  nab mtr mts mtrs

local col = 1

foreach j of local subindex{
	local row = 1
	foreach i of local assumptions{
		
			display "mate_`j'0_`i'"
			
			engine_cbound highIQ_kid highIQ_dad highIQ_mum , ///
				bound("mate_`j'0_`i'") replications(199)  

			matrix mRES[`row', `col'] = __lb_pc
			matrix mRES[`row', `col'+1] = __ub_pc
			matrix mRES[`row'+1, `col'] = __lbci_pc
			matrix mRES[`row'+1, `col'+1] = __ubci_pc
			local row = `row'+ 2
	}
	local col = `col' + 2

}

matrix list mRES


local outcomedir "/Users/user/Dropbox/Econometrics/youthIntelligence/outcomes/"

matrix v = mRES
matrix rownames v= NAB . MTR . MTS . MTRS . 
matrix colnames v= lb ub lb ub
outtable using "`outcomedir'conditional_bounds_fullsample" , mat(v) nobox ///
	caption("Partial identification regions, parental IQ on child IQ; full sample") replace
	

	
	
	
cd "/Users/user/Dropbox/Econometrics/youthIntelligence/code/smaller_datasets"
local outcomedir "/Users/user/Dropbox/Econometrics/youthIntelligence/outcomes/"

local datasets data_disengage_0.dta data_disengage_1.dta data_wealthy_0.dta data_wealthy_1.dta
			   
local outputfiles cbnd_parent_diseng0 cbnd_parent_diseng1 cbnd_wealthy_0 cbnd_wealth_1
local caption1 "Partial identification regions, parental IQ on child IQ; disengaged parents"
local caption2 "Partial identification regions, parental IQ on child IQ; engaged parents"
local caption3 "Partial identification regions, parental IQ on child IQ; poor parents"
local caption4 "Partial identification regions, parental IQ on child IQ; wealthy parents"

			   
local B= 199		   
local w = 1

foreach k of local datasets{

    * Load appropriate dataset
    local outputfile : word `w' of `outputfiles'
    
	use "`k'", clear
	
	// Normal loop
	
		
	local subindex f s
	local assumptions nab mtr mts mtrs

	matrix mRES = J(8,4,.)

	local col = 1

	foreach j of local subindex{
	local row = 1
	foreach i of local assumptions{
		
			display "mate_`j'0_`i'"
			
			engine_cbound highIQ_kid highIQ_dad highIQ_mum , ///
				bound("mate_`j'0_`i'") replications(199)  

			matrix mRES[`row', `col'] = __lb_pc
			matrix mRES[`row', `col'+1] = __ub_pc
			matrix mRES[`row'+1, `col'] = __lbci_pc
			matrix mRES[`row'+1, `col'+1] = __ubci_pc
			local row = `row'+ 2
	}
	local col = `col' + 2

	}
	matrix list mRES 

	local outcomedir "/Users/user/Dropbox/Econometrics/youthIntelligence/outcomes/"

	matrix v = mRES
	matrix rownames v= NAB . MTR . MTS . MTRS . 
	matrix colnames v= lb ub lb ub
	outtable using "`outcomedir'`outputfile'" , mat(v) nobox ///
		caption(`caption`w'') replace

	
	// end normal loop
	
	local w = `w'+1
}
