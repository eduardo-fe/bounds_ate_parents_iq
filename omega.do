/*
This code, together with boundsPeerEffects_v3 is the one I used
to produce the results in the Ecnomic Letters paper (Revisions -PUblished)
*/

set matsize 10000



mata:
mata drop eigenvalue_adjust()
mata drop phat()
mata drop rho_adjust()
end

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

program drop npbpe
program drop m_sqrt
program drop engine
program drop npbound

program npbpe, rclass

	syntax varname(numeric max=3) [if] [in] [, REGressors(varlist) AT(numlist) BW(numlist)]
	tempvar  Y Df Dm 
	
	tokenize `varlist'
			
	gen `Y' = `1'
	gen `Df' = `2'
	gen `Dm' = `3'
	
	
	marksample touse 
	qui {
	if ("`regressors'" == ""){
		sum `Y' if `touse'
		local Ey=r(mean) 

		sum `Df' if `touse' 
		local pf1 = r(mean)
		local pf0 = 1-`pf1'
			
		sum `Dm' if `touse' 
		local pm1 = r(mean)
		local pm0 = 1-`pm1'
		
		gen __aux = (`Df' == 1 & `Dm' == 1)
		sum __aux if `touse'  
		local p11 = r(mean)
		drop __aux  
		
		gen __aux = (`Df' == 1 & `Dm' == 0)
		sum __aux if `touse'
		local p10 = r(mean)
		drop __aux
		
		gen __aux = (`Df' == 0 & `Dm' == 1)
		sum __aux if `touse'
		local p01 = r(mean)
		drop __aux
		
		gen __aux = (`Df' == 0 & `Dm' == 0)
		sum __aux if `touse'
		local p00 = r(mean)
		drop __aux
		
		
		local P11 = `p11'*100
		local P10 = `p10'*100
		local P01 = `p01'*100
		local P00 = `p00'*100
		local Q11 = 100-`P11'
		local Q10 = 100-`P10'
		local Q01 = 100-`P01'
		local Q00 = 100-`P00'
		
		/* This is different from the same file in folder "code" */
		sum `Y' if `touse'
		local a = r(min)
	
		sum `Y' if `touse'
		local b = r(max)
					
		sum `Y' if `Df' == 1 & `Dm' == 1 & `touse'
		local Ey11 = r(mean)
		
		sum `Y' if `Df' == 0 & `Dm' == 1 & `touse'
		local Ey01 = r(mean)
		
		sum `Y' if `Df' == 1 & `Dm' == 0 & `touse'
		local Ey10 = r(mean)

		sum `Y' if `Df' == 0 & `Dm' == 0 & `touse'
		local Ey00 = r(mean)
		
		sum `Y' if `Df' == 1 & `touse'
		local Eyf1 = r(mean)
	
		sum `Y' if `Df' == 0 & `touse'
		local Eyf0 = r(mean)
		
		sum `Y' if `Dm' == 1 & `touse'
		local Eym1 = r(mean)
	
		sum `Y' if `Dm' == 0 & `touse'
		local Eym0 = r(mean)
	}
		 
	else {
		liRacine  `Y' `regressors' if `touse', at(`at')   
		local Ey=r(ytilde)
		 display "`Ey'"

		liRacine `Df' `regressors' if `touse', at(`at') 
		local pf1 = r(ytilde)
		local pf0 = 1-`pf1'
		
		liRacine `Dm' `regressors' if `touse', at(`at')
		local pm1 = r(ytilde)
		local pm0 = 1-`pm1'
		
		 
		
		gen __aux = (`Df' == 1 & `Dm' == 1)
		liRacine __aux `regressors' if `touse', at(`at')
		local p11 = r(ytilde)
		drop __aux
			
		gen __aux = (`Df' == 1 & `Dm' == 0)
		liRacine __aux `regressors' if `touse', at(`at')
		local p10 = r(ytilde)
		drop __aux
		
		gen __aux = (`Df' == 0 & `Dm' == 1)
		liRacine __aux `regressors' if `touse', at(`at')
		local p01 = r(ytilde)
		drop __aux
		
		gen __aux = (`Df' == 0 & `Dm' == 0)
		liRacine __aux `regressors' if `touse', at(`at')
		local p00 = r(ytilde)
		drop __aux
		
		 
		
		local P11 = `p11'*100
		local P10 = `p10'*100
		local P01 = `p01'*100
		local P00 = `p00'*100
		local Q11 = 100-`P11'
		local Q10 = 100-`P10'
		local Q01 = 100-`P01'
		local Q00 = 100-`P00'
		
		/* This is different from the code in the folder "code"*/
		sum `Y' if `touse'
		local a = r(min)
	
		sum `Y' if `touse'
		local b = r(max)
		
		
		
		liRacine `Y' `regressors' if `Df' == 1 & `Dm' == 1 & `touse', at(`at')
		local Ey11 = r(ytilde)
		
		liRacine `Y' `regressors' if `Df' == 0 & `Dm' == 1 & `touse', at(`at')
		local Ey01 = r(ytilde)
		
		liRacine `Y' `regressors' if `Df' == 1 & `Dm' == 0 & `touse', at(`at')
		local Ey10 = r(ytilde)

		liRacine `Y' `regressors' if `Df' == 0 & `Dm' == 0 & `touse', at(`at')
		local Ey00 = r(ytilde)
		
		liRacine `Y' `regressors' if `Df' == 1 & `touse', at(`at')
		local Eyf1 = r(ytilde)
	
		liRacine `Y' `regressors' if `Df' == 0 & `touse', at(`at')
		local Eyf0 = r(ytilde)
		
		liRacine `Y' `regressors' if `Dm' == 1 & `touse', at( `at')
		local Eym1 = r(ytilde)
	
		liRacine `Y' `regressors' if `Dm' == 0 & `touse', at(`at')
		local Eym0 = r(ytilde)
	
	 }
		 
	}	 
	
	
			
			
	//local methods zr cmtr cmtrg  ots mts cmtrMts  cmtrOts cmtrgOts otsMts cMtrOtsMts
	local methods zr cmtr   mts cmtrMts  
	foreach j of local methods {
		if("`j'" == "zr") {
			 display "ZR"
			
			local `j'Y00_l = `Ey00'*`p00'+`a'*`p10' +`a'*`p01'+`a'*`p11'
			local `j'Y00_u = `Ey00'*`p00'+`b'*`p10' +`b'*`p01'+`b'*`p11'
		
			local `j'Y11_l = `Ey11'*`p11'+ `a'*`p10' +`a'*`p01' +`a'*`p00'
			local `j'Y11_u = `Ey11'*`p11'+ `b'*`p10' +`b'*`p01' +`b'*`p00'
			
			local `j'Y10_l = `Ey10'*`p10'+ `a'*`p11'+`a'*`p01' + `a'*`p00'
			local `j'Y10_u = `Ey10'*`p10'+ `b'*`p11'+`b'*`p01' + `b'*`p00'
			
			local `j'Y01_l = `Ey01'*`p01'+ `a'*`p11'+`a'*`p10' + `a'*`p00'
			local `j'Y01_u = `Ey01'*`p01'+ `b'*`p11'+`b'*`p10' + `b'*`p00'
			
			matrix m`j'atef_l = ``j'Y10_l' - ``j'Y00_u'
			matrix m`j'atef_u = ``j'Y10_u' - ``j'Y00_l'
			
			matrix m`j'atem_l = ``j'Y01_l' - ``j'Y00_u'
			matrix m`j'atem_u = ``j'Y01_u' - ``j'Y00_l'
			
			matrix m`j'atej_l = ``j'Y11_l' - ``j'Y00_u'
			matrix m`j'atej_u = ``j'Y11_u' - ``j'Y00_l'
			
			matrix m`j'ateif_l = ``j'Y11_l' - ``j'Y01_u'
			matrix m`j'ateif_u = ``j'Y11_u' - ``j'Y01_l'
			
			matrix m`j'ateim_l = ``j'Y11_l' - ``j'Y10_u'
			matrix m`j'ateim_u = ``j'Y11_u' - ``j'Y10_l'
			
		 
			
		}
		else if("`j'" == "cmtr") {
			 
			display "CMTR"
			local `j'Y00_u = `Ey' 
			local `j'Y00_l = `Ey00'*`p00'+`a'*(1-`p00')
		
			local `j'Y11_u = `Ey11'*`p11'+ `b'*(1-`p11')
			local `j'Y11_l = `Ey'
			
			local `j'Y10_u = `Ey10'*`p10' +`Ey11'*`p11'+`b'*`p01' + `b'*`p00'
			local `j'Y10_l = `Ey10'*`p10' +`Ey00'*`p00' +`a'*`p01' + `a'*`p11'
			
			local `j'Y01_u = `Ey01'*`p01' +`Ey11'*`p11' +`b'*`p10'+`b'*`p00'
			local `j'Y01_l = `Ey01'*`p01' +`Ey00'*`p00'+`a'*`p10' + `a'*`p11'
			
			matrix m`j'atef_l = ``j'Y10_l' - ``j'Y00_u'
			matrix m`j'atef_u = ``j'Y10_u' - ``j'Y00_l'
			
			matrix m`j'atem_l = ``j'Y01_l' - ``j'Y00_u'
			matrix m`j'atem_u = ``j'Y01_u' - ``j'Y00_l'
			
			matrix m`j'atej_l = ``j'Y11_l' - ``j'Y00_u'
			matrix m`j'atej_u = ``j'Y11_u' - ``j'Y00_l'
			
			matrix m`j'ateif_l = ``j'Y11_l' - ``j'Y01_u'
			matrix m`j'ateif_u = ``j'Y11_u' - ``j'Y01_l'
			
			matrix m`j'ateim_l = ``j'Y11_l' - ``j'Y10_u'
			matrix m`j'ateim_u = ``j'Y11_u' - ``j'Y10_l'
			
		
			
			
			
			
		}
		else if ("`j'" == "mts"){
			 
			
			display "MTS."
			local `j'Y00_u = `Ey00'*`p00' + `b'*`p10'+ `b'*`p01'+ `b'*`p11' 
			local `j'Y00_l = `Ey00'
							
			local `j'Y11_u = `Ey11' 
			local `j'Y11_l = `Ey11'*`p11' + `a'*`p10'+ `a'*`p01'+ `a'*`p00'
			
			local `j'Y10_u = `Ey10'*(`p10'+`p00')+`b'*`p01' + `b'*`p11'
			local `j'Y10_l = `Ey10'*(`p10'+`p11')+`a'*`p01' + `a'*`p00'
				
			local `j'Y01_u = `Ey01'*(`p01'+`p00')+`b'*`p10' + `b'*`p11'
			local `j'Y01_l = `Ey01'*(`p01'+`p11')+`a'*`p10' + `a'*`p00'
			
			matrix m`j'atef_l = ``j'Y10_l' - ``j'Y00_u'
			matrix m`j'atef_u = ``j'Y10_u' - ``j'Y00_l'
			
			matrix m`j'atem_l = ``j'Y01_l' - ``j'Y00_u'
			matrix m`j'atem_u = ``j'Y01_u' - ``j'Y00_l'
			
			matrix m`j'atej_l = ``j'Y11_l' - ``j'Y00_u'
			matrix m`j'atej_u = ``j'Y11_u' - ``j'Y00_l'
			
			matrix m`j'ateif_l = ``j'Y11_l' - ``j'Y01_u'
			matrix m`j'ateif_u = ``j'Y11_u' - ``j'Y01_l'
			
			matrix m`j'ateim_l = ``j'Y11_l' - ``j'Y10_u'
			matrix m`j'ateim_u = ``j'Y11_u' - ``j'Y10_l'
			
			  
		}
		else if ("`j'" == "cmtrMts"){
			display "CMTR +MTR "
			
			local `j'Y00_u = `Ey' 
			local `j'Y00_l = `Ey00'
			local `j'Y11_u = `Ey11'
			local `j'Y11_l = `Ey'
			local `j'Y10_u = `Ey10'*(`p10'+`p00') + `b'*`p01' + `Ey11'*`p11'
			local `j'Y10_l = `Ey10'*(`p10'+`p11') + `a'*`p01' + `Ey00'*`p00'
			local `j'Y01_u = `Ey01'*(`p01'+`p00') + `b'*`p10'+`Ey11'*`p11' 
			local `j'Y01_l = `Ey01'*(`p01'+`p11') + `a'*`p10' + `Ey00'*`p00'
			 
			
			matrix m`j'atef_l = ``j'Y10_l' - ``j'Y00_u'
			matrix m`j'atef_u = ``j'Y10_u' - ``j'Y00_l'
			
			matrix m`j'atem_l = ``j'Y01_l' - ``j'Y00_u'
			matrix m`j'atem_u = ``j'Y01_u' - ``j'Y00_l'
			
			matrix m`j'atej_l = ``j'Y11_l' - ``j'Y00_u'
			matrix m`j'atej_u = ``j'Y11_u' - ``j'Y00_l'
			
			matrix m`j'ateif_l = ``j'Y11_l' - ``j'Y01_u'
			matrix m`j'ateif_u = ``j'Y11_u' - ``j'Y01_l'
			
			matrix m`j'ateim_l = ``j'Y11_l' - ``j'Y10_u'
			matrix m`j'ateim_u = ``j'Y11_u' - ``j'Y10_l'
			  
			 
		}
		
			return matrix `j'atef_l=m`j'atef_l
			return matrix `j'atef_u=m`j'atef_u
			return matrix `j'atem_l=m`j'atem_l
			return matrix `j'atem_u=m`j'atem_u
			return matrix `j'atej_l=m`j'atej_l
			return matrix `j'atej_u=m`j'atej_u
			return matrix `j'ateif_l=m`j'ateif_l
			return matrix `j'ateif_u=m`j'ateif_u
			return matrix `j'ateim_l=m`j'ateim_l
			return matrix `j'ateim_u=m`j'ateim_u
			
		local `j'atef_l = ``j'Y10_l' - ``j'Y00_u'
		local `j'atef_u = ``j'Y10_u' - ``j'Y00_l'
			
		local `j'atem_l = ``j'Y01_l' - ``j'Y00_u'
		local `j'atem_u = ``j'Y01_u' - ``j'Y00_l'
			
		local `j'atej_l = ``j'Y11_l' - ``j'Y00_u'
		local `j'atej_u = ``j'Y11_u' - ``j'Y00_l'
			
		local `j'ateif_l = ``j'Y11_l' - ``j'Y01_u'
		local `j'ateif_u = ``j'Y11_u' - ``j'Y01_l'
			
		local `j'ateim_l = ``j'Y11_l' - ``j'Y10_u'
		local `j'ateim_u = ``j'Y11_u' - ``j'Y10_l'
			 
		display "---------------------------------------------"
		display ""
		display "Y_00:  " _column(16) "[" %9.4f ``j'Y00_l' "," %9.4f ``j'Y00_u' "]"
		display ""
		display "Y_10:  " _column(16) "[" %9.4f ``j'Y10_l' "," %9.4f ``j'Y10_u' "]"
		display ""
		display "Y_01:  " _column(16) "[" %9.4f ``j'Y01_l' "," %9.4f ``j'Y01_u' "]"
		display ""
		display "Y_11:  " _column(16) "[" %9.4f ``j'Y11_l' "," %9.4f ``j'Y11_u' "]"
		display ""
		display "ATE_f:  " _column(16) "[" %9.4f ``j'atef_l' "," %9.4f ``j'atef_u' "]"	
		display ""
		display "ATE_m:  " _column(16) "[" %9.4f ``j'atem_l' "," %9.4f ``j'atem_u' "]"	
		display ""
		display "ATE_j:  " _column(16) "[" %9.4f ``j'atej_l' "," %9.4f ``j'atej_u' "]"	
		display ""
		display "ATE_if:  " _column(16) "[" %9.4f ``j'ateif_l' "," %9.4f ``j'ateif_u' "]"	
		display ""
		display "ATE_im:  " _column(16) "[" %9.4f ``j'ateim_l' "," %9.4f ``j'ateim_u' "]"	
		display "---------------------------------------------"
		display ""
			
		

	}
	
	 
end



program m_sqrt, rclass

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
	 
	qui npbpe `Y' `Df' `Dm'
	
	display "we want `bound'_u"
     
	matrix __mRES_u =  r(`bound'_u) 
	matrix __mRES_l =  r(`bound'_l)
	
    forvalues b = 1/`B' {
        preserve
        bsample
        qui npbpe `Y' `Df' `Dm'
		
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


program engine, rclass

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

	qui npbpe `Y'  `Df' `Dm'
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

program npbound, rclass

	syntax varname(numeric max=2) [if] [, Positive]	
	
	tempvar  Y D 
	
	tokenize `varlist'
			
	gen `Y' = `1'
	gen `D' = `2'
	
		
	marksample touse

	qui{
	sum `Y' // Note here we compute the max in the full sample, not touse
	local ymin = r(min)
	local ymax = r(max)

	sum `Y' if `D' & `touse'
	local y_1 = r(mean)

	sum `Y' if !`D' & `touse'
	local y_0 = r(mean)

	sum `D' if  `touse'
	local p_1 = r(mean)
	local p_0 = 1-`p_1'
	
	}

	local methods nab mtr mts mtrs
	foreach j of local methods {

		if("`j'" == "nab") {
			display "No Assumption Bound."
			local `j'Y0_l = `y_0' * `p_0' + `ymin'*`p_1'
			local `j'Y0_u = `y_0' * `p_0' + `ymax'*`p_1'
			local `j'Y1_l = `y_1' * `p_1' + `ymin'*`p_0'
			local `j'Y1_u = `y_1' * `p_1' + `ymax'*`p_0'
			local `j'ate_l = ``j'Y1_l' - ``j'Y0_u'
			local `j'ate_u = ``j'Y1_u' - ``j'Y0_l'
		
			
			 
		}
		else if ("`j'" == "mtr"){
			
			display "Monotone Treatment Response"
			if ("`positive'" != "") {
				local `j'Y0_l = `y_0' * `p_0' + `ymin'*`p_1'
				local `j'Y0_u = `y_0' * `p_0' + `y_1'*`p_1'
				local `j'Y1_l = `y_1' * `p_1' + `y_0'*`p_0'
				local `j'Y1_u = `y_1' * `p_1' + `ymax'*`p_0'
			
			}
			else {
			
				local `j'Y0_l = `y_0' * `p_0' + `y_1'*`p_1'
				local `j'Y0_u = `y_0' * `p_0' + `ymax'*`p_1'
				local `j'Y1_l = `y_1' * `p_1' + `ymin'*`p_0'
				local `j'Y1_u = `y_1' * `p_1' + `y_0'*`p_0'
			}
			local `j'ate_l = ``j'Y1_l' - ``j'Y0_u'
			local `j'ate_u = ``j'Y1_u' - ``j'Y0_l'
			
			 
		}
		else if ("`j'" == "mts"){
			display "Monotone Treatment Selection"
			
			if ("`positive'" != "") {
				
				
			local `j'Y0_l = `y_0' 
				local `j'Y0_u = `y_0' * `p_0' + `ymax'*`p_1'
				local `j'Y1_l = `y_1' * `p_1' + `ymin'*`p_0'
				local `j'Y1_u = `y_1'
			}
			else {
			local `j'Y0_l = `y_0' * `p_0' + `ymin'*`p_1'
				local `j'Y0_u = `y_0'
				local `j'Y1_l = `y_1' 
				local `j'Y1_u = `y_1' * `p_1' + `ymax'*`p_0'
				
			}
			local `j'ate_l = ``j'Y1_l' - ``j'Y0_u'
			local `j'ate_u = ``j'Y1_u' - ``j'Y0_l'
			
			 
		}
		else if ("`j'" == "mtrs"){
			display "Monotone Treatment Selection and Response"
			if ("`positive'" != "") {
				local `j'Y0_l = `y_0' 
				local `j'Y0_u = `y_0' * `p_0' + `y_1'*`p_1'
				local `j'Y1_l = `y_1' * `p_1' + `y_0'*`p_0'
				local `j'Y1_u = `y_1'
			}
			else{
				local `j'Y0_l = `y_0' * `p_0' + `y_1'*`p_1'
				local `j'Y0_u = `y_0'
				local `j'Y1_l = `y_1' 
				local `j'Y1_u = `y_1' * `p_1' + `y_0'*`p_0'
			}
			local `j'ate_l = ``j'Y1_l' - ``j'Y0_u'
			local `j'ate_u = ``j'Y1_u' - ``j'Y0_l'
			 
			 
		}
		
		
		display "---------------------------------------------"
		display ""
		display "Y_0:  " _column(16) "[" %9.4f ``j'Y0_l' "," %9.4f ``j'Y0_u' "]"
		display ""
		display "Y_1:  " _column(16) "[" %9.4f ``j'Y1_l' "," %9.4f ``j'Y1_u' "]"
		display ""
		display "ATE:  " _column(16) "[" %9.4f ``j'ate_l' "," %9.4f ``j'ate_u' "]"	
		display "---------------------------------------------"
		display ""
			
		

	}

local moments Y0 Y1 ate
foreach j of local methods {

	foreach k of local moments {
	
		return scalar `j'`k'_l = ``j'`k'_l'
		return scalar `j'`k'_u = ``j'`k'_u'
	
	}

}
	
end



* Data Analysis

 
local B = 199
 
local subindex f m j if im
local assumptions1 zrate cmtrate  mtsate cmtrMtsate 

matrix mRES = J(10,8,.)


local col = 1

foreach i of local assumptions1{
local row = 1
	foreach j of local subindex{
		display "`i'`j'"
		 		
		engine highIQ_kid highIQ_dad highIQ_mum , bound("`i'`j'") replications(`B') 
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
matrix rownames v= ate_f . ate_s . ate_j . ate_if . ate_is  . 
matrix colnames v= lb_zb ub_zb lb_CMTR ub_CMTR ///
	lb_MTS ub_MTS lb_CMTRS ub_CMTRS 
outtable using "`outcomedir'bounds_all" , mat(v) nobox ///
	caption("Partial identification regions, parental IQ on child IQ; full sample") replace

	
	
	
	
* Define parameters
cd "/Users/user/Dropbox/Econometrics/youthIntelligence/code/smaller_datasets"
local outcomedir "/Users/user/Dropbox/Econometrics/youthIntelligence/outcomes/"

local statuses disengaged no_disengaged
local datasets data_disengage_0.dta data_disengage_1.dta data_wealthy_0.dta data_wealthy_1.dta
			   
local outputfiles bnd_parent_diseng0 bnd_parent_diseng1 bnd_wealthy_0 bnd_wealth_1
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
	
		
	local subindex f m j if im
	local assumptions1 zrate cmtrate  mtsate cmtrMtsate 

	matrix mRES = J(10,8,.)

	local col = 1

	foreach i of local assumptions1{
	local row = 1
		foreach j of local subindex{
			display "`i'`j'"
					
			engine ravScore_sdHigh memScore_sd_dadHigh memScore_sd_mumHigh , ///
				bound("`i'`j'") replications(`B') 
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
	matrix rownames v= ate_f . ate_s . ate_j . ate_if . ate_is  . 
	matrix colnames v= lb_zb ub_zb lb_CMTR ub_CMTR ///
		lb_MTS ub_MTS lb_CMTRS ub_CMTRS 
	outtable using "`outcomedir'`outputfile'" , mat(v) nobox ///
		caption(`caption`w'') replace

	
	
	
	// end normal loop
	
	local w = `w'+1
}

		
		
		
////////////////////////////////////////////////////////////////////////*
/*		
		
		
			   
local B = 199
local subindex f m j if im
local assumptions1 zrate cmtrate mtsate cmtrMtsate
 
* Loop over statuses
local w = 1
foreach k of local datasets{

    * Load appropriate dataset
    local outputfile : word `w' of `outputfiles'
    local caption : word `w' of `captions'
    
	use "`k'", clear
	count
    matrix mRES = J(10,8,.)

    local col = 1
    foreach i of local assumptions1 {
        local row = 1
        foreach j of local subindex {
            //display "`i'`j'"

            qui engine ravScore_sdHigh memScore_sd_dadHigh memScore_sd_mumHigh , bound("`i'`j'") replications(`B')  

            matrix mRES[`row', `col'] = __lb_pc
            matrix mRES[`row', `col'+1] = __ub_pc
            matrix mRES[`row'+1, `col'] = __lbci_pc
            matrix mRES[`row'+1, `col'+1] = __ubci_pc
            local row = `row' + 2
        }
        local col = `col' + 2
    }

    * Set up matrix labels
    matrix v = mRES
    matrix rownames v = ate_f . ate_s . ate_j . ate_if . ate_is .  
    matrix colnames v = lb_zb ub_zb lb_CMTR ub_CMTR lb_MTS ub_MTS lb_CMTRS ub_CMTRS
					
    * Output table
	display "To here"
	display "`outputfile'"
    outtable using "`outcomedir'`outputfile'", mat(v) nobox    replace
	
	local w= `w' + 1
	display `w'

}



	
	
	
* Parent disengaged

use "/Users/user/Dropbox/Econometrics/youthIntelligence/code/smaller_datasets/data_disengage_0.dta", clear


local B = 99

local subindex f m j if im
local assumptions1 zrate cmtrate  mtsate cmtrMtsate 
 
matrix mRES = J(10,8,.)


local col = 1

foreach i of local assumptions1{
local row = 1
	foreach j of local subindex{
		display "`i'`j'"
		 		
		engine ravScore_sdHigh memScore_sd_dadHigh memScore_sd_mumHigh , bound("`i'`j'") replications(`B')  

		matrix mRES[`row', `col'] = __lb_pc
		matrix mRES[`row', `col'+1] = __ub_pc
		matrix mRES[`row'+1, `col'] = __lbci_pc
		matrix mRES[`row'+1, `col'+1] = __ubci_pc
		local row = `row'+ 2
	}
	local col = `col' + 2

}


local outcomedir "/Users/user/Dropbox/Econometrics/youthIntelligence/outcomes/"

matrix v = mRES
matrix rownames v= ate_f . ate_s . ate_j . ate_if . ate_is . ate_f . ate_s . ate_j . ate_if . ate_is . 
matrix colnames v= lb_zb ub_zb lb_CMTR ub_CMTR ///
	lb_MTS ub_MTS lb_CMTRS ub_CMTRS 
outtable using "`outcomedir'bounds_parent_diseng1" , mat(v) nobox ///
	caption("Partial identification regions, parental IQ on child IQ; full sample") replace

 

 
 * Parent no disengaged
 
 
use "/Users/user/Dropbox/Econometrics/youthIntelligence/code/smaller_datasets/data_disengage_0.dta", clear


local B = 99

local subindex f m j if im
local assumptions1 zrate cmtrate  mtsate cmtrMtsate 

matrix mRES = J(10,8,.)


local col = 1

foreach i of local assumptions1{
local row = 1
	foreach j of local subindex{
		display "`i'`j'"
		 		
		engine ravScore_sdHigh memScore_sd_dadHigh memScore_sd_mumHigh , bound("`i'`j'") replications(`B')  

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
matrix rownames v= ate_f . ate_s . ate_j . ate_if . ate_is . ate_f . ate_s . ate_j . ate_if . ate_is . 
matrix colnames v= lb_zb ub_zb lb_CMTR ub_CMTR ///
	lb_MTS ub_MTS lb_CMTRS ub_CMTRS
outtable using "`outcomedir'bounds_parent_diseng0" , mat(v) nobox ///
	caption("Partial identification regions, parental IQ on child IQ; full sample") replace


* WEalthy family

 
use "/Users/user/Dropbox/Econometrics/youthIntelligence/code/smaller_datasets/data_wealthy_1.dta", clear


local B = 99

local subindex f m j if im
local assumptions1 zrate cmtrate  mtsate  cmtrMtsate 

matrix mRES = J(10,8,.)


local col = 1

foreach i of local assumptions1{
local row = 1
	foreach j of local subindex{
		display "`i'`j'"
		 		
		engine ravScore_sdHigh memScore_sd_dadHigh memScore_sd_mumHigh , bound("`i'`j'") replications(`B')  

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
matrix rownames v= ate_f . ate_s . ate_j . ate_if . ate_is . ate_f . ate_s . ate_j . ate_if . ate_is . 
matrix colnames v= lb_zb ub_zb lb_CMTR ub_CMTR ///
	lb_MTS ub_MTS lb_CMTRS ub_CMTRS
outtable using "`outcomedir'bounds_parent_wealthy1" , mat(v) nobox ///
	caption("Partial identification regions, parental IQ on child IQ; full sample") replace


*poor family


use "/Users/user/Dropbox/Econometrics/youthIntelligence/code/smaller_datasets/data_wealthy_0.dta", clear


local B = 99
local alph = 0.9

local subindex f m j if im
local assumptions1 zrate cmtrate mtsate cmtrMtsate 

matrix mRES = J(20,8,.)


local col = 1

foreach i of local assumptions1{
local row = 1
	foreach j of local subindex{
		display "`i'`j'"
		 		
		engine ravScore_sdHigh memScore_sd_dadHigh memScore_sd_mumHigh , bound("`i'`j'") replications(`B')  

		matrix mRES[`row', `col'] = __lb_pc
		matrix mRES[`row', `col'+1] = __ub_pc
		matrix mRES[`row'+1, `col'] = __lbci_pc
		matrix mRES[`row'+1, `col'+1] = __ubci_pc
		local row = `row'+ 2
	}
	local col = `col' + 2

}
matrix list mRES

local subindex f m j if im
local assumptions2  cmtrMtsate cmtrOtsate motsMtsate mcMtrOtsMtsate  //

local col = 1

foreach i of local assumptions2{
local row = 11
	foreach j of local subindex{
		display "`i'`j'"
		
		engine ravScore_sdHigh memScore_sd_dadHigh memScore_sd_mumHigh , bound("`i'`j'") replications(`B')  


		matrix mRES[`row', `col'] = __lb_pc
		matrix mRES[`row', `col'+1] = __ub_pc
		matrix mRES[`row'+1, `col'] = __lbci_pc
		matrix mRES[`row'+1, `col'+1] = __ubci_pc
		local row = `row'+ 2
	}
	local col = `col' + 2

}


local outcomedir "/Users/user/Dropbox/Econometrics/youthIntelligence/outcomes/"

matrix v = mRES
matrix rownames v= ate_f . ate_s . ate_j . ate_if . ate_is . ate_f . ate_s . ate_j . ate_if . ate_is . 
matrix colnames v= lb_zb_MTRMTS ub_zb_MTRMTS lb_MTR_MTROTS ub_MTR_MTROTS ///
	lb_oTS_MTSOTS ub_oTS_MTSOTS lb_MTS_MTROTSMTS ub_MTS_MTROTSMTS 
outtable using "`outcomedir'bounds_parent_wealthy0" , mat(v) nobox ///
	caption("Partial identification regions, parental IQ on child IQ; full sample") replace




/*

program bindate, rclass

	syntax varname(numeric max=2) [if] [, Positive]	
	
	tokenize `varlist'
			
	gen __Y = `1'
	gen __D = `2'
	gen __fulldata = __Y!=. & __D!=.
	
	levelsof __D, local(_lvls)
	mata: bindate()
	
end

mata: 
mata drop bindate()
end

mata:
real void bindate() {
    
	real matrix Y, D, fullset, lvls,  mean_y_given_z, p_z, p_z_lt, p_z_gt, p_z_t,  mean_y_t, Nates
    real scalar y_min, y_max, total_obs, k, Nates_i
    real matrix lower_bound, upper_bound
      
	fullset = st_data(., "__fulldata")  // Full dataset
    Y = st_data(., "__Y", "__fulldata")  // Outcome variable
    D = st_data(., "__D", "__fulldata")  // Treatment or grouping variable

    // Step 2: Get distinct levels of D
    lvls = uniqrows(D)

	k= rows(lvls)
	k

 // Step 3: Calculate probabilities P(D = s), P(D < t), P(D > t)
    p_z = J(k, 1, 0)  // Initialize P(D = s)
	p_z_lt = J(k, 1, 0)  // Initialize P(D < t)
    p_z_gt = J(k, 1, 0)  // Initialize P(D > t)
	
	total_obs = rows(D)
	
	
    
    for (i = 1; i <= k; i++) {
		
		
        // P(D = s)
        p_z[i] = sum(D :== lvls[i]) / total_obs

        // Cumulative probabilities
        p_z_lt[i] = sum(D :< lvls[i]) / total_obs
        p_z_gt[i] = sum(D :> lvls[i]) / total_obs
    }
	
	// Step 4: Calculate mean(Y | D = s)
    mean_y_given_z = J(rows(lvls), 1, 0)  // Initialize mean(Y | D = s)
    for (i = 1; i <= rows(lvls); i++) {
        mean_y_given_z[i] = mean(select(Y, D :== lvls[i]))
    }

    // Step 5: Define y_min and y_max
    y_min = min(Y)
    y_max = max(Y)

	
    // Step 6: Calculate Manski bounds
    lower_bound = J(rows(lvls), 1, 0)  // Initialize lower bound
    upper_bound = J(rows(lvls), 1, 0)  // Initialize upper bound
	lower_bound_mtsr = J(rows(lvls), 1, 0)  // Initialize lower bound
    upper_bound_mtsr = J(rows(lvls), 1, 0)  // Initialize upper bound
    for (t = 1; t <= rows(lvls); t++) {
        // Lower bound
		
		p_z_t = p_z[1::t]  // Probabilities P(z = s) for s <= t
		mean_y_t = mean_y_given_z[1::t]  // Means E[y | z = s) for s <= t

        lower_bound[t] = sum(mean_y_t :* p_z_t) + y_min * p_z_gt[t]
		lower_bound_mtsr[t] = sum(mean_y_t :* p_z_t) +  mean_y_given_z[t]* p_z_gt[t]

        // Upper bound
		p_z_t = p_z[t::k]  // Probabilities P(z = s) for s >= t
		mean_y_t = mean_y_given_z[t::k]  // Means E[y | z = s) for s >= t

    // Compute the sum for current t
    
        upper_bound[t] = sum(mean_y_t :* p_z_t)  + y_max * p_z_lt[t]
		upper_bound_mtsr[t] = sum(mean_y_t :* p_z_t)  + mean_y_given_z[t] * p_z_lt[t]
    }
	lower_bound
	upper_bound
	lower_bound_mtsr
	upper_bound_mtsr
    // Step 7: Output results
    st_matrix("lower_bound", lower_bound)  // Save lower bounds to Stata matrix
    st_matrix("upper_bound", upper_bound)  // Save upper bounds to Stata matrix
	
	N = rows(lower_bound)
	Nates = N*(N-1)/2
	lower_ATE = J(Nates, 3, .) // Matrix to store lower bounds for ATEs
	upper_ATE = J(Nates, 3, .) // Matrix to store upper bounds for ATEs

// Loop over all pairs (z, z') with z > z'
	Nates_i = 1
	for (z = 2; z <= N; z++) {
		for (z_prime = 1; z_prime < z; z_prime++) {
			// Compute bounds for E(Y(z) - Y(z'))
			lower_ATE[Nates_i,3] = lower_bound[z] - upper_bound[z_prime]
			upper_ATE[Nates_i,3] = upper_bound[z] - lower_bound[z_prime]
			lower_ATE[Nates_i,1] =z
			lower_ATE[Nates_i,2] =z_prime
			upper_ATE[Nates_i,1] =z
			upper_ATE[Nates_i,2] =z_prime
			
			Nates_i =Nates_i+1
		}
	}

	// Output the results
	lower_ATE
	upper_ATE
	
	
	// Loop over all pairs (z, z') with z > z'
	Nates_i = 1
	for (z = 2; z <= N; z++) {
		for (z_prime = 1; z_prime < z; z_prime++) {
			// Compute bounds for E(Y(z) - Y(z'))
			lower_ATE[Nates_i,3] = lower_bound_mtsr[z] - upper_bound_mtsr[z_prime]
			upper_ATE[Nates_i,3] = upper_bound_mtsr[z] - lower_bound_mtsr[z_prime]
			lower_ATE[Nates_i,1] =z
			lower_ATE[Nates_i,2] =z_prime
			upper_ATE[Nates_i,1] =z
			upper_ATE[Nates_i,2] =z_prime
			
			Nates_i =Nates_i+1
		}
	}

	// Output the results
	lower_ATE
	upper_ATE
	
		
    }
end


