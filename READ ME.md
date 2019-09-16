READ ME

***CLTP predictions:
./CLPT.exe + power spectrum file

**Output 

***Xi_r_CLPT.dat
| Column | Value                                 
|--------+---------------------------------------|
|      1 | $r$ (in $\mathrm{Mpc\ h}^{-1}$)       
|      2 | \xi_m_f0
|      3 | \xi_m_f1 	       
|      4 | \xi_m_f2
|      5 | \xi_m_f1f1 	      
|      6 | \xi_m_f1f2 		
|      7 | \xi_m_f2f2 		

***V_12_CLPT.dat
| Column | Value                                 
|--------+---------------------------------------|
|      1 | $r$ (in $\mathrm{Mpc\ h}^{-1}$)       
|      2 | \v12_f0
|      3 | \v12_f1 	       
|      4 | \v12_f2
|      5 | \v12_f1f1 	      
|      6 | \v12_f1f2 			 

***Sigma_12_CLPT.dat
| Column | Value                                 
|--------+---------------------------------------|
|      1 | $r$ (in $\mathrm{Mpc\ h}^{-1}$)       
|      1 | $r$ (in $\mathrm{Mpc\ h}^{-1}$)       
|      2 | \s12_par_f0
|      3 | \s12_par_f1 	       
|      4 | \s12_par_f2
|      5 | \s12_par_f1f1 	      
|      6 | \s12_perp_f0
|      7 | \s12_perp_f1 	       
|      8 | \s12_perp_f2
|      9 | \s12_perp_f1f1 	      	 


** Input

***./GSRSD.exe + Parameter file 

**Parameters files 
smin
smax
b 									
f 				
f1
f2
sig_shift					
powerspectrum file


**Output 

***Xi_m.dat
| Column | Value                                 
|--------+---------------------------------------|
|      1 | $r$ (in $\mathrm{Mpc\ h}^{-1}$)       
|      2 | $\xi_m										) 		 

***V12.dat
| Column | Value                                 
|--------+---------------------------------------|
|      1 | $r$ (in $\mathrm{Mpc\ h}^{-1}$)       
|      2 | $\V12(without factor (f*b/pi*pi) 		 

***Psipar.dat
| Column | Value                                 
|--------+---------------------------------------|
|      1 | $r$ (in $\mathrm{Mpc\ h}^{-1}$)       
|      2 | $\Psipar(without factor (f*f/2*pi*pi) 		 

***Psiper.dat
| Column | Value                                 
|--------+---------------------------------------|
|      1 | $r$ (in $\mathrm{Mpc\ h}^{-1}$)       
|      2 | $\Psiper(without factor (f*f/2*pi*pi) 		 


***multipole.dat and multipole_CLPT.dat :
| Column | Value                                 
|--------+---------------------------------------|
|      1 | $r$ (in $\mathrm{Mpc\ h}^{-1}$)       
|      2 | $\xi_0											 
|      3 | $\xi_2 		                               				
|      4 | $\xi_4			                     	       
	
** Xi_s.dat and Xi_s_CLPT.dat
| Column | Value                                 |
|--------+---------------------------------------|
|      1 | $r$ (in $\mathrm{Mpc\ h}^{-1}$)       |
|      2 | $\s_par (in $\mathrm{Mpc\ h}^{-1}$) 	 
|      3 | $\s_per (in $\mathrm{Mpc\ h}^{-1}$)	     				
|      4 | $\xi_s			                     	       	
	
