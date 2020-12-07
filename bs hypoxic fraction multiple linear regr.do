tempname sim
postfile `sim' p_dose coef_dose p_hf coef_hf p_ic50anox coef_ic50anox p_vdt coef_vdt using resout,replace
use datacp506-nomissingER.dta, clear
gen hf_simul=.

set output error
set seed 17082020

local hf1=17.79
local sd1=6.58
local hf2=16.89
local sd2=4.75
local hf3=5.91  // MDA-MB-231 exp 1
local sd3=0.92 // MDA-MB-231 exp 1
local hf4=18.33
local sd4=3.4   //avarage of all types
local hf5=11.83
local sd5=1.5
local hf8=10.2	
local sd8=2.08
local hf9=5.61	
local sd9=3.27
local hf10=3.79	
local sd10=2.51
local hf11=17.55
local sd11=2.03
local hf12=7.64	
local sd12=3.97
local hf13=25.65
local sd13=6.28
local hf14=22.61  // H1650 exp 1	
local sd14=3.68  // H1650 exp 1
local hf15=7.97	
local sd15=3.66
local hf16=8.43	
local sd16=4.18
local hf17=0
local sd17=0
local hf18=0
local sd18=0
local hf19=23.87	// H1650 exp 2
local sd19=5.11  // H1650 exp2
local hf20=6.11  // MDA-MB-231 exp 2
local sd20=1.1  // MDA-MB-231 exp 2

gen p_dose=.
gen coef_dose=.
gen p_hf=.
gen coef_hf=.
gen p_ic50anox=.
gen coef_ic50anox=.
gen p_vdt=.
gen coef_vdt=.

forvalues global=1/500 {
	replace hf_sim=`hf1' + `sd1'* invnorm(uniform()) if tu_id==1
	replace hf_sim=`hf2' + `sd2'* invnorm(uniform()) if tu_id==2
	replace hf_sim=`hf3' + `sd3'* invnorm(uniform()) if tu_id==3
	replace hf_sim=`hf4' + `sd4'* invnorm(uniform()) if tu_id==4
	replace hf_sim=`hf5' + `sd5'* invnorm(uniform()) if tu_id==5
	replace hf_sim=`hf8' + `sd8'* invnorm(uniform()) if tu_id==8
	replace hf_sim=`hf9' + `sd9'* invnorm(uniform()) if tu_id==9
	replace hf_sim=`hf10' + `sd10'* invnorm(uniform()) if tu_id==10
	replace hf_sim=`hf11' + `sd11'* invnorm(uniform()) if tu_id==11
	replace hf_sim=`hf12' + `sd12'* invnorm(uniform()) if tu_id==12
	replace hf_sim=`hf13' + `sd13'* invnorm(uniform()) if tu_id==13
	replace hf_sim=`hf14' + `sd14'* invnorm(uniform()) if tu_id==14
	replace hf_sim=`hf15' + `sd15'* invnorm(uniform()) if tu_id==15
	replace hf_sim=`hf16' + `sd16'* invnorm(uniform()) if tu_id==16
	replace hf_sim=`hf17' + `sd17'* invnorm(uniform()) if tu_id==17
	replace hf_sim=`hf18' + `sd18'* invnorm(uniform()) if tu_id==18
	replace hf_sim=`hf19' + `sd19'* invnorm(uniform()) if tu_id==19
	replace hf_sim=`hf20' + `sd20'* invnorm(uniform()) if tu_id==20
	
	replace hf_sim=0 if hf_simul<0
	regress SGD2 abs_cum_dose hf_sim ic50anox vdt
	
	local t = _b[abs_cum_dose]/_se[abs_cum_dose]
	replace p_dose =2*ttail(e(df_r),abs(`t'))
	replace coef_dose=_b[abs_cum_dose]
	
	local t = _b[hf_sim]/_se[hf_sim]
	replace p_hf =2*ttail(e(df_r),abs(`t'))
	replace coef_hf=_b[hf]
	
	local t = _b[ic50anox]/_se[ic50anox]
	replace p_ic50anox =2*ttail(e(df_r),abs(`t'))
	replace coef_ic50anox=_b[ic50anox]
	
	local t = _b[vdt]/_se[vdt]
	replace p_vdt =2*ttail(e(df_r),abs(`t'))
	replace coef_vdt=_b[vdt]
	
	post `sim' (p_dose) (coef_dose) (p_hf) (coef_hf) (p_ic50anox) (coef_ic50anox) (p_vdt) (coef_vdt)
	
}

postclose `sim'

