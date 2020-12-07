# Selectively Targeting Tumor Hypoxia by CP-506, a Next-Generation Hypoxia-Activated Prodrug

Alexander M.A. van der Wiel1†, Victoria Jackson-Patel2,3†, Raymon Niemans1†, Ala Yaromina1, Emily Lui2, Damiënne Marcus1, Alexandra M. Mowday1,2, Natasja G. Lieuwes1, Rianne Biemans1, Xiaojing Lin2,3, Arthur Jochems1, Amir Ashoorzadeh2,3, Robert F. Anderson2,3, Kevin O. Hicks2,3, Matthew Bull2,3, Maria M. Abbatista2,3, Christopher P. Guise2,3, Arne Heyerick4, Morwena Solivio5, Silvia Balbo5, Jeff B. Smaill2,3, Jan Theys1, Ludwig J. Dubois1‡, Adam V. Patterson2,3‡, and Philippe Lambin1‡.
1 The D-Lab and The M-Lab, Department of Precision Medicine, GROW – School for Oncology and Developmental Biology, Maastricht University, Maastricht, the Netherlands
2 Auckland Cancer Society Research Centre, School of Medical Sciences, University of Auckland, Auckland, Private Bag 92019, Auckland 1142, New Zealand
3 Maurice Wilkins Centre for Molecular Biodiscovery, University of Auckland, Auckland, Private Bag 92019, Auckland 1142, New Zealand
4 Convert Pharmaceuticals, Liège, Belgium
5 Masonic Cancer Center, University of Minnesota, Minneapolis, USA
† Indicates equal contribution; ‡ Indicates senior co-authorship


## Abstract
Hypoxia is a pervasive feature of human tumors that is acknowledged as a major impediment to treatment success. Hypoxia-activated prodrugs (HAP) are a promising class of antineoplastic agents that can selectively eliminate hypoxic tumor cells. The present study evaluated the hypoxia-selectivity and antitumor activity of CP-506, a DNA alkylating HAP with favorable properties. Radiolytic methods revealed CP-506 to undergo oxygen-reversible one-electron reduction; multi-electron reduction products were readily formed in the absence of oxygen. Cellular metabolism was inhibited completely by oxygen concentrations above 1 µM. CP-506 was activated under anoxia by human diflavin oxidoreductases but was resistant to aerobic activation by human aldo-keto reductase 1C3. CP-506 demonstrated cytotoxicity selectively in hypoxic 2D and 3D cell cultures, with normoxic/anoxic cytotoxicity IC50 ratios up to 203. In vivo, the antitumor effects of CP-506 were selective for hypoxic tumor cells and causally related to tumor oxygenation. CP-506 effectively decreased the hypoxic fraction of tumor models and inhibited the growth of a wide range of xenografts, but only when hypoxia was present. A multivariate regression analysis revealed baseline tumor hypoxia and in vitro sensitivity to CP-506 significantly correlate with treatment response. Taken together, our results demonstrate the hypoxia-specific cytotoxicity and broad antitumor effects of CP-506.

##In this paper: 
To test whether HF and sensitivity of tumor cell lines to CP-506 defined as the anoxic IC50 value influenced ER or SGD, a database was constructed containing 434 animals (Supplementary Table 4). Only models for which HF and ER were available were included resulting in 381 observations. Outcome data were fit to multiple linear regression with different independent variables: HF (%), anoxic IC50 (?M), absolute cumulative exposure (mg), and mean volume doubling time (VDT, days). Absolute cumulative exposure (further referred to as CP-506 dose) was calculated as the relative dose (mg/kg) times the total number of injections the animal receives times the bodyweight (kg) of the animal at the start of treatment. Since HF could not be determined in the same animal as an outcome variable but was rather determined in sentinel animals in an accompanying experiment, HF was drawn from normal distribution with mean HF and SD estimated for each tumor type in parallel (vide supra). Then, HF was randomly assigned to each tumor and a multivariate linear regression analysis was performed (STATA/IC 11.1). This procedure was repeated 500 times. The number of significant associations was counted and the direction of the effect (positive or negative) was recorded for each parameter in the model. The nomogram to predict ER was obtained using R (version 4.0.2).

### Code (do file) for simulations is STATA11:
Bootstrapping HF to fit multiple linear regression
. cd C:\01Ala\CP-506\modelling
. do simul
//*
	tempname sim
	postfile `sim' p_dose coef_dose p_hf coef_hf p_ic50anox 	coef_ic50anox p_vdt coef_vdt using resout,replace
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
	replace hf_sim=`hf1' + `sd1'* invnorm(uniform()) if 	tu_id==1
	replace hf_sim=`hf2' + `sd2'* invnorm(uniform()) if 	tu_id==2
	replace hf_sim=`hf3' + `sd3'* invnorm(uniform()) if 	tu_id==3
	replace hf_sim=`hf4' + `sd4'* invnorm(uniform()) if 	tu_id==4
	replace hf_sim=`hf5' + `sd5'* invnorm(uniform()) if 	tu_id==5
	replace hf_sim=`hf8' + `sd8'* invnorm(uniform()) if 	tu_id==8
	replace hf_sim=`hf9' + `sd9'* invnorm(uniform()) if 	tu_id==9
	replace hf_sim=`hf10' + `sd10'* invnorm(uniform()) if 	tu_id==10
	replace hf_sim=`hf11' + `sd11'* invnorm(uniform()) if 	tu_id==11
	replace hf_sim=`hf12' + `sd12'* invnorm(uniform()) if 	tu_id==12
	replace hf_sim=`hf13' + `sd13'* invnorm(uniform()) if 	tu_id==13
	replace hf_sim=`hf14' + `sd14'* invnorm(uniform()) if 	tu_id==14
	replace hf_sim=`hf15' + `sd15'* invnorm(uniform()) if 	tu_id==15
	replace hf_sim=`hf16' + `sd16'* invnorm(uniform()) if 	tu_id==16
	replace hf_sim=`hf17' + `sd17'* invnorm(uniform()) if 	tu_id==17
	replace hf_sim=`hf18' + `sd18'* invnorm(uniform()) if 	tu_id==18
	replace hf_sim=`hf19' + `sd19'* invnorm(uniform()) if 	tu_id==19
	replace hf_sim=`hf20' + `sd20'* invnorm(uniform()) if 	tu_id==20
	
	replace hf_sim=0 if hf_simul<0
	regress ER2 abs_cum_dose hf_sim ic50anox vdt
	
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
	
	post `sim' (p_dose) (coef_dose) (p_hf) (coef_hf) 	(p_ic50anox) (coef_ic50anox) (p_vdt) (coef_vdt)
	
	}

	postclose `sim'

*//
. tempname sim
. postfile `sim' p_dose coef_dose p_hf coef_hf p_ic50anox coef_ic50anox p_vdt coef_vdt using resout,replace
. use datacp506-nomissingER.dta, clear
. gen hf_simul=.
(400 missing values generated)
. set output error
end of do-file
. use resout, clear
. save results_bs, replace
file results_bs.dta saved
