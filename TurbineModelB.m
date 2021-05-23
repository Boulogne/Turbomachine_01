function [po4,Wt,To4]=TurbineModelB(Wc,po3,To3,eff_pt,eff_m,gamma_g)
Wt=Wc/eff_m;
To4=To3-Wt/Cpg;
po4=po3*(To3/To4)^(-(gamma_g)/(eff_pt*(gamma_g-1)));
