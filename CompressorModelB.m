function [To2,Wc,po2]=CompressorModelB(po1,To1,piM_c,eff_pc,mdot_c,gamma_a,Cpa)
po2 = po1*piM_c;
To2=To1*(po2/po1)^((gamma_a-1)/(eff_pc*gamma_a));

Wc=mdot_c*Cpa*(To2-To1);