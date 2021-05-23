function [pout,To4,mdot_F,gamma_g]=Combustionchamber(poin,delta_p_b,mdot_A,TAin,Tcc,eff_cc,Cpg,Ro)

To4=Tcc;

a=10;
b=22;

W_0=16.0e-3;
W_C=12.0e-3;
W_H=1.0e-3;
W_N=14.0067e-3;
W_C02=W_C+2*W_0; %kg/mol
W_H2O=2*W_H+W_0;
W_02=2*W_0; %Molecular mass 
W_N2=2*W_N;
W_F=a*W_C+b*W_H;
W_air=(W_02+3.76*W_N2)/(1+3.76);
R_air=Ro/W_air;
gamma_a=Cpa/(Cpa-R_air);

R_H20=Ro/W_H2O;
R_C02=Ro/W_C02;
R_02=Ro/W_02;
R_N2=Ro/W_N2;
R_F=Ro/W_F;

fst=(a+b/4)*(W_02+3.76*W_N2)/W_F;

% Cp_CO2=R*(2.401+8.735e3*T-6.607e6*T^2+2.002e9*T^3);
% Cp_H20=R*(4.070-1.108e3*T+4.152e6*T^2-2.964e9*T^3+0.807e12*T^4);
% Cp_02=R*(3.626-1.878e3*T+7.055e6*T^2-6.764e9*T^3+2.156e12*T^4);
% Cp_N2=R*(3.675-1.208e3*T+2.324e6*T^2-0.632e9*T^3-0.226e12*T^4);

%Tcc
% h_C02=-393.520 + Int_Cp_CO2;
% h_H20=-241.820 + Int_Cp_H20;
% h_02=0 + Int_Cp_02;
% h_N2=0 + Int_Cp_N2;
% h_02_in=0 + Int_Cp_02_in;
% h_N2_in=0 + Int_Cp_N2_in;
h_F_in=-294366;

CoeffhCO2L=[0.2356773520e+01 0.8984596770e-02 -0.7123562690e-05 0.2459190220e-08 -0.1436995480e-12 -0.4837196970e+05];
CoeffhCO2H=[0.3857460290e+01 0.4414370260e-02 -0.2214814040e-05 0.5234901880e-09 -0.4720841640e-13 -0.4875916600e+05];
CoeffhH20L=[0.4198640560e+01 -0.2036434100e-02 0.6520402110e-05 -0.5487970620e-08 0.1771978170e-11 -0.3029372670e+05];
CoeffhH20H=[0.3033992490e+01 0.2176918040e-02 -0.1640725180e-06 -0.9704198700e-10 0.1682009920e-13 -0.3000429710e+05];
CoeffhN2L=[0.3298677000e+01 0.1408240400e-02 -0.3963222000e-05 0.5641515000e-08 -0.2444854000e-11 -0.1020899900e+04];
CoeffhN2H=[0.2926640000e+01 0.1487976800e-02 -0.5684760000e-06 0.1009703800e-09 -0.6753351000e-14 -0.9227977000e+03];
Coeffh02L=[0.3782456360e+01 -0.2996734160e-02 0.9847302010e-05 -0.9681295090e-08 0.3243728370e-11 -0.1063943560e+04];
Coeffh02H=[0.3282537840e+01 0.1483087540e-02 -0.7579666690e-06 0.2094705550e-09 -0.2167177940e-13 -0.1088457720e+04];

if Tain<1000
h_02_in=R_02*TAin*(Coeffh02L(1)+(Coeffh02L(2)/2)*TAin+(Coeffh02L(3)/3)*TAin^2+(Coeffh02L(4)/4)*TAin^3+(Coeffh02L(5)/5)*TAin^4+(Coeffh02L(6))/TAin);
h_N2_in=R_N2*TAin*(CoeffhN2L(1)+(CoeffhN2L(2)/2)*TAin+(CoeffhN2L(3)/3)*TAin^2+(CoeffhN2L(4)/4)*TAin^3+(CoeffhN2L(5)/5)*TAin^4+(CoeffhN2L(6))/TAin);
end

if Tain>1000
h_02_in=R_02*TAin*(Coeffh02H(1)+(Coeffh02H(2)/2)*TAin+(Coeffh02H(3)/3)*TAin^2+(Coeffh02H(4)/4)*TAin^3+(Coeffh02H(5)/5)*TAin^4+(Coeffh02H(6))/TAin);
h_N2_in=R_N2*TAin*(CoeffhN2H(1)+(CoeffhN2H(2)/2)*TAin+(CoeffhN2H(3)/3)*TAin^2+(CoeffhN2H(4)/4)*TAin^3+(CoeffhN2H(5)/5)*TAin^4+(CoeffhN2H(6))/TAin);
end

if Tcc < 1000
    h_CO2=R_C02*Tcc*(CoeffhCO2L(1)+(CoeffhCO2L(2)/2)*Tcc+(CoeffhCO2L(3)/3)*Tcc^2+(CoeffhCO2L(4)/4)*Tcc^3+(CoeffhCO2L(5)/5)*Tcc^4+(CoeffhCO2L(6))/Tcc);
    h_H20=R_H20*Tcc*(CoeffhH20L(1)+(CoeffhH20L(2)/2)*Tcc+(CoeffhH20L(3)/3)*Tcc^2+(CoeffhH20L(4)/4)*Tcc^3+(CoeffhH20L(5)/5)*Tcc^4+(CoeffhH20L(6))/Tcc);
    h_N2=R_N2*Tcc*(CoeffhN2L(1)+(CoeffhN2L(2)/2)*Tcc+(CoeffhN2L(3)/3)*Tcc^2+(CoeffhN2L(4)/4)*Tcc^3+(CoeffhN2L(5)/5)*Tcc^4+(CoeffhN2L(6))/Tcc);
    h_02=R_02*Tcc*(Coeffh02L(1)+(Coeffh02L(2)/2)*Tcc+(Coeffh02L(3)/3)*Tcc^2+(Coeffh02L(4)/4)*Tcc^3+(Coeffh02L(5)/5)*Tcc^4+(Coeffh02L(6))/Tcc);
end
if Tcc > 1000
    h_CO2=R_C02*Tcc*(CoeffhCO2H(1)+(CoeffhCO2H(2)/2)*Tcc+(CoeffhCO2H(3)/3)*Tcc^2+(CoeffhCO2H(4)/4)*Tcc^3+(CoeffhCO2H(5)/5)*Tcc^4+(CoeffhCO2H(6))/Tcc);
    h_H20=R_H20*Tcc*(CoeffhH20H(1)+(CoeffhH20H(2)/2)*Tcc+(CoeffhH20H(3)/3)*Tcc^2+(CoeffhH20H(4)/4)*Tcc^3+(CoeffhH20H(5)/5)*Tcc^4+(CoeffhH20H(6))/Tcc);
    h_N2=R_N2*Tcc*(CoeffhN2H(1)+(CoeffhN2H(2)/2)*Tcc+(CoeffhN2H(3)/3)*Tcc^2+(CoeffhN2H(4)/4)*Tcc^3+(CoeffhN2H(5)/5)*Tcc^4+(CoeffhN2H(6))/Tcc);
    h_02=R_02*Tcc*(Coeffh02H(1)+(Coeffh02H(2)/2)*Tcc+(Coeffh02H(3)/3)*Tcc^2+(Coeffh02H(4)/4)*Tcc^3+(Coeffh02H(5)/5)*Tcc^4+(Coeffh02H(6))/Tcc);
end    
    
Lambda=(h_F_in-a*h_CO2-(b/2)*h_H20 + (a+b/4)*h_02)/((a+b/4)*((h_02-h_02_in)+3.76*(h_N2-h_N2_in)));
pout=poin-delta_p_b;

%Actual fuel air ratio
W_output=((a*W_CO2+(b/2)*W_H20+(Lambda-1)*(a+b/4)*W_02+3.76*Lambda*(a+b/4)*W_N2)/(a+b/2+(Lambda-1)*(a+b/4)+3.76*Lambda*(a+b/4)));
R_g=Ro/W_output;
gamma_g=Cpg/(Cpg-R_g);

fth=1/(Lambda*fst);

f=fth/eff_cc;

mdot_F=f*mdot_A;   