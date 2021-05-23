%Exercise1
% Input sata 
pa=1e5 ; %[Pa]
Ta=298.15; %[K]
To3=1600; %[K]
Fuel ='C10H22(l)';
Tfin = 25 ; %[°C]
% Isentropic overall efficiencies
eff_c=0.85;
eff_t=0.9;
% Mechanical efficiency 
eff_m=0.98;
%Combustion efficiency
eff_cc=0.97;
% Heat exchanger (recuparator) efficiency 
eff_hx=0.85;

delta_p_53=0.06*1e5; %[Pa]
delta_p_o25=0.; %[Pa] % Condition isobar
delta_p_o46=0.; %[Pa] % Condition isobar

Cpa=1005.;
Cpg=1150.;

Fuel ='C10H22(l)';
a=10; % Number of carbon
b=22; % Number of Hydrogen 
% n-Decane 
Cp=296;             %[J/molK]
LHV=44.17e6;   %[J/kg]
H0=-294366;     %[J/mol]

W_0=16.0e-3;%[kg/mol]
W_C=12.0e-3;%[kg/mol]
W_H=1.0e-3;%[kg/mol]
W_N=14.0067e-3;%[kg/mol]
W_C02=W_C+2*W_0;  %[kg/mol]
W_H2O=2*W_H+W_0;%[kg/mol]
W_02=2*W_0; %[kg/mol] 
W_N2=2*W_N;%[kg/mol]
W_F=a*W_C+b*W_H;%[kg/mol]
W_air=(W_02+3.76*W_N2)/(1+3.76);%[kg/mol]
R_air=Ro/W_air;%[kg/mol]
gamma_a=Cpa/(Cpa-R_air);
%mass flow rate 
mdot=220;

Sin=3; %[m^2] SA
S1=3.75;%[m^2] S1

Po1=pa;
To1=Ta;


%% Comp 1-2

% Model A
Po2 = Po1*piM_c;
To2S=To1*(Po2/Po1)^((gamma_a-1)/(gamma_a));
To2=To1+(To2S-To1)/eff_c;
wc=mdot*Cpa*(To2-To1);

%% Without Heat exchanger entry 2 - 5
Po5=Po2 ;
To5=To2;

%% Combustion chamber 5 - 3
%[Po4,To4,mdot_F,gamma_g]=Combustionchamber(Po3,delta_p_b,mdot_A,TAin,Tcc,eff_cc,Cpg,Ro)
Tcc=To3;
TAin=Ta;
poin=Po5;
delta_p_b=delta_p_53;
a=10;
b=22;

[Po3,To3,mdot_F,gamma_g,mdot_g,R_g,f]=Combustionchamber(Po5,delta_p_b,TAin,Tcc,eff_cc,Cpg,Ro,Cpa,mdot)


%% Turbine 3 - 4
%[po5,Wtc,To5]=TurbineModelB(Wc,po4,To4,eff_pt,eff_m);
wt=wc/eff_m;
To4=To3-wt/(mdot_g*Cpg);
To4S=To3-(To3-To4)/eff_t;
Po4=Po3*(To3/To4S)^(-(gamma_g)/(gamma_g-1));

%%
Po6=Po4-delta_p_o46;
To6=To4

%% Output Postprocessus

fprintf('Stagnation conditions for the turbofan\n ')
fprintf('              Tt(K)                    Pt(Pa)\n ')
fprintf('1 =       % 6.2f              %6.3f \n',To1,Po1*1e-5);
fprintf('2 =       % 6.2f              %6.3f \n',To2,Po2*1e-5);
fprintf('5 =       % 6.2f              %6.3f \n',To5,Po5*1e-5);
fprintf('3 =       % 6.2f              %6.3f \n',To3,Po3*1e-5);
fprintf('4 =       % 6.2f              %6.3f \n',To4,Po4*1e-5);
fprintf('6 =       % 6.2f              %6.3f \n',To6,Po6*1e-5);
SFC=mdot_F/Fs;%[kg/(N.s)]
SFCg=mdot_F*1e6/Fs; %[g/(kN.s)]
fprintf('\n Specific fuel consumption  SFC = % 6.3e  kg/(N.s)  \n',SFC);
fprintf(' Specific fuel consumption  SFC = % 6.3f  kg/(N.h)  \n',3600*SFC);
fprintf(' Specific fuel consumption  SFC = % 6.3f  g/(kN.s)  \n',SFCg);
%% Air properties
function [R_air,gamma_a] = propertyair(Cpa,Ro)
%PROPERTYAIR Summary of this function goes here

W_0=16.0e-3;%[kg/mol]
W_C=12.0e-3;%[kg/mol]
W_H=1.0e-3;%[kg/mol]
W_N=14.0067e-3;%[kg/mol]
W_02=2*W_0; %[kg/mol] 
W_N2=2*W_N;%[kg/mol]
W_air=(W_02+3.76*W_N2)/(1+3.76);%[kg/mol]
R_air=Ro/W_air;%[kg/mol]
gamma_a=Cpa/(Cpa-R_air);
end

%% Input model B
function [To1,Po1,Toext,Poext,Toin,Poin,V1,Vin] = InputmodelB(Ta,Pa,Va,Cpa,gamma_a,R_air,mdot,Sin,S1,eff_d)

Toext=Ta+Va^2/(2*Cpa);
Poext=Pa*(Toext/Ta)^(gamma_a/(gamma_a-1));
%Exterior to duct entry (ext - in)
Toin=Toext;
Poin=Poext;
diff_Vin=1;
mindiff=1e-6;
Vin=250;

 while diff_Vin>mindiff
     Tin=Toin-(Vin^2)/(2*Cpa);
     Pin=Poin*(Tin/Toin)^(gamma_a/(gamma_a-1));
     rho_in=Pin/(R_air*Tin);
     VinN=mdot/(rho_in*Sin);
     diff_Vin=abs(VinN-Vin); 
     Vin=VinN;
 end

%Duct entry to engine face (in-1)
To1=Toin;
To1S=Tin+eff_d*(To1-Tin);
Po1=Pin*(To1S/Tin)^(gamma_a/(gamma_a-1));
diff_V1=1;
V1=Vin;

while diff_V1>mindiff
    T1=To1-(V1)^2/(2*Cpa);
    P1=Po1*(T1/To1)^(gamma_a/(gamma_a-1));
    rho1=P1/(R_air*T1); 
    V1N=mdot/(rho1*S1);
    diff_V1=abs(V1N-V1);
    V1=V1N;
end

end

%% Compressor model B
function [To_out,Wc,po_out]=CompressorModelB(po_in,To_in,piM_c,eff_pc,mdot,gamma,Cp)
po_out = po_in*piM_c;
To_out=To_in*(po_out/po_in)^((gamma-1)/(eff_pc*gamma));
Wc=mdot*Cp*(To_out-To_in);
end
%% Combustion chamber 

function [Po_out,To_out,mdot_F,gamma_g,mdot_g,R_g,f]=Combustionchamber(po_in,delta_p_b,TAin,Tcc,eff_cc,Cpg,Ro,Cpa,mdot)

To_out=Tcc;
Po_out=po_in-delta_p_b;

a=10;
b=22;
W_0=16.0e-3;%Molecular mass 
W_C=12.0e-3;
W_H=1.0e-3;
W_N=14.0067e-3;
W_CO2=W_C+2*W_0; %kg/mol
W_H20=2*W_H+W_0;
W_02=2*W_0; 
W_N2=2*W_N;
W_F=a*W_C+b*W_H;
W_air=(W_02+3.76*W_N2)/(1+3.76);
R_air=Ro/W_air;
gamma_a=Cpa/(Cpa-R_air);

  h_F_in=-294366;

 
CoeffhCO2L=[0.2356773520e+01 0.8984596770e-02 -0.7123562690e-05 0.2459190220e-08 -0.1436995480e-12 -0.4837196970e+05];
CoeffhCO2H=[0.3857460290e+01 0.4414370260e-02 -0.2214814040e-05 0.5234901880e-09 -0.4720841640e-13 -0.4875916600e+05];
CoeffhH20L=[0.4198640560e+01 -0.2036434100e-02 0.6520402110e-05 -0.5487970620e-08 0.1771978170e-11 -0.3029372670e+05];
CoeffhH20H=[0.3033992490e+01 0.2176918040e-02 -0.1640725180e-06 -0.9704198700e-10 0.1682009920e-13 -0.3000429710e+05];
CoeffhN2L=[0.3298677000e+01 0.1408240400e-02 -0.3963222000e-05 0.5641515000e-08 -0.2444854000e-11 -0.1020899900e+04];
CoeffhN2H=[0.2926640000e+01 0.1487976800e-02 -0.5684760000e-06 0.1009703800e-09 -0.6753351000e-14 -0.9227977000e+03];
Coeffh02L=[0.3782456360e+01 -0.2996734160e-02 0.9847302010e-05 -0.9681295090e-08 0.3243728370e-11 -0.1063943560e+04];
Coeffh02H=[0.3282537840e+01 0.1483087540e-02 -0.7579666690e-06 0.2094705550e-09 -0.2167177940e-13 -0.1088457720e+04];

if TAin<1000
h_02_in=Ro*TAin*(Coeffh02L(1)+(Coeffh02L(2)/2)*TAin+(Coeffh02L(3)/3)*TAin^2+(Coeffh02L(4)/4)*TAin^3+(Coeffh02L(5)/5)*TAin^4+(Coeffh02L(6))/TAin);
h_N2_in=Ro*TAin*(CoeffhN2L(1)+(CoeffhN2L(2)/2)*TAin+(CoeffhN2L(3)/3)*TAin^2+(CoeffhN2L(4)/4)*TAin^3+(CoeffhN2L(5)/5)*TAin^4+(CoeffhN2L(6))/TAin);
end

if TAin>1000
h_02_in=Ro*TAin*(Coeffh02H(1)+(Coeffh02H(2)/2)*TAin+(Coeffh02H(3)/3)*TAin^2+(Coeffh02H(4)/4)*TAin^3+(Coeffh02H(5)/5)*TAin^4+(Coeffh02H(6))/TAin);
h_N2_in=Ro*TAin*(CoeffhN2H(1)+(CoeffhN2H(2)/2)*TAin+(CoeffhN2H(3)/3)*TAin^2+(CoeffhN2H(4)/4)*TAin^3+(CoeffhN2H(5)/5)*TAin^4+(CoeffhN2H(6))/TAin);
end

if Tcc < 1000
    h_CO2=Ro*Tcc*(CoeffhCO2L(1)+(CoeffhCO2L(2)/2)*Tcc+(CoeffhCO2L(3)/3)*Tcc^2+(CoeffhCO2L(4)/4)*Tcc^3+(CoeffhCO2L(5)/5)*Tcc^4+(CoeffhCO2L(6))/Tcc);
    h_H20=Ro*Tcc*(CoeffhH20L(1)+(CoeffhH20L(2)/2)*Tcc+(CoeffhH20L(3)/3)*Tcc^2+(CoeffhH20L(4)/4)*Tcc^3+(CoeffhH20L(5)/5)*Tcc^4+(CoeffhH20L(6))/Tcc);
    h_N2=Ro*Tcc*(CoeffhN2L(1)+(CoeffhN2L(2)/2)*Tcc+(CoeffhN2L(3)/3)*Tcc^2+(CoeffhN2L(4)/4)*Tcc^3+(CoeffhN2L(5)/5)*Tcc^4+(CoeffhN2L(6))/Tcc);
    h_02=Ro*Tcc*(Coeffh02L(1)+(Coeffh02L(2)/2)*Tcc+(Coeffh02L(3)/3)*Tcc^2+(Coeffh02L(4)/4)*Tcc^3+(Coeffh02L(5)/5)*Tcc^4+(Coeffh02L(6))/Tcc);
end
if Tcc > 1000
    h_CO2=Ro*Tcc*(CoeffhCO2H(1)+(CoeffhCO2H(2)/2)*Tcc+(CoeffhCO2H(3)/3)*Tcc^2+(CoeffhCO2H(4)/4)*Tcc^3+(CoeffhCO2H(5)/5)*Tcc^4+(CoeffhCO2H(6))/Tcc);
    h_H20=Ro*Tcc*(CoeffhH20H(1)+(CoeffhH20H(2)/2)*Tcc+(CoeffhH20H(3)/3)*Tcc^2+(CoeffhH20H(4)/4)*Tcc^3+(CoeffhH20H(5)/5)*Tcc^4+(CoeffhH20H(6))/Tcc);
    h_N2=Ro*Tcc*(CoeffhN2H(1)+(CoeffhN2H(2)/2)*Tcc+(CoeffhN2H(3)/3)*Tcc^2+(CoeffhN2H(4)/4)*Tcc^3+(CoeffhN2H(5)/5)*Tcc^4+(CoeffhN2H(6))/Tcc);
    h_02=Ro*Tcc*(Coeffh02H(1)+(Coeffh02H(2)/2)*Tcc+(Coeffh02H(3)/3)*Tcc^2+(Coeffh02H(4)/4)*Tcc^3+(Coeffh02H(5)/5)*Tcc^4+(Coeffh02H(6))/Tcc);
end

Lambda=(h_F_in-a*h_CO2-(b/2)*h_H20 + (a+(b/4))*h_02)/((a+(b/4))*((h_02-h_02_in)+3.76*(h_N2-h_N2_in)));

fst=(a+b/4)*(W_02+3.76*W_N2)/W_F;
fth=1/(Lambda*fst);
f=fth/eff_cc;

%Actual fuel air ratio
W_output=((a*W_CO2+(b/2)*W_H20+(Lambda-1)*(a+b/4)*W_02+3.76*Lambda*(a+b/4)*W_N2)/(a+b/2+(Lambda-1)*(a+b/4)+3.76*Lambda*(a+b/4)));
R_g=Ro/W_output;
gamma_g=Cpg/(Cpg-R_g);

mdot_F=f*mdot;

mdot_g=mdot+mdot_F;

end

%% Turbine Model B
function [po_out,Wt,To_out]=TurbineModelB(Wc,po_in,To_in,eff_pt,eff_m,gamma,Cp,mdot)
Wt=Wc/(eff_m);
To_out=To_in-Wt/(mdot*Cp);
po_out=po_in*(To_in/To_out)^(-(gamma)/(eff_pt*(gamma-1)));
end

%% Propelling nozzle
function [To_out,Po_out,A_out,V_out,P_out,T_out,M_out] = Propellantnozzle(To_in,Po_in,gamma,R,Cp,eff_nc,mdot,Pa)
%PROPELLANT NOZZLE Summary of this function goes here
To_out=To_in;
Po_out=Po_in;

M_out=1; % Chocked flow / Critical flow
T_out=To_out/(1+gamma*R/(2*Cp));
V_out=M_out*sqrt(gamma*R*T_out);
TS_out=To_out*(1-1/eff_nc)+T_out/eff_nc;
Pcr_out=Po_out*(TS_out/To_in)^(gamma/(gamma-1)); 
P_out=Pcr_out;

if Pcr_out<Pa   %Pa % Then the flow is unchocked  
    P_out=Pa;
    TS_out=To_out*(P_out/Po_in)^(-gamma/(gamma-1));
    T_out=To_out-(To_out-TS_out)*eff_nc;
    V_out=sqrt((To_out-T_out)*2*Cp);
    M_out=V_out/sqrt(gamma*R*T_out);
end
To_out =T_out * (1 + (gamma- 1)/2 * M_out^2); 
Po_out = P_out * (1 + (gamma - 1)/2 * M_out^2)^(gamma/(gamma-1));

rho_out=P_out/(R*T_out);
A_out=mdot/(rho_out*V_out);
end
