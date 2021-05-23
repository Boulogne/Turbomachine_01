% Exercise 2: Twin-Spool Turbofan 
clear;
clc;
%input data
Ro=8.314462618;
Pa=0.3*1e5 ;                        %[Pa]
Ta=-45+273.15;                   %[K]
Va=250;                                %[m/s]
Stagnationpratio=0.99;
piM_f=1.8;
piM_c=16;
B=5;                                     % Bypass ratio 
To4=1600.;   %Combustion temperature         %[K]
T_fin =273.16+25;                                             %[K]


eff_pf=0.9;     %Fan polytropic efficiency
eff_pc=0.9;    %Compressor polytropic efficiency
eff_pt=0.92;   %Turbine polytropic efficiency
eff_d=0.9;      %
eff_nc=0.96;%Discharge isentropic efficiency for the core and bypass air (total to static)
eff_nb=eff_nc;
eff_m=0.97;     %Mechanical efficiency (for the two spools)
eff_cc=0.97;    %Combustion efficiency 

delta_p34=0.05; %[bar]
%total mass flow rate 
mdot=220.; %[kg/s]

Cpa=1005.;   %Heat capacity for the air   %[J/kgK]
Cpg=1150.;    %Heat capacity for the gaz  %[J/kgK]

Fuel ='C10H22(l)';
a=10; % Number of carbon
b=22; % Number of Hydrogen 
% n-Decane 
Cp=296;             %[J/molK]
LHVf=44.17e6;   %[J/kg]
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

mdot_h=mdot/(B+1); %[kg/s]
mdot_c=mdot-mdot_h;%[kg/s]
Va=250;
Sin=3; %[m^2] SA
S1=3.75;%[m^2] S1

% Intake (ext - 1 ) (stagnation pressure ratio as input data po1/poext)
% Model A 
%[To1,po1]=inputModelA(Pa,Ta,Va,Stagnationpratio,Cp,gamma_a);

%% INTAKE MODEL B
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



%% Fan 1-2
%[To2,Wf,Po2]=CompressorModelB(Po1,To1,piM_f,eff_pf,mdot_c,gamma_a,Cpa);

Po2 = Po1*piM_f;
To2=To1*(Po2/Po1)^((gamma_a-1)/(eff_pf*gamma_a));
Wf=mdot*Cpa*(To2-To1);

%% Comp 2-3
%[To3,Wc,Po3]=CompressorModelB(Po2(1),To2(1),piM_c,eff_pc,mdot_c,gamma_a,Cpa);

Po3 = Po2*piM_c;
To3=To2*(Po3/Po2)^((gamma_a-1)/(eff_pc*gamma_a));
Wc=mdot_c*Cpa*(To3-To2);


%% Combustion chamber 3 - 4
%[Po4,To4,mdot_F,gamma_g]=Combustionchamber(Po3,delta_p_b,mdot_A,TAin,Tcc,eff_cc,Cpg,Ro)
Tcc=To4;
Tain=Ta;
TAin=Ta;
poin=Po3;
delta_p_b=delta_p34;
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

R_H20=Ro/W_H2O;
R_C02=Ro/W_C02;
R_02=Ro/W_02;
R_N2=Ro/W_N2;
R_F=Ro/W_F;

fst=(a+b/4)*(W_02+3.76*W_N2)/W_F;


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
h_02_in=Ro*TAin*(Coeffh02L(1)+(Coeffh02L(2)/2)*TAin+(Coeffh02L(3)/3)*TAin^2+(Coeffh02L(4)/4)*TAin^3+(Coeffh02L(5)/5)*TAin^4+(Coeffh02L(6))/TAin);
h_N2_in=Ro*TAin*(CoeffhN2L(1)+(CoeffhN2L(2)/2)*TAin+(CoeffhN2L(3)/3)*TAin^2+(CoeffhN2L(4)/4)*TAin^3+(CoeffhN2L(5)/5)*TAin^4+(CoeffhN2L(6))/TAin);
end

if Tain>1000
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

Lambda=(h_F_in-a*h_CO2-(b/2)*h_H20 + (a+b/4)*h_02)/((a+b/4)*((h_02-h_02_in)+3.76*(h_N2-h_N2_in)));
Po4=poin-delta_p_b;
clc;  
%Actual fuel air ratio
W_output=((a*W_CO2+(b/2)*W_H20+(Lambda-1)*(a+b/4)*W_02+3.76*Lambda*(a+b/4)*W_N2)/(a+b/2+(Lambda-1)*(a+b/4)+3.76*Lambda*(a+b/4)));
R_g=Ro/W_output;
gamma_g=Cpg/(Cpg-R_g);
fst=(a+b/4)*(W_02+3.76*W_N2)/W_F;
fth=1/(fst);%*Lambda
f=fth/eff_cc;
Lambda=1/(fst*f);
mdot_F=f*mdot_c; 
% mdot_F=0.5;
mdot_g=mdot_c+mdot_F;

%% High pressure Turbine 4 - 5
%[po5,Wtc,To5]=TurbineModelB(Wc,po4,To4,eff_pt,eff_m);
Wt_HPT=Wc/eff_m;
To5=To4-Wt_HPT/(mdot_g*Cpg);
Po5=Po4*(To4/To5)^(-(gamma_g)/(eff_pt*(gamma_g-1)));


%% Low Pressure Turbine 5 - 6
%[po6,Wtf,To6]=TurbineModelB(Wf,po5,To5,eff_pt,eff_m);
Wt_LPT=Wf/eff_m;
To6=To5-Wt_LPT/(mdot_g*Cpg);
Po6=Po5*(To5/To6)^(-(gamma_g)/(eff_pt*(gamma_g-1)));

%% Propellant nozzle 6 - 7
To7=To6;
Po7=Po6;

M7=1; % Chocked flow / Critical flow
T7=To6/(1+gamma_g*R_g/2*Cpg);
V7=M7*sqrt(gamma_g*R_g*T7);
T7S=To6*(1-1/eff_nc)+T7/eff_nc;
P7cr=Po6*(T7S/To6)^(gamma_g/(gamma_g-1)); 
P7=P7cr;

if P7cr<Pa   %Pa % Then the flow is unchocked  
    P7=Pa;
    T7S=To6*(P7/Po6)^(gamma_g/(gamma_g-1));
    T7=To6-(To6-T7S)*eff_nc;
    V7=sqrt((To6-T7)*2*Cpg);
    M7=V7/sqrt(gamma_g*R_g*T7);
end
% To7 =T7 * (1 + (gamma_g- 1)/2 * M7^2); %issues
% Po7 = P7 * (1 + (gamma_g - 1)/2 * M7^2)^(gamma_g/(gamma_g-1));

rho7=P7/(R_g*T7);
A7=mdot_g/(rho7*V7);

%% NOZZLE - 8
To8=To2;
Po8=Po2;

M8=1; % Chocked flow / Critical flow
T8=To2/(1+gamma_a*R_air/2*Cpa);
V8=M8*sqrt(gamma_a*R_air*T8);
T8S=To2*(1-1/eff_nb)+T8/eff_nb;
P8cr=Po2*(T8S/To2)^(gamma_a/(gamma_a-1)); 
P8=P8cr;

if P8cr<Pa   %Pa % Then the flow is unchocked  
    P8=Pa;
    T8S=To2*(P8/Po2)^(gamma_a/(gamma_a-1));
    T8=To2-(To2-T8S)*eff_nb;
    V8=sqrt((To2-T8)*2*Cpa);
    M8=V8/sqrt(gamma_a*R_air*T8);
end
%  To8 =T8 * (1 + (gamma_a- 1)/2 * M8^2);  % Issues
%  Po8 = P8 * (1 + (gamma_a - 1)/2 * M8^2)^(gamma_a/(gamma_a-1));

rho8=P8/(R_air*T8);
A8=mdot_g/(rho8*V8);




%% Output Postprocessus


fprintf('Stagnation conditions for the turbofan\n ')
fprintf('              Tt(K)                    Pt(Bar)\n ')
fprintf('a =       % 6.2f              %6.2f \n',Toext,Poext*1e-5);
fprintf('A =       % 6.2f              %6.2f \n',Toin,Poin*1e-5);
fprintf('1 =       % 6.2f              %6.2f \n',To1,Po1*1e-5);
fprintf('2 =       % 6.2f              %6.2f \n',To2,Po2*1e-5);
fprintf('3 =       % 6.2f              %6.2f \n',To3,Po3*1e-5);
fprintf('4 =       % 6.2f              %6.2f \n',To4,Po4*1e-5);
fprintf('5 =       % 6.2f              %6.2f \n',To5,Po5*1e-5);
fprintf('6 =       % 6.2f              %6.2f \n',To6,Po6*1e-5);
fprintf('7 =       % 6.2f              %6.2f \n',To7,Po7*1e-5);
fprintf('8 =       % 6.2f              %6.2f \n',To8,Po8*1e-5);

%Total thrust 
Fs=(mdot_c+mdot_F)*V7-mdot*Va -Pa*(A8+A7)+P8*A8+P7*A7; %[N]
MassFlowVerification = mdot_h+mdot_g-mdot-mdot_F;
%Specific fuel consumption
SFC=mdot_F/Fs; %[kg/(N.s)]
SFCg=mdot_F*1e6/Fs; %[g/(kN.s)]
%Global efficiencies
%Thermal 
eff_th=((mdot_h*V8^2)/2+((mdot_c+mdot_F)*V7^2)/2-(mdot*Va^2)/2)/(mdot_F*LHVf);
%Propulsion
eff_p=Fs*Va/((mdot_h*V8^2)/2+((mdot_c+mdot_F)*V7^2)/2-(mdot*Va^2)/2);
% Overall efficiency 
eff_o=eff_th*eff_p;%Fs*Va/(mdot*LHVf);

fprintf('\n Tota thrust Fs = % 6.2f  N  \n',Fs);
fprintf(' Tota thrust Fs = % 6.2f  kN  \n',Fs*1e-3);
fprintf('\n Specific fuel consumption  SFC = % 6.2e  kg/(N.s)  \n',SFC);
fprintf(' Specific fuel consumption  SFC = % 6.2f  g/(kN.s)  \n',SFCg);
fprintf('\n\n Thermal efficiency      Propulsive efficiency       Overall efficiency  \n');
fprintf(' % 6.2f                                        %6.3f                              %6.3f               \n',eff_th,eff_p,eff_o);