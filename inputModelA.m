function [To1,po1]=inputModelA(pext,Text,Vext,Stagnationpratio,Cp,gamma_a)

Toext = Text + (Vext^2)/(2*Cp);
poext = pext*(Toext/Text)^(gamma_a/(gamma_a-1));
To1 = Toext;
po1 = poext*Stagnationpratio;
