clc, clear all
% Stoichiometric parameters
YH=0.758; % Yield for heterotrophic biomass
YA=0.24; % Yield for autotrophic biomass
frXI=0.08; %Fraction of biomass yielding to participate products
iNBM=0.068;% g Ng?1 COD Mass of biomass per mass of COD in biomass
iNXI=0.06; %g Ng?1 COD Mass of biomass per mass of COD in products formed by biomass decay

%Kinetic parameters
uH=4.0; %day?1 Maximum specific growth rate for heterotrophic biomass
bH= 0.94; %day?1 Decay rate coefficient for heterotrophic biomass
KS =10.0; %g CODm?3 Half-saturation coefficient for heterotrophic biomass
KOH= 0.2; %gO2 m?3 Oxygen half-saturation coefficient for heterotrophic biomass
KNO=0.5; %gNm?3 Nitrate half-saturation coefficient for denitrifying heterotrophic biomass
uA =0.5; %day?1 Maximum specific growth rate for autotrophic biomass
bA= 0.05; %day?1 Decay rate coefficient for autotrophic biomass
KNHA= 1.0; %gNm?3 Ammonia half-saturation coefficient for autotrophic biomass
KOA= 0.4; %gO2 m?3 Oxygen half-saturation coefficient for autotrophic biomass
nNOg= 0.8; %Correction factor for ?H under anoxic conditions
nNOh= 0.8; %Correction factor for ?h under anoxic conditions
kh =3.0; %day?1 Maximum specific hydrolysis rate
KX =0.1; %Half-saturation coefficient for hydrolysis of slowly biodegradable substrate
ka =0.05; %m3 g?1 CODday?1 Ammonification rate

%Settling parameters
fns=2.45*10^(-3); %Non-settable fraction of the effluent suspended solids


%Component conc. in feed and reactor- REPLACE WITH REAL DATA
xin=[343, 300, 350, 300, 285, 315, 320, 360, 344, 318, 333].*1000 ; %mg/m3 %Component concentration vector infeed
xat=[202, 200, 150, 116, 129, 197, 156, 144, 133, 140, 202].*1000; %mg/m3 %Component conc. vector in reactor
%Unknowns- MEGH HELP -at=aeration tank
Ssat=9999; 
Soat=9999;
XBHat=9999;
XBAat=9999;
XSat=9999;
XNDat=9999;
SNOat=9999;
SNHat=9999;
SNDat=9999;

%Kinetic rates
p1=(uH*Ssat*Soat*XBHat)/((KS+Ssat)*(KOH+Soat));
p2=(uH*Ssat*KOH*SNOat*nNOg*XBHat)/((KS+Ssat)*(KOH+Soat)*(KNO+SNOat));
p3=(uA*SNHat*Soat*XBAat)/((KNHA+SNHat)*(KOA+Soat));
p4=bH*XBHat;
p5=bA*XBAat;
p6=ka*SNDat*XBHat;
p7=((XBHat*kh*(XSat/XBHat))/(KX+(XSat/XBHat)))*((Soat/(KOH+Soat))...
    +(nNOh*KOH*SNOat)/((KOH+Soat)*(KNO+SNOat)));
p8=((XBHat*kh*(XNDat/XBHat))/(KX+(XSat/XBHat)))*((Soat/(KOH+Soat))...
    +(nNOh*KOH*SNOat)/((KOH+Soat)*(KNO+SNOat))); 

%Model state variables
Si=0;%g COD/m3
Ss=(-1/YH)*(p1+p2)+p7; %g COD/m3
Xi=frXI*(p4+p5);%g COD/m3
Xs=((1-frXI)*(p4+p5)) - p7;%g COD/m3
Xbh=p1+p2 -p4; %g COD/m3
Xba= p3-p5;%g COD/m3
Sno= ((-1*p2)*((1-YH)/(2.86*YH))) + (p3/YA); %gN/m3
Snh= (-1*iNBM*(p1+p2))-((iNBM+(1/YA))*p3) + p6;%gN/m3
Snd= -p6 + p8;%gN/m3
Xnd=((iNBM-(frXI*iNXI))*(p4+p5)) - p8; %gN/m3
So=((-1*p1)*((1-YH)/YH))- ((p3*(4.57-YA))/YA);%gO2/m3

%Apparent Reaction rates r
r=[Si,Ss, Xi, Xs, Xbh, Xba, Sno, Snh, Snd, Xnd, So];

% Process Configuration
%Dry Weather Conditions
Qin=3800; %m3/d
CODin = 343; %mg/L
TNin=33; %mg/L after primary treatment
Vat = 2050; %m3 volume of aeration tank
nt = 3; %# of turbines - Mechanical surface aerators
turbinep = 30; %kW power in each turbine
power = nt*turbinep; %kW
kLa=4.5; %h^-1 oxygen rate coefficient
Aset= 855; %m2 %Setling tank area
Hset = 2.8; %m %Settling tank depth
Qrs = 7600; %m3/d recirculation rate of solids to aeration tank via settlign tank
Qws = 75; %m3/d removed from system via settling tank
Nc =11; % # of cycles - does not vary day to day. DEFAULT AERATION
atime = 12.75; %h/d of cumulted aeration time

%AAS process model ASM1 - Activated Sludge Model No. 1 ASSUMPTIONS
%Total alkalinity not included
%CSTR Reactor assumption - completely stirred
%inert particulate material from influent and biomass decay = Xi
%Particulate material in effluent = fns*concentration in bioreactor
%Sosat = saturation conc. of dissolved oxygen
%So = dissolved oxygen conc.
%kLa = coefficient of oxygen transfer
%f(1) & f(2)= ODES for on-off sequence respectively 

Din=Qin/Vat; %d^-1
Drs=Qrs/Vat; %d^-1
Dws=Qws/Vat; %d^-1
gamma= ((1-fns)*(Din-Dws))/(Drs+Dws);
sc=[1,2,7,8,9]; %soluble components
pc=[3,4,5,6,10]; %particulate components
doi=[11]; %dissolved oxygen
f1=[]; f2=[]; f111=[]; f112=[];
for i=1:Nc 
    if length(find(sc==i))==1
        f1=[f1, Din*(xin(i)-xat(i)) + r(i)];
        f2=[f2, Din*(xin(i)-xat(i)) + r(i)];
    elseif length(find(pc==i))==1
        f1=[f1, Din*(xin(i)-xat(i)) + gamma*Drs*xat(i)+ r(i)];
        f2=[f2, Din*(xin(i)-xat(i)) + gamma*Drs*xat(i)+ r(i)];
    elseif length(find(doi==i))==1
        Sosat=0.6; %MEGH CHANGE THESE
        Ao = kLa*(Sosat-So);
        f111=[f111, Din*(xin(i)-xat(i)) +r(i) + Ao]; %air on
        f112=[f112, Din*(xin(i)-xat(i)) +r(i)]; %air off
end
end

%OUTPUTS
%f1, f2, f111, f112