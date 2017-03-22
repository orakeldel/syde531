function [ ft, ff, total_tn ] = asm1(Qin, COD, blower_conf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



Q_avg = 3050;       % cubic m per d
COD_avg = 343;      % mg/L
TIN_avg = 33;       % mg/L
Qin_f = @(t) Q_avg*Qin(t);
COD_f = @(t) COD_avg*COD(t);

V_reactor = 2050;   % cubic m
Q_rs = 7600;        % cubic m per day
Q_w = 75;           % cubic m per day

Q_factor = @(t) (Qin_f(t)+Q_rs)./(Q_rs + Q_w);

Din = @(t) Qin_f(t)/V_reactor;
Drs = Q_rs/V_reactor;
Dws = Q_w/V_reactor;

BIG_T = [];
BIG_F = [];

opts = odeset('RelTol',1e-7,'AbsTol',1e-7);

for i=1:length(blower_conf(:,1))
    
    status = blower_conf(i,3);
    start_t = blower_conf(i,1);
    end_t = blower_conf(i,2);
    tspan  = [start_t end_t];
    
    if i==1 
        y0 = [17.98 2.27 2120.15 79.55 2238.65 115.18 0.02 9.70 0.14 6.29 0]; % mg/L
    else
        y0 = BIG_F(end,:);
    end
    
    if status == 1
        [t, f] = ode45(@(t, f) react_blower_on(t,f, COD_f, TIN_avg, Din, Drs, Dws), tspan, y0, opts);
    else
        [t, f] = ode45(@(t, f) react_blower_off(t,f, COD_f, TIN_avg, Din, Drs, Dws), tspan, y0, opts);
    end
    
    BIG_T = [BIG_T; t];
    BIG_F = [BIG_F; f];
end

ft = BIG_T;
ff = BIG_F;
fns=2.45*10^(-3); %Non-settable fraction of the effluent suspended solids
iNBM=0.068; % g Ng?1 COD Mass of biomass per mass of COD in biomass
iNXI=0.06;  %g Ng?1 COD Mass of biomass per mass of COD in products formed by biomass decay


total_tn = BIG_F(:,11) + BIG_F(:,8) + BIG_F(:,9) + ...
            fns*(BIG_F(:, 10) + iNBM*(BIG_F(:, 6) + BIG_F(:, 5)) + iNXI*BIG_F(:, 3));

end

function [r,p] = calculate_rxn(f, t_at)

% reaction-rate parameters
% -------------------------

YH=0.758;   % Yield for heterotrophic biomass
YA=0.24;    % Yield for autotrophic biomass
frXI=0.08;  %Fraction of biomass yielding to participate products
iNBM=0.068; % g Ng?1 COD Mass of biomass per mass of COD in biomass
iNXI=0.06;  %g Ng?1 COD Mass of biomass per mass of COD in products formed by biomass decay

% Kinetic parameters
% -------------------------
uH=4.0;     %day?1 Maximum specific growth rate for heterotrophic biomass
bH= 0.94;   %day?1 Decay rate coefficient for heterotrophic biomass
KS =10.0;   %g CODm?3 Half-saturation coefficient for heterotrophic biomass
KOH= 0.2;   %gO2 m?3 Oxygen half-saturation coefficient for heterotrophic biomass
KNO=0.5;    %gNm?3 Nitrate half-saturation coefficient for denitrifying heterotrophic biomass
uA =0.5;    %day?1 Maximum specific growth rate for autotrophic biomass
bA= 0.05;   %day?1 Decay rate coefficient for autotrophic biomass
KNHA= 1.0;  %gNm?3 Ammonia half-saturation coefficient for autotrophic biomass
KOA= 0.4;   %gO2 m?3 Oxygen half-saturation coefficient for autotrophic biomass
nNOg= 0.8;  %Correction factor for ?H under anoxic conditions
nNOh= 0.8;  %Correction factor for ?h under anoxic conditions
kh =3.0;    %day?1 Maximum specific hydrolysis rate
KX =0.1;    %Half-saturation coefficient for hydrolysis of slowly biodegradable substrate
ka =0.05;   %m3 g?1 CODday?1 Ammonification rate
KNHH = 0.05;   %

% Settling parameters
% ------------------------
fns=2.45*10^(-3); %Non-settable fraction of the effluent suspended solids

% ODE model
% ----------------

Si = f(1);
Ss = f(2);
Xi = f(3);
Xs = f(4);
Xbh = f(5);
Xba = f(6);
Sno = f(7);
Snh = f(8);
Snd = f(9);
Xnd = f(10);
So = f(11);

% Base Reaction Matrix
% -------------------------
p = [];
p(1) = uH * (Ss)/(KS + Ss) * ((So)/(KOH + So)) * (Snh/(KNHH+Snh)) * Xbh;
p(2) = uH * ((Ss)/(KS + Ss)) * ((KOH)/(KOH + So))* (Snh/(KNHH+Snh)) * ((Sno)/(KNO + Sno)) * nNOg * Xbh;
p(3) = uA * ((Snh)/(KNHA + Snh)) * ((So)/(KOA + So)) * Xba;
p(4) = bH*Xbh;
p(5) = bA*Xba;
p(6) = ka*Snd*Xbh;
s_no_factor = (So/(KOH + So)) + nNOh * (KOH/(KOH + So)) * (Sno/(KNO + Sno)) * Xbh;
p(7) = kh * (Xs/Xbh)/(KX + (Xs/Xbh)) * s_no_factor;
p(8) = kh * (Xnd/Xbh)/(KX + (Xnd/Xbh)) * s_no_factor;

% sanity test
%p(isnan(p))=0;

% Finalized Reaction rates
% -----------------------------

r = [];
r(1) = 0;
r(2) = ((-1/YH) * (p(1) + p(2))) + p(7);
r(3) = frXI * (p(4) + p(5));
r(4) = (1-frXI) * (p(4) + p(5)) - p(7);
r(5) = p(1) + p(2) - p(4);
r(6) = p(3) - p(5);
r(7) = -((1-YH)/(2.86*YH))*p(2) + (1/YA)*p(3);
r(8) = -iNBM*(p(1) + p(2)) - (iNBM + 1/YA)*p(3) + p(6);
r(9) = p(8) - p(6);
r(10) = (iNBM - frXI*iNXI)*(p(4) + p(5)) - p(8);
r(11) = -((1-YH)/(YH))*p(1) - ((4.57-YA)/YA) * p(3);

end

function dfdt = react_blower_off(t,f, COD, TN, Din, Drs, Dws)

% WWTP Flow Parameters
% -----------------------
fns =2.45*10^(-3);
gamma = ((1 - fns)*(Din(t) - Dws))/(Drs + Dws);

% Previous state values
% ------------------------

Si = f(1);
Ss = f(2);
Xi = f(3);
Xs = f(4);
Xbh = f(5);
Xba = f(6);
Sno = f(7);
Snh = f(8);
Snd = f(9);
Xnd = f(10);
So = f(11);

% Reaction Rates 
% --------------------------
[r,p] = calculate_rxn(f, t);

% Incoming Influent 
% --------------------------
X_in = [];

X_in(1) = 0.05*COD(t);
X_in(2) = 0.35*COD(t);
X_in(7) = 0*TN;
X_in(8) = 0.66*TN;
X_in(9) = 0.02*TN;

X_in(3) = 0.10*COD(t);
X_in(4) = 0.35*COD(t);
X_in(5) = 0.15*COD(t);
X_in(6) = 0.00*COD(t);
X_in(10) = 0.32*TN;

X_in(11) = 0;

% ODE Equations
% --------------
D_in = Din(t);
dfdt = zeros(11,1);
dfdt(1)  = D_in * (X_in(1) - Si) + r(1);
dfdt(2)  = D_in * (X_in(2) - Ss) + r(2);
dfdt(7)  = D_in * (X_in(7) - Sno) + r(7);
dfdt(8)  = D_in * (X_in(8) - Snh) + r(8);
dfdt(9)  = D_in * (X_in(9) - Snd) + r(9);

dfdt(3)  = D_in * (X_in(3) - Xi)  + gamma*Drs*Xi + r(3);
dfdt(4)  = D_in * (X_in(4) - Xs)  + gamma*Drs*Xs + r(4);
dfdt(5)  = D_in * (X_in(5) - Xbh)  + gamma*Drs*Xbh + r(5);
dfdt(6)  = D_in * (X_in(6) - Xba)  + gamma*Drs*Xba + r(6);
dfdt(10) = D_in * (X_in(10) - Xnd) + gamma*Drs*Xnd + r(10);
dfdt(11) = D_in * (X_in(11) - Xnd) + r(11); 

% So_max = 10.0; % g/m^3l
% Kla = 4.5*24;
% Ao = Kla*(So_max - X_in(11));
% fdt(11) = Din * (X_in(11) - Xnd) + r(11) + Ao;

end

function dfdt = react_blower_on(t,f, COD, TN, Din, Drs, Dws)

% WWTP Flow Parameters
% -----------------------
fns = 2.45*10^(-3);
gamma = ((1 - fns)*(Din(t) - Dws))/(Drs + Dws);

% Previous state values
% ------------------------

Si = f(1);
Ss = f(2);
Xi = f(3);
Xs = f(4);
Xbh = f(5);
Xba = f(6);
Sno = f(7);
Snh = f(8);
Snd = f(9);
Xnd = f(10);
So = f(11);

% Reaction Rates 
% --------------------------
[r,p] = calculate_rxn(f, t);

% Incoming Influent 
% --------------------------
X_in = [];

X_in(1) = 0.05*COD(t);
X_in(2) = 0.35*COD(t);
X_in(7) = 0*TN;
X_in(8) = 0.66*TN;
X_in(9) = 0.02*TN;

X_in(3) = 0.10*COD(t);
X_in(4) = 0.35*COD(t);
X_in(5) = 0.15*COD(t);
X_in(6) = 0.00*COD(t);
X_in(10) = 0.32*TN;

X_in(11) = 0;

% ODE Equations
% --------------
Din_t = Din(t);
dfdt = zeros(11,1);
dfdt(1)  = Din_t * (X_in(1) - Si) + r(1);
dfdt(2)  = Din_t * (X_in(2) - Ss) + r(2);
dfdt(7)  = Din_t * (X_in(7) - Sno) + r(7);
dfdt(8)  = Din_t * (X_in(8) - Snh) + r(8);
dfdt(9)  = Din_t * (X_in(9) - Snd) + r(9);

dfdt(3)  = Din_t * (X_in(3) - Xi)  + gamma*Drs*Xi + r(3);
dfdt(4)  = Din_t * (X_in(4) - Xs)  + gamma*Drs*Xs + r(4);
dfdt(5)  = Din_t * (X_in(5) - Xbh)  + gamma*Drs*Xbh + r(5);
dfdt(6)  = Din_t * (X_in(6) - Xba)  + gamma*Drs*Xba + r(6);
dfdt(10) = Din_t * (X_in(10) - Xnd) + gamma*Drs*Xnd + r(10);

So_max = 10.0; % g/m^3l
Kla = 4.5*24; % hr-1
Ao = Kla*(So_max - So);
dfdt(11) = Din_t * (X_in(11) - So) + r(11) + Ao;

end
