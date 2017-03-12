function [ dd ] = asm1(dt,start_time,end_time, data_in, blower_orgs)
%ASM1 Summary of this function goes here
%   Input:
%       dt: timestep in days

time_vector = start_time:dt:end_time;
interpolated_data = interpolate_inflow(data_in);

y0 = zeros(11,1);

tspan = [0.01, 1];

opts = odeset('RelTol',1e-2,'AbsTol',1e-2);

y0(1) = get_incoming(interpolated_data, 'Si', 0, 13);
y0(2) = get_incoming(interpolated_data, 'Ss', 0, 13);
y0(3) = get_incoming(interpolated_data, 'Xi', 0, 13);
y0(4) = get_incoming(interpolated_data, 'Xs', 0, 13);
y0(5) = get_incoming(interpolated_data, 'Xbh', 0, 13);
y0(6) = get_incoming(interpolated_data, 'Xba', 0, 13);
y0(7) = get_incoming(interpolated_data, 'Sno', 0, 13);
y0(8) = get_incoming(interpolated_data, 'Snh', 0, 13);
y0(9) = get_incoming(interpolated_data, 'Snd', 0, 13);
y0(10) = get_incoming(interpolated_data, 'Xnd', 0, 13);
y0(11) = get_incoming(interpolated_data, 'So', 0, 13);


[t, f] = ode45(@(t, f) react(t,f,blower_orgs, interpolated_data, 13.9), tspan, y0, opts)

finalized_Si = f(:,1);
finalized_Ss = f(:,2);
finalized_Sno = f(:,7);
finalized_Snh = f(:,8);
finalized_Snd = f(:,9);

finalized_Xi = f(:,3);
finalized_Xs = f(:,4);
finalized_Xbh = f(:,5);
finalized_Xba = f(:,6);
finalized_Xnd = f(:,10);
finalized_So = f(:,11);

fns=2.45*10^(-3); %Non-settable fraction of the effluent suspended solids
iNBM=0.068; % g Ng?1 COD Mass of biomass per mass of COD in biomass
iNXI=0.06;  %g Ng?1 COD Mass of biomass per mass of COD in products formed by biomass decay


total_tn = finalized_Sno + finalized_Snh + finalized_Snd + ...
            fns*(finalized_Xnd + iNBM*(finalized_Xba + finalized_Xbh) + iNXI*finalized_Xi);

dd = f;

end

function interpolated_data = interpolate_inflow(data)
% perform interpolation

    fitted_data = {};
    ft = fittype('smoothingspline');
    for i=2:length(data(1,:))
        [xData, yData] = prepareCurveData( data(:,1), data(:,i));
        [fitresult, gof] = fit( xData, yData, ft );
        fitted_data{i-1} = fitresult;
    end
    interpolated_data = fitted_data;
    
end


function X_incoming = get_incoming(interpolated_data, type, t, t_end)
    
    tt = {'T','Si','Ss','Xi','Xs','Xbh','Xba','Xp','So','Sno','Snh','Snd','Xnd','Salk','Q'};

    % cycle the profile if the request time is beyond influent data
    while t > t_end
       t = t - t_end;
    end
    
    % get all influent value at any time
    values = [];
    for i=1:length(interpolated_data)
        y = interpolated_data{i};
        values = [values, y(t)];
    end
    
    % note we don't count T therefore we need to make sure
    Index = find(not(cellfun('isempty', strfind(tt, type))))-1;
    X_incoming = values(Index);

end

function y = is_blower_on(t)


start_t = [0, 0.039, 0.13, 0.217, 0.341, 0.45, 0.56, 0.64, 0.72, 0.81, 0.85, 0.92];
end_t = [0.10, 0.094, 0.19, 0.302, 0.405, 0.51, 0.58, 0.66, 0.76, 0.83, 0.89, 1];

time = [start_t', end_t'];
val = 0;
for i=1:length(time(:,1))
    if t >= time(i,1) && t <= time(i,2)
        val = 1;
    end
end
    
y = val;
end


function dfdt = react(t,f, Blower, interpolated_data, t_end)
    
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
p(1) = uH * (Ss)/(KS + Ss) * ((So)/(KOH + So)) * Xbh;
p(2) = uH * ((Ss)/(KS + Ss)) * ((KOH)/(KOH + So)) * ((Sno)/(KNO + Sno)) * nNOg * Xbh;
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

% WWTP Flow Parameters
% -----------------------

Qin = get_incoming(interpolated_data, 'Q', t, t_end);
Vat = 2050;
Qws = 75;
Qrs = 7600;
Din = Qin/Vat;
Drs = Qrs/Vat;
Dws = Qws/Vat;

gamma = ((1 - fns)*(Din - Dws))/(Drs + Dws);

% Incoming Influent 
% --------------------------

X_in = [];

X_in(1) = get_incoming(interpolated_data, 'Si', t, t_end);
X_in(2) = get_incoming(interpolated_data, 'Ss', t, t_end);
X_in(7) = get_incoming(interpolated_data, 'Sno', t, t_end);
X_in(8) = get_incoming(interpolated_data, 'Snh', t, t_end);
X_in(9) = get_incoming(interpolated_data, 'Snd', t, t_end);

X_in(3) = get_incoming(interpolated_data, 'Xi', t, t_end);
X_in(4) = get_incoming(interpolated_data, 'Xs', t, t_end);
X_in(5) = get_incoming(interpolated_data, 'Xbh', t, t_end);
X_in(6) = get_incoming(interpolated_data, 'Xba', t, t_end);
X_in(10) = get_incoming(interpolated_data, 'Xnd', t, t_end);
X_in(11) = get_incoming(interpolated_data, 'So', t, t_end);


% ODE Equations
% --------------

dfdt = zeros(11,1);
dfdt(1)  = Din * (X_in(1) - Si) + r(1);
dfdt(2)  = Din * (X_in(2) - Ss) + r(2);
dfdt(7)  = Din * (X_in(7) - Sno) + r(7);
dfdt(8)  = Din * (X_in(8) - Snh) + r(8);
dfdt(9)  = Din * (X_in(9) - Snd) + r(9);

dfdt(3)  = Din * (X_in(3) - Xi)  + gamma*Drs*Xi + r(3);
dfdt(4)  = Din * (X_in(4) - Xs)  + gamma*Drs*Xs + r(4);
dfdt(5)  = Din * (X_in(5) - Xbh)  + gamma*Drs*Xbh + r(5);
dfdt(6)  = Din * (X_in(6) - Xba)  + gamma*Drs*Xba + r(6);
dfdt(10) = Din * (X_in(10) - Xnd) + gamma*Drs*Xnd + r(10);

So_max = 10.0; % g/m^3l
Kla = 4.5*24;
Ao = Kla*(So_max - X_in(11));

% check if blower is on?
t_new = t;
while 1 < t_new
    t_new = t_new - 1;
end

if is_blower_on(t_new)
    dfdt(11) = Din * (X_in(11) - Xnd) + r(11) + Ao;
else
    dfdt(11) = Din * (X_in(11) - Xnd) + r(11);    
end



end


