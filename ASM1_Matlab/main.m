clc, clear, close all

load dryinfluent.mat

t_in = DRYINFLUENT(:,1);        % days
Si_in = DRYINFLUENT(:,2);       % g/m3
Ss_in = DRYINFLUENT(:,3);       % g/m3
Xi_in = DRYINFLUENT(:,4);       % g/m3
Xs_in = DRYINFLUENT(:,5);       % g/m3
Xbh_in = DRYINFLUENT(:,6);      % g/m3
Xba_in = DRYINFLUENT(:,7);      % g/m3
Xp_in = DRYINFLUENT(:,8);       % g/m3
So_in = DRYINFLUENT(:,9);       % g/m3
Sno_in = DRYINFLUENT(:,10);     % g/m3
Snh_in = DRYINFLUENT(:,11);     % g/m3
Snd_in = DRYINFLUENT(:,12);     % g/m3
Xnd_in = DRYINFLUENT(:,13);     % g/m3
Salk_in = DRYINFLUENT(:,14);    % g/m3
Q_in = DRYINFLUENT(:,15);       % m3/day


fns=2.45*10^(-3); %Non-settable fraction of the effluent suspended solids
iNBM=0.068; % g Ng?1 COD Mass of biomass per mass of COD in biomass
iNXI=0.06;  %g Ng?1 COD Mass of biomass per mass of COD in products formed by biomass decay


total_tn = Sno_in + Snh_in + Snd_in + ...
            fns*(Xnd_in + iNBM*(Xba_in + Xbh_in) + iNXI*Xi_in);
          
[data]=asm1(0.1,0,1,DRYINFLUENT, 0);
