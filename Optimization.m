function Optimization
clc, close all
%load('ASM1.m');
% Parameters
Nc = 15;
Nd = 2;            % Number of days to consider
t0 = 0;             % initial time (s)
tf = 3600*24*Nd;    % final time (s)
l0 = (tf-t0)/(Nd*Nc);
t_on_min = 15*60;   % Min ontime 15 minutes (2.3.2)
t_on_max = 3600*2;  % Max ontime is 2 hours (2.3.2)
t_off_min = 15*60;  % Min offtime 15 minutes (2.3.2)
t_off_max = 3600*2; % Max offtime is 2 hours (2.3.2)

TN_max = 10;      % mgL^-1 of total nitrogen allowed
% Define Influent Flow Equation
Q_in = @(t) 1+(-0.32*cos(2*1*pi*t) - 0.18*sin(2*1*pi*t)) + ...
            (0.23*cos(2*2*pi*t) - 0.01*sin(2*2*pi*t)) + ...
           (-0.06*cos(2*3*pi*t) - 0.01*sin(2*3*pi*t));

% Define Influent COD Equation
COD = @(t) 1+(0.24*cos(2*1*pi*t) - 0.20*sin(2*1*pi*t)) + ...
            (-0.09*cos(2*2*pi*t) + 0.07*sin(2*2*pi*t)) + ...
           (0.04*cos(2*3*pi*t) - 0.02*sin(2*3*pi*t));


% To optimize
% ak = On time in cycle k (k=1..N)


% min 1/(tc^(Nc*Nd) - t0) sum(ak)
%optim_coeff = 1/(Nc*Nd); %<- Do not care for this in minimization...

fun = @(x) sum(sum(x));

x0 = repmat(3600*24/Nc/2, 1, Nc*Nd); % Start with all ak in half the time open


% No linear constraints
A = [];
b = [];
Aeq = [];
beq = [];

% Lower and upper bound of ak
lb = repmat(t_on_min, 1, Nc);
ub = repmat(t_on_max, 1, Nc);

%function y = FTNforIntegral(t)
    % integral from t0 to tf = [max(0, FTN(x) - TN_max)]^2dt
    % FTN(x(t)) = SNoat(t) + SNHat(t) + SNDat(t) + fns* (XNDat(t) + iNBM*(XBHat(t) + XBAat(t)) + iNXI*XIat(t))
%    t = mod(t, 3600*24); % We only have info for one day...
%    ftn = SNoat(t) + SNHat(t) + SNDat(t) + fns* (XNDat(t) + iNBM*(XBHat(t) + XBAat(t)) + iNXI*XIat(t));
%    ftn = ftn - repmat(TN_max,size(ftn));
%    y = max(ftn .^ 2 , 0);
% end

%ftnintfun = @FTNforIntegral;
lastTN = [];
counter = 0;
function [c,ceq] = nonlinearconstraints(x)
    
    % Bounds for l0-ak = off time
    offconstr_min = repmat(t_off_min, size(x)) - (repmat(l0, size(x)) - x);
    offconstr_max = (repmat(l0, size(x)) - x) - repmat(t_off_max, size(x));
    % constraining maximal total nitrogen
    tic
    [t, f, tn] = asm1(Q_in, COD, divide_on_off(x, Nc));
    toc
    lastTN = tn;
    counter = counter + 1
    nitrogenconstr = max(0,(max(tn(t>1)) - TN_max));
    %nitrogenconstr = 0;
    c = [offconstr_min,offconstr_max];
    ceq = [nitrogenconstr];
end

nonlcon = @nonlinearconstraints;


function stop = outfun(x, optimValues, state)
    % Preparing x for plotting
    plotXandTN(x, lastTN, Nc, t0, tf);    
    save(strcat('step',int2str(optimValues.iteration),'_',Nc,'_',Nd,'.mat'),'x');
end

options = optimoptions('fmincon','Display', 'iter', 'MaxIter', 10, 'Algorithm', 'sqp','OutputFcn', @outfun);

x = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options)

xx = 1:size(x,2);

subplot(2,1,1);
plot(lastT,lastTN);
subplot(2,1,2);
plot(xx,x);

end