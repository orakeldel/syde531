function Optimization
clc, close all

% Scenario related parameters:
s_count = 5;                            % Number of scenarios
s_pr = fspecial('gaussian', [1,5], 1);  % Probabilies of differeent scenarios
s_error = [-0.2, -0.1, 0, 0.1, 0.2];    % error of asm in scenario
theta = 0.5

%load('ASM1.m');
% Parameters
Nc = 27;
Nd = 30;            % Number of days to consider
t0 = 0;             % initial time (s)
tf = 3600*24*Nd;    % final time (s)
l0 = (tf-t0)/(Nd*Nc);
t_on_min = 15*60;   % Min ontime 15 minutes (2.3.2)
t_on_max = 3600*2;  % Max ontime is 2 hours (2.3.2)
t_off_min = 15*60;  % Min offtime 15 minutes (2.3.2)
t_off_max = 3600*2; % Max offtime is 2 hours (2.3.2)

TN_max = 12;      % mgL^-1 of total nitrogen allowed
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




x0 = [repmat(3600*24/Nc/3*2, 1, Nc*Nd) zeros(1,s_count*2)]; % Start with all ak in half the time open


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
    xs = x(1:end-s_count*2); % actual values
    sfsu = reshape(x(end-s_count*2+1:end),[2,5]);
    
    % Bounds for l0-ak = off time
    offconstr_min = repmat(t_off_min, size(xs)) - (repmat(l0, size(xs)) - xs);
    offconstr_max = (repmat(l0, size(xs)) - xs) - repmat(t_off_max, size(xs));
    onconstr_min = repmat(t_on_min, size(xs)) - xs;
    onconstr_max = xs - repmat(t_on_max, size(xs));
    % constraining maximal total nitrogen
    %tic
    [t, f, tn] = asm1(Q_in, COD, divide_on_off(xs, Nc));
    %toc
    lastTN = tn;
    counter = counter + 1;
    
    reptn = repmat(tn, [1, s_count]);
    reperr = repmat(s_error, [size(tn,1) 1]);
    tnerr = reptn .* (1 + reperr);
    maxtn = max(tnerr(t>1,:));
    ceq = max(zeros(size(s_error)),maxtn - repmat(TN_max, size(maxtn))) + sfsu(1,:) - sfsu(2,:);
    %nitrogenconstr = 0;
    c = [offconstr_min,offconstr_max,onconstr_min,onconstr_max];
end

nonlcon = @nonlinearconstraints;


function stop = outfun(x, optimValues, state)
    % Preparing x for plotting
    plotXandTN(x, lastTN, Nc, t0, tf);    
    save(strcat('step',int2str(optimValues.iteration),'_',int2str(Nc),'_',int2str(Nd),'_1.mat'),'x');
    stop = false;
end

function y = objective(x)
    % min 1/(tc^(Nc*Nd) - t0) sum(ak)
    %optim_coeff = 1/(Nc*Nd); %<- Do not care for this in minimization...
    xs = x(1:end-s_count*2); % actual values
    csf = 10;
    csu = 1000;
    sfsu = reshape(x(end-s_count*2+1:end),[2,5]);
    z = sfsu'*[csf;csu];
    zbar = s_pr*z;
    z2bar = s_pr*(z.^2);
    varz = z2bar-zbar^2;
    
    y = sum(xs) + zbar + theta*sqrt(varz);
end


% Set zero as lower bound to avoid getting negative sf/su

% This function does not ensure constraints in intermediat solutions, so
% has to run till the end, and that takes forever and some.
% options = optimoptions('fmincon','Display', 'iter-detailed', 'MaxIter', 10, 'Algorithm', 'sqp','OutputFcn', @outfun, 'DiffMinChange', 1);
% x = fmincon(@(x) objective(x), x0, A, b, Aeq, beq, zeros(size(x0)), [], nonlcon, options)

% This function ensures the constraints, but seams not to do anything :(
options = setoptimoptions('AlwaysHonorConstraints','All','Display', 'iter', 'MaxIter', 10);
[sol, fval, exitflag, output] = minimize(@(x) objective(x), x0, A, b, Aeq, beq, zeros(size(x0)), [], nonlcon, options)

save(strcat('final','_',int2str(Nc),'_',int2str(Nd),'.mat'),'sol','fval','exitflag','output');




xs = x(1:end-s_count*2); % actual values
[t, f, tn] = asm1(Q_in, COD, divide_on_off(xs, Nc));
plotXandTN(x, tn, Nc, t0, tf);    
save(strcat('final_',int2str(Nc),'_',int2str(Nd),'.mat'),'x');

end