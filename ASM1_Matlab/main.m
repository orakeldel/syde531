clc, clear, close all

% Define Influent Flow Equation
Q_in = @(t) 1+(-0.32*cos(2*1*pi*t) - 0.18*sin(2*1*pi*t)) + ...
            (0.23*cos(2*2*pi*t) - 0.01*sin(2*2*pi*t)) + ...
           (-0.06*cos(2*3*pi*t) - 0.01*sin(2*3*pi*t));

% Define Influent COD Equation
COD = @(t) 1+(0.24*cos(2*1*pi*t) - 0.20*sin(2*1*pi*t)) + ...
            (-0.09*cos(2*2*pi*t) + 0.07*sin(2*2*pi*t)) + ...
           (0.04*cos(2*3*pi*t) - 0.02*sin(2*3*pi*t));

       


start_t = [0, 0.039, 0.13, 0.217, 0.341, 0.45, 0.56, 0.64, 0.72, 0.81, 0.85, 0.92];
end_t = [0.010, 0.094, 0.19, 0.302, 0.405, 0.51, 0.58, 0.66, 0.76, 0.83, 0.89, 1];

blower_conf = [];

for i=1:length(start_t)
    blower_conf = [blower_conf; start_t(i) end_t(i) 1];
    if i < length(start_t)
        blower_conf = [blower_conf;end_t(i) start_t(i+1) 0];
    end
end
org_blow = blower_conf;
for i=1:3
    new_b = org_blow(:,1:2)+i;
    new_s = org_blow(:,3);
    new_b = [new_b, new_s];
    blower_conf = [blower_conf; new_b];
end



[t,f, tn] = asm1(Q_in, COD, blower_conf);
figure
subplot(3,1,1)
plot(t,tn,'-k')
ylabel('\bfTN (mg/L)')
xlabel('\bfDays')
grid on
hold off

subplot(3,1,2)
fplot(Q_in, [0 blower_conf(end,2)],'-k')
hold on
fplot(COD, [0 blower_conf(end,2)],'--k')
ylim([0 2.5])
legend('Inflow (m^3/day)','COD (mg/L)')
ylabel('\bfWeighted Function (-)')
xlabel('\bfDays')
grid on

subplot(3,1,3)
ylabel('\bfS_o (mg/L)')
xlabel('\bfDays')
max_y = max(f(:,11)) - 0.10* max(f(:,11));
x = [blower_conf(1,1), blower_conf(1,1), blower_conf(1,2), blower_conf(1,2)];
y = [0 max_y max_y 0];
patch(x,y,'b')
hold on
for i=2:length(blower_conf(:,1))
    if blower_conf(i,3) == 1
        x = [blower_conf(i,1), blower_conf(i,1), blower_conf(i,2), blower_conf(i,2)];
        y = [0 max_y max_y 0];
        patch(x,y,'b')
        alpha(0.05) 
    end
end
plot(t,f(:,11),'-k')
ylim([0 max(f(:,11))])
grid on


