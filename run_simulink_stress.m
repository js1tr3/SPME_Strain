% Initialize
clc;
clear all;
% Choose Simulink Model
Model = 'simulink_SPME_stress'
%% Model Setup
% Simulation Time
simtime = 60*60;
% Initial SOC
SOC_0 = 0.1;
% Input Current in Amps
Iin = -5*5;
% Positive for Discharge
% Switch
%  1 for Constant Current
%  2 for Pulses
%  3 for CC-CV
%  4 for Manual Input
sw = 3;
% exp_sw = -1;
% Manual Current Input
for k = 1:simtime+1
    t1(k) = k-1;
    if(t1(k))<30*60
        I1(k) = 0;
    else
        I1(k) = 0;
    end
    
end
I_man =[t1',I1'];
% I_man =[t_py',u_py'];
% I_man =[time',53*C_rate'];

% CV Voltage Upper Limit
Vlim = 4.2; % Volts
Vmax = 4.2;
% CV Current Limit
Ilim = -2.5; % Amps
% CV Charging Gain
CV_gain = 50;
KI = 10;
Kaw = 1;
% Ambient Temperature (Celcius)
CC.Ta = 25.3;
% Inital Temperature  (Celcius)
CC.T0 = 25.3;
%% Run Simulink
tic

parameters


% Iin = -5*Capacity;

simOut = sim(Model);
elapsed_time = toc


%%
    
u = simOut.u;
simOut.tout = [0:1:length(u)-1]';
t = simOut.tout;

%%
Vt = simOut.Vt;
de_t_t = simOut.de_t_t;
expansion = simOut.exp;
Temp = simOut.Temp;

cs_p = simOut.cs_p;
cs_n = simOut.cs_n;
css_p= cs_p(:,end);
css_n= cs_n(:,end);
csavg_n = simOut.csavg_n;
csavg_p = simOut.csavg_p;
SOC_p = simOut.SOC_p;
SOC_n = simOut.SOC_n;
ce = simOut.ce;
ce_p = ce(:,end);
ce_n = ce(:,1);
sigma_r = simOut.sigma_r;
sigma_t = simOut.sigma_t;
sigma_h = simOut.sigma_h;
phie = simOut.phie;
eta_pl = simOut.eta_pl;
j_pl = simOut.j_pl;
j_n = simOut.j_n;
Li_loss_rate = simOut.Li_loss_rate;
Li_loss = simOut.Li_loss;
Q_loss = simOut.Q_loss;

% plot_plating
% plot_slider_pl
