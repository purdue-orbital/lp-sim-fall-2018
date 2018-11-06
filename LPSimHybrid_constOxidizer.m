% Lumped Parameter Simulation
%here keep mo_ox constant in order to keep chamber pressure and regression
%rate up
clear all;

%% Inputs:
delta_t = 0.1;
t_max = 10;                           % time in seconds when the simulation ends
t = 0:delta_t:t_max;                    
OtoF_init = 1.5;                        % initial OtoF, changes with time          
N=1;                                    % number of ports

R_init = 0.015;                         % initial port radius[m]
M_fuel = 30.870912;                     % total mass of the fuel[kg]
rho_f = 1.066*1000;                     % fuel density [kg/m^3]
Vol = M_fuel/rho_f;                     % Volume of the fuel [m^3]
m_fuel_init = 0.174996;                 % Initial fuel mass flow [kg/s]
m_ox=OtoF_init*m_fuel_init;             % Initial oxidizer mass flow [kg/s]

% CEA inputs
fuel = string('HTPB');                  % 'HTPB' for Al add , string('AL(cr)')
fuel_mix = [1];                         % 1 for Al add [0.8, 0.2]
oxidizer = 'O2';                        % ["H2O(L)" "H2O2(L)"]
ox_mix = 1;                             % [0.05 0.95]
temp = 298;                             % [298 298]
eps = 180;                              % expansion ratio (as set by design)

p_c_init = 20;                           % Initial chamber pressure [Bar]

A_t=3.724791*10^(-4);                   % Throat area in [m^2]

a = 0.117;                              % regression rate coefficient r=a*Gox^n
n = 0.956;                              % regression rate coefficient r=a*Gox^n
Units = 10;                             % '10' if GOx is in cm units, '1' of in m

% calculate the length of the grain based on the surface area that is (with the fuel density)
% needed to achieve the initial mass flow
L=m_fuel_init/(0.001*a*(m_ox/N/(pi*R_init^2*Units))^n*2*pi*R_init*N*rho_f);


%% Calculations
% Initialize arrays
R = zeros(1,length(t));
r = zeros(1,length(t)); 
m_dot_fuel = zeros(1,length(t));
m_dot_ox = zeros(1,length(t));
m_dot_tot = zeros(1,length(t));
T = zeros(1,length(t));
Isp = zeros(1,length(t));
cstar = zeros(1,length(t));
p_c = zeros(1,length(t));
G_Ox= zeros(1,length(t));
OtoF=zeros(1,length(t));
m_sum = zeros(1,length(t));

%% First time step

m_dot_fuel(1) = m_fuel_init;                        % [kg/s]
OtoF(1)=m_ox/m_dot_fuel(1);                         % [kg/s]
m_dot_tot(1) = m_ox + m_dot_fuel(1);                % total propellant flow [kg/s]
R(1) = R_init;                                      % [m]
G_Ox(1) = m_ox/N/(pi*R(1)^2)/Units;                 % Oxidizer mass flux
r(1) = a*(G_Ox(1))^n/1000;                          % regression rate [m/s]
% callCEA with momentary state in chamber
[T(1),Isp(1),cstar(1)]=callCEA(p_c_init,OtoF(1),eps,fuel, fuel_mix, oxidizer, ox_mix, temp);
p_c(1) = m_dot_tot(1)*cstar(1)/A_t*10^(-6);         % calculate combustion pressure for next step
m_sum(1) = m_dot_fuel(1)*delta_t;                   % total fuel burnt so far

for i=2:length(t)
    'Time :'
    i*delta_t
    R(i) = R(i-1)+r(i-1)*delta_t;                   % Radius [m]
    G_Ox(i) = m_ox/N/(pi*R(i)^2)/Units;             % [kg/(s*m^2)]
    r(i) = a*(G_Ox(i))^n/1000;                      % [m] 
    m_dot_fuel(i) = 2*N*pi*rho_f*R(i)*L*r(i);                             
    OtoF(i)=m_ox/m_dot_fuel(i);
    m_sum(i) = m_sum(i-1) + delta_t*m_dot_fuel(i);
    m_dot_tot(i) = m_ox + m_dot_fuel(i);            % total mass flow over time
    [T(i),Isp(i),cstar(i)]=callCEA(p_c(i-1),OtoF(i),eps,fuel, fuel_mix, oxidizer, ox_mix, temp);
    p_c(i)=m_dot_tot(i)*cstar(i)/A_t*10^(-6);
    if m_sum(i)>M_fuel
        t_burn = t(i);
        t_max=t_burn;
        break;
    end
end

Thrust = m_dot_tot.*Isp;                         % Thrust profile over time [N]

Thrust_avg = trapz(m_sum(1:i),Thrust(1:i))/m_sum(i)
ISP_avg = trapz(m_sum(1:i),Isp(1:i))/m_sum(i)
p_avg= trapz(m_sum(1:i),p_c(1:i)*10)/m_sum(i)
T_avg = sum(T)/(t_max/delta_t+1)
%% Plots

figure(1)
plot(t,m_dot_ox,'g');
hold on
plot(t,m_dot_fuel,'b');
plot(t,m_dot_tot,'r');
title('Mass flows')
xlabel('Time [s]')
ylabel('Mass flow [kg/s]')

figure(2)
plot(t,Thrust);
title('Thrust')
xlabel('Time [s]')
ylabel('Thrust [N]')

figure(3)
plot(t,Isp);
title('ISP')
xlabel('Time [s]')
ylabel('ISP [m/s]')

figure(4)
plot(t,p_c*10);
title('Chamber Pressure')
xlabel('Time [s]')
ylabel('Pressure [bar]')

figure(5)
plot(t,T);
title('Temperature')
xlabel('Time [s]')
ylabel('Temperature [K]')

figure(6)
plot(t,r*1000);
title('Regression rate')
xlabel('Time [s]')
ylabel('Regression rate [mm/s]')

figure(7)
plot(t,OtoF);
title('O/F')
xlabel('Time [s]')
ylabel('O/F')

figure(8)
plot(t(1:i),m_sum(1:i));
title('Fuel Mass Burned')
xlabel('Time [s]')
ylabel('mass [kg]')

function [T_c, ISP, cstar] = callCEA (p_c, OF, eps,fuel, fuel_mix, oxidizer, ox_mix, temp)

CEA_RUN = true;
CEA_SAVE_FILE = 'cea_PurdueLumped.mat';

inp = containers.Map;
inp('type') = 'eq';              % Sets the type of CEA calculation
inp('p') = p_c;                % Chamber pressure
inp('p_unit') = 'bar';              % Chamber pressure units
inp('o/f') = OF;              % Mixture ratio
inp('sup') = eps;               % Supersonic area ratios
inp('fuel') = fuel;% Fuel name from thermo.inp
inp('fuel_wt%') = fuel_mix;
inp('ox') = oxidizer;% Fuel name from thermo.inp
inp('ox_wt%') = ox_mix; 
inp('ox_t') = temp;               % Ox inlet temperature
inp('ox_t_unit') = 'K';
inp('file_name') = 'CEAcalls.inp';% Input/output file name

if CEA_RUN
    data = cea_rocket_run(inp);     % Call the CEA MATLAB code
    save(CEA_SAVE_FILE, 'data');
else
    load(CEA_SAVE_FILE);
end

data_eq = data('eq');
ISP=squeeze(data_eq('isp'));
ISP=ISP(2);
cstar=squeeze(data_eq('cstar'));
cstar=cstar(2);
T_c=squeeze(data_eq('t'));
T_c = T_c(1);
delete 'CEAcalls1.inp'
delete 'CEAcalls1.out'
delete 'cea_PurdueLumped.mat'
end