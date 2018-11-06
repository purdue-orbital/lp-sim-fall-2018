% Purdue Orbital: HTPB 
% NEEDS CEA FILES TO RUN
clear;
clc;

%%Constants
g0 = 9.81;                              %m/s^2    
p_a = 101325;                           %atmospheric pressure (0 if vacuum, 101325 at sea level)

%Inputs
thrust = 1000;                          % N
N=1;                                    % number of ports 
R_init = 0.03;                          % initial port radius[m]
rho_f = 900;                            % fuel density [kg/m^3]
p_c= 20;                                % chamber pressure [bar]  10bar = 1MPa

%CEA inputs
fuel =  'C32H66';                       %[ "HTPB" "AL(cr)"];  % 'HTPB';
fuel_mix =  1;                          %Fuel mixture ratio by weight (Vector if more than one component)
fuel_temp = 298;                        %Fuel temperature (K) changed if using LH2, etc
oxidizer = 'N2O';                      % Oxidizer Defenition
ox_mix = 1;                             % Oxidizer mixture ratio by weight (Vector if more than one component)
ox_temp = 298;                          % Oxidizer temperature (K) changed if using LOx, etc
epsilon = 18;                           % expansion ratio (as set by design) (low if sea level, high if vacuum)

%looping sets
delta_t = 1;                           %time step (s) lower time step, more accuracy
t_max = 20;                            % time in seconds when the simulation ends (our test will be 10 seconds)
t = 0:delta_t:t_max;                   %time vector
OtoF_init = 5.7;                       % initial OtoF, changes with time (can be changed to find optimal)           
%N2O O/F = 7.9  
%H2O2 O/F = 6.5
%H2O2 Paraffin O/F = 5.7
%GOX O/F = 2.1

%Regression Constants
a = 0.034;                            % [m^(1+2n) kg^(-n) s^(n-1)] regression rate coefficient r=a*Gox^n (empirical data from research)
n = 0.96;                             % regression rate coefficient r=a*Gox^n (empirical data from research)

%IMPORTANT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%If the output regression rate units from the empirical regression constants are in mm/s... 
%SET FLAG TO 1 
%Otherwise, flag will remain 0.

flag = 1; 

if flag == 1
    a = a/1000;
end 

%First run of CEA
[p_e,T_c, ISP, epsilon_eq, MW, o_f_1, gamma, c_star, cf_eq] = callCEAfirst (p_c, OtoF_init, fuel, fuel_mix, fuel_temp, oxidizer, ox_mix, ox_temp,epsilon);

% thrust and burn time
isp = ISP / g0;                         %convert CEA ISP from m/s to s
exit_vel = isp * g0;                    %exit velocity [m/s]

%pe = .87607e5;
pe = p_e;
A1 = (2*gamma^2/(gamma-1)) * (2/(gamma+1))^((gamma+1)/(gamma-1));
A2 = 1 - (pe(1)/(p_c*10^5))^((gamma-1)/gamma);
A3 = (pe(1)/(p_c*10^5)-p_a/(p_c*10^5))*epsilon;

CF(1) = sqrt(A1*A2) + A3; 
Isp(1) = CF(1)*c_star/g0;


A_t = thrust/(CF*(p_c*10^5));
m_dot_total = (p_c*10^5)*A_t/(0.95*c_star);

% Mass Flow                        
m_dot_ox = (m_dot_total*o_f_1) / (1+o_f_1);   %Solved using substitution
m_dot_fuel = m_dot_ox / o_f_1;

p_c_SI = p_c * 10^5;                    %pressure in Pa from bar
A_e= A_t*epsilon_eq(1,2);               %exit area [m^2] (epsilon(1,1) is 0)
d_t = sqrt(4 * A_t / pi);               %throat diameter [m]
d_e = sqrt(4 * A_e / pi);               %exit diameter [m]
A_p = pi*R_init^2;
m_flux_ox = m_dot_ox/(N*(pi*R_init^2))                         % Oxidizer mass flux [kg/m^2-s]
r_dot = a*(m_flux_ox)^n;

L = (m_dot_fuel/N)/(2*pi*R_init*rho_f*r_dot);

%CEA Loop from LPSim
m_dot_fuel(1) = m_dot_fuel                        % [kg/s]
OtoF(1)=o_f_1;                         % [kg/s]
m_dot_total(1) = m_dot_total;                 % total propellant flow [kg/s]
R(1) = R_init;                                           % [m]
m_flux_ox(1) = m_flux_ox;                         % Oxidizer mass flux
r_dot(1) = r_dot;                          % regression rate [m/s]
B1 = 2*pi*N*rho_f*L*a*(m_dot_ox/(pi*N))^n;
B2 = a*(2*n+1)*(m_dot_ox/(pi*N))^n;
B3 = R_init^(2*n+1);
x = (1-2*n)/(1+2*n);
m_dot_fuel(1) = B1*(B2*t(1)+B3)^x;
thrust(1) = (m_dot_ox(1)+m_dot_fuel(1))*Isp(1)*g0

% callCEA with momentary state in chamber
[p_e(1), T(1),isp(1),cstar(1)]=callCEA(p_c,OtoF(1),epsilon,fuel, fuel_mix, oxidizer, ox_mix, ox_temp);
p_c(1) = p_c;         % calculate combustion pressure for next step
m_sum(1) = m_dot_fuel(1)*delta_t;     % total fuel burnt so far


for i=2:length(t)
    Time = i*delta_t
    R(i) = R(i-1)+r_dot(i-1)*delta_t;                   % Radius [m]
    m_flux_ox(i) = m_dot_ox/(N*(pi*R(i)^2));             % [kg/(s*m^2)]
    r_dot(i) = a*(m_flux_ox(i))^n; 
    
    % [m]
    B1 = 2*pi*N*rho_f*L*a*(m_dot_ox/(pi*N))^n;
    B2 = a*(2*n+1)*(m_dot_ox/(pi*N))^n;
    B3 = R_init^(2*n+1);
    x = (1-2*n)/(1+2*n);
    m_dot_fuel(i) = B1*(B2*Time+B3)^x;
    OtoF(i)=m_dot_ox/m_dot_fuel(i);
    m_sum(i) = m_sum(i-1) + delta_t*m_dot_fuel(i);
    m_dot_total(i) = m_dot_ox + m_dot_fuel(i);            % total mass flow over time
    
    [p_e(i), T(i),isp(i),cstar(i), MW(i), gamma(i)]=callCEA(p_c(i-1),OtoF(i),epsilon,fuel, fuel_mix, oxidizer, ox_mix, ox_temp);
    p_c(i)= (m_dot_total(i)*cstar(i)/(A_t))*10^(-5);
    
    A1 = (2*gamma(i)^2/(gamma(i)-1)) * (2/(gamma(i)+1))^((gamma(i)+1)/(gamma(i)-1));
    A2 = 1 - (p_e(i)/(p_c(i)*10^5))^((gamma(i)-1)/gamma(i));
    A3 = (p_e(i)/(p_c(i)*10^5)-p_a/(p_c(i)*10^5))*epsilon;
    CF(i) = sqrt(A1*A2) + A3;  
    Isp(i) = CF(i)*cstar(i)/g0;
    thrust(i) = (m_dot_ox+m_dot_fuel(i))*Isp(i)*g0;
end

%% Plots

figure(1)
%plot(t,m_dot_ox);
plot(t,m_dot_fuel);
title('Fuel Mass flows')
xlabel('Time [s]')
ylabel('Mass flow [kg/s]')

figure(2)
plot(t,m_dot_total);
title('Total Mass flows')
xlabel('Time [s]')
ylabel('Mass flow [kg/s]')

figure(3)
plot(t,thrust);
title('Thrust')
xlabel('Time [s]')
ylabel('Thrust [N]')

figure(4)
plot(t,Isp);
title('ISP')
xlabel('Time [s]')
ylabel('ISP [s]')
axis([0 t_max 200 400])

figure(5)
plot(t,p_c);
title('Pressure')
xlabel('Time [s]')
ylabel('Pressure [bar]')

figure(6)
plot(t,T);
title('Temperature')
xlabel('Time [s]')
ylabel('Temperature [K]')

r_dotmm = r_dot * 1000;  %Converts regression rate to mm/s
figure(7)
plot(t,r_dotmm);
title('Regression rate')
xlabel('Time [s]')
ylabel('Regression rate [mm/s]')

figure(8)
plot(t,OtoF);
title('O/F')
xlabel('Time [s]')
ylabel('O/F')

figure(9)
plot(t(1:i),m_sum(1:i));
title('Fuel Mass Burned')
xlabel('Time [s]')
ylabel('mass [kg]')


%%CEA Calls
function [p_e,t_c, ISP, epsilon_eq, MW, o_f, gamma, c_star, cf_eq] = callCEAfirst (p_c, OF, fuel, fuel_mix, fuel_temp, oxidizer, ox_mix, ox_temp,epsilon)
    CEA_RUN = true;                             %initializes program
    CEA_SAVE_FILE = 'cea_PurdueOrbital.mat';    %save file

    inp = containers.Map;                       %format
    inp('type') = 'eq';                         % Sets the type of CEA calculation
    inp('p') = p_c;                             % Chamber pressure
    inp('p_unit') = 'bar';                      % Chamber pressure units
    inp('o/f') = OF;                            % Mixture ratio
    inp('sup') = epsilon;                       % Supersonic area ratios (use sup because we aren't perfectly expanded)
    %inp('pip') = 20;                           % Pressure ratios
    inp('fuel') = fuel;                         % Fuel name from thermo.inp
    inp('fuel_wt%') = fuel_mix;
    inp('fuel_t') = fuel_temp;
    inp('ox') = oxidizer;                       % Fuel name from thermo.inp
    inp('ox_wt%') = ox_mix; 
    inp('ox_t') = ox_temp;                          % Ox inlet temperature
    inp('ox_t_unit') = 'K';
    inp('file_name') = 'test.inp';   % Input/output file name

    if CEA_RUN
        data = cea_rocket_run(inp);     % Call the CEA MATLAB code
        save(CEA_SAVE_FILE, 'data');
    else
        load(CEA_SAVE_FILE);
    end

    data_eq = data('eq');
    cf_eq=squeeze(data_eq('cf'));       %thrust coefficient
    epsilon_eq = squeeze(data_eq('ae/at'));     %expansion ratio
    t_c = squeeze(data_eq('t'));        %combustion temperature
    o_f = squeeze(data_eq('o/f'));      %oxidizer to fuel ratio
    MW=squeeze(data_eq('m'));            %molecular weight
    gamma=squeeze(data_eq('gammas'));   %specific heat ratio
    c_star = squeeze(data_eq('cstar')); %characteristic velocity
    ISP=squeeze(data_eq('isp')); 
    p = squeeze(data_eq('p'));
    ISP=ISP(2);                         %Isp

    
    cf_eq=max(cf_eq);                   
    [M,I] = max(ISP);
    p_e = p(2);
    t_c = t_c (I,1);
    MW=MW(I);
    gamma=gamma(I, 1);
    c_star = c_star(I, 1);
end

function [p_e, T_c, ISP, cstar, MW, gamma] = callCEA (p_c, OF, eps,fuel, fuel_mix, oxidizer, ox_mix, ox_temp)

    CEA_RUN = true;
    CEA_SAVE_FILE = 'cea_orbital.mat';

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
    inp('ox_t') = ox_temp;               % Ox inlet temperature
    inp('ox_t_unit') = 'K';
    inp('file_name') = 'test.inp';% Input/output file name

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
    cf_eq=squeeze(data_eq('cf'));
    esp_eq = squeeze(data_eq('ae/at'));
    o_f = squeeze(data_eq('o/f'));
    c_f = squeeze(data_eq('cf'));
    MW=squeeze(data_eq('m'));
    gamma=squeeze(data_eq('gammas'));
    c_star = squeeze(data_eq('cstar'));
    p = squeeze(data_eq('p'));

    cf_eq=max(cf_eq);
    [M,I] = max(ISP);
    MW=MW(I);
    gamma=gamma(I, 1);
    p_e = p(2);
    c_star = c_star(I, 1);
end