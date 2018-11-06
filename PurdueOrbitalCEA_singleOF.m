% Purdue Orbital: CEA
% NEEDS CEA FILES TO RUN
% THIS CEA RUNTHROUGH IS SET UP TO RUN FOR A GIVEN O/F RATIO NOT A RANGE OF
% VALUES

CEA_RUN = true;
CEA_SAVE_FILE = 'cea_PurdueOrbital.mat';

inp = containers.Map;
inp('type') = 'eq';          % Sets the type of CEA calculation
inp('p') = 50;               % Chamber pressure
inp('p_unit') = 'atm';       % Chamber pressure units
inp('o/f') = 2           % Mixture ratio
inp('sup') = 16;             % Supersonic area ratios
%inp('pip') = 20;            % Pressure ratios
inp('fuel') = 'HTPB';       % Fuel name from thermo.inp
inp('fuel_t') = 298.15;       % Fuel inlet temperature
inp('fuel_t_unit') = 'K';    % Fuel Temp Units
inp('ox') = 'O2(L)';         % Ox name from thermo.inpj
inp('ox_t') = 90.17;         % Ox inlet temperature
inp('ox_t_unit') = 'K';
inp('file_name') = 'CEA_Defenition.inp';% Input/output file name

if CEA_RUN
    data = cea_rocket_run(inp);     % Call the CEA MATLAB code
    save(CEA_SAVE_FILE, 'data');
else
    load(CEA_SAVE_FILE);
end

data_eq = data('eq');
cf_eq=squeeze(data_eq('cf'));
eps_eq = squeeze(data_eq('ae/at'));
Isp2 = squeeze(data_eq('isp'));
t_c = squeeze(data_eq('t'));
o_f = squeeze(data_eq('o/f'));
c_f = squeeze(data_eq('cf'));
m=squeeze(data_eq('m'));
gamma=squeeze(data_eq('gammas'));
c_star = squeeze(data_eq('cstar'));

cf_eq=max(cf_eq);
[M,I] = max(Isp2);
Isp = M;
t_c = t_c (I,1);
m=m(I);
gamma=gamma(I, 1);
c_star = c_star(I, 1);
Isp1 = Isp(1)/9.81;

%% Outputs
fprintf('Isp Stage 1 (s): %f\n', Isp1);
fprintf('T_c (K): %f\n', t_c);
fprintf('Molecular Weight(g/mol): %f\n', m);
fprintf('Gamma: %f\n', gamma);
fprintf('O/F ideal: %f\n', o_f);
fprintf('C Star(m/s): %f\n', c_star);
fprintf('Cf : %f\n\n',cf_eq );
