clear; close all; clc;

% grog1401
% boum0704
% gagr1618


%% 1) PARAMETRES DU PROBLEME
% -------------------------------------------------------------------------
% Commun
Rs    = 70;          % [ohms] Resistance serie
d     = 1;           % [m] distance
c     = 340;         % [m/s] vitesse du son
rho   = 1;           % [kg/m^3] densite air approx

% Haut-parleur #1 (Fostex)
HP1.Re = 6.4;        % [ohms]   Resistance bobine
HP1.Le = 0.051e-3;   % [H]      Inductance bobine
HP1.Bl = 10.8;       % [Tm]     Facteur de couplage
HP1.Mm = 13.3e-3;    % [kg]     Masse mobile
HP1.Rm = 0.50;       % [kg/s]   Amortissement
HP1.Km = 935;        % [N/m]    Raideur
HP1.Sm = 0.0201;     % [m^2]    Surface membrane

% Haut-parleur #2 (SEAS)
HP2.Re = 3.0;        
HP2.Le = 0.05e-3;    
HP2.Bl = 4.2;        
HP2.Mm = 10.5e-3;
HP2.Rm = 0.48;        
HP2.Km = 400;        
HP2.Sm = 0.0222;     

% Pour la tension
U_arduino = 3.3;     % [V]  amplitude max
Imax_lim  = 0.05;    % [A]  (50 mA) limite courant

%% 2) Fonction de transferts
% -------------------------------------------------------------------------
% Electrique + Mecanique
H_elec_mec_Fostex = tf([0, HP1.Mm, HP1.Rm, HP1.Km], [(HP1.Mm*HP1.Le), (HP1.Mm*HP1.Re + HP1.Mm*Rs + HP1.Rm*HP1.Le), (HP1.Rm*HP1.Re + HP1.Rm*Rs + HP1.Km*HP1.Le + HP1.Bl*HP1.Bl), (HP1.Km*HP1.Re + HP1.Km*Rs)]);
H_elec_mec_SEAS   = tf([0, HP2.Mm, HP2.Rm, HP2.Km], [(HP2.Mm*HP2.Le), (HP2.Mm*HP2.Re + HP2.Mm*Rs + HP2.Rm*HP2.Le), (HP2.Rm*HP2.Re + HP2.Rm*Rs + HP2.Km*HP2.Le + HP2.Bl*HP2.Bl), (HP2.Km*HP2.Re + HP2.Km*Rs)]);


% Mecanique
H_mec_Fostex =  tf([0,0,HP1.Bl],[HP1.Mm, HP1.Rm, HP1.Km]);
H_mec_SEAS   =  tf([0,0,HP2.Bl],[HP2.Mm, HP2.Rm, HP2.Km]);

% Accoustique
H_acc_Fostex = tf([(rho*HP1.Sm) / (2*pi()*d), 0, 0],[1], "InputDelay", d/c);
H_acc_SEAS   = tf([(rho*HP2.Sm) / (2*pi()*d), 0, 0],[1], "InputDelay", d/c);

% Globale
H_global_Fostex = H_elec_mec_Fostex * H_acc_Fostex;
H_global_SEAS   = H_elec_mec_SEAS * H_acc_SEAS;

dt = 0.01;
t = 0:dt:10;
ha_SEAS = impulse(minreal(H_global_SEAS),t);
ha_Fostrex = impulse(minreal(H_global_Fostrex),t);
figure('Name',"Réponse à l'impulsion")
subplot(2,1,1);
plot(h_SEAS,t);
plot(h_Fostrex);


