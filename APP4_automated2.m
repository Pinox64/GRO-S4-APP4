clear; close all; clc;

%% ------------------------------------------------------------------------
% 1) PARAMÈTRES DU PROBLÈME
% -------------------------------------------------------------------------
Rs    = 70;    % [ohms] Résistance série
d     = 1;     % [m]    Distance
c     = 340;   % [m/s]  Vitesse du son
rho   = 1;     % [kg/m^3] Densité de l'air approx.

% Haut-parleur #1 (Fostex)
HP1.Re = 6.4;         % [ohms]   Résistance bobine
HP1.Le = 0.051e-3;    % [H]      Inductance bobine
HP1.Bl = 10.8;        % [Tm]     Facteur de couplage
HP1.Mm = 13.3e-3;     % [kg]     Masse mobile
HP1.Rm = 0.50;        % [kg/s]   Amortissement
HP1.Km = 935;         % [N/m]    Raideur
HP1.Sm = 0.0201;      % [m^2]    Surface membrane

% Haut-parleur #2 (SEAS)
HP2.Re = 3.0;         
HP2.Le = 0.05e-3;     
HP2.Bl = 4.2;         
HP2.Mm = 10.5e-3;     
HP2.Rm = 0.48;        
HP2.Km = 400;         
HP2.Sm = 0.0222;      

% Pour la tension
U_arduino = 3.3;      % [V]  amplitude max
Imax_lim  = 0.05;     % [A]  (50 mA) limite courant

PWM_Freq   = 100;                 % [Hz] fréquence du PWM
PWM_Period = 1 / PWM_Freq;        % [s]  Période du PWM

%% ------------------------------------------------------------------------
% 2) FONCTIONS DE TRANSFERT
% -------------------------------------------------------------------------
s = tf('s');

% a) Sous-système mécanique (déplacement X(s) / courant I(s))
H_mec_Fostex = minreal( HP1.Bl / (HP1.Mm*s^2 + HP1.Rm*s + HP1.Km) );
H_mec_SEAS   = minreal( HP2.Bl / (HP2.Mm*s^2 + HP2.Rm*s + HP2.Km) );

% b) Sous-système électrique (I(s) / U(s))
%    En tenant compte de l’impédance bobine + force contre-électromotrice + Rs
H_elec_Fostex = minreal( 1 / ( (HP1.Le + HP1.Bl*H_mec_Fostex)*s + (HP1.Re + Rs) ) );
H_elec_SEAS   = minreal( 1 / ( (HP2.Le + HP2.Bl*H_mec_SEAS)*s + (HP2.Re + Rs) ) );

% c) Sous-système acoustique (P(s) / X(s)) avec délai d/c
H_acc_Fostex = minreal( tf( [(rho*HP1.Sm)/(2*pi*d), 0, 0], 1, "InputDelay", d/c ) );
H_acc_SEAS   = minreal( tf( [(rho*HP2.Sm)/(2*pi*d), 0, 0], 1, "InputDelay", d/c ) );

% d) Fonction de transfert globale HP1 et HP2 (P(s) / U(s))
H_global_Fostex = minreal( series( H_acc_Fostex, series(H_mec_Fostex, H_elec_Fostex) ) );
H_global_SEAS   = minreal( series( H_acc_SEAS,   series(H_mec_SEAS,   H_elec_SEAS)   ) );

% Récupération des polynômes pour usage éventuel
[global_num_Fostex, global_den_Fostex] = tfdata(H_global_Fostex, 'v');
[global_num_SEAS,   global_den_SEAS  ] = tfdata(H_global_SEAS,   'v');

[e_num_Fostex, e_den_Fostex] = tfdata(H_elec_Fostex, 'v');  % i(t) / u(t) Fostex
[e_num_SEAS,   e_den_SEAS  ] = tfdata(H_elec_SEAS,   'v');  % i(t) / u(t) SEAS

% Pour vérifier rapidement (diagrammes de Bode)
figure('Name',"Diagrammes de Bode");
subplot(2,2,1); bode(H_elec_Fostex,'b', H_elec_SEAS,'r',{1,1e5}); grid on; title("Électrique");
subplot(2,2,2); bode(H_mec_Fostex, 'b', H_mec_SEAS, 'r',{1,1e5}); grid on; title("Mécanique");
subplot(2,2,3); bode(H_global_Fostex,'b', H_global_SEAS,'r',{1,1e5}); grid on; title("Global");
subplot(2,2,4); bode(H_acc_Fostex, 'b', H_acc_SEAS, 'r',{1,1e5}); grid on; title("Acoustique");

%% ------------------------------------------------------------------------
% 3) RÉPONSE À UNE IMPULSION (10 ms à 3,3 V) PAR LA MÉTHODE DES RÉSIDUS
% -------------------------------------------------------------------------
% On veut la réponse dans le domaine du temps pour:
%   - la pression p(t) = H_global(s) * U(s)
%   - le courant i(t)  = H_elec(s)   * U(s)
%
% Ici, on n'utilise plus impulse() ni impulse_rational(), mais on fait le
% développement en fractions partielles ("residue") puis on convertit
% par exponentielles, et enfin on fait la convolution avec le signal d’entrée
% qui n’est pas un Dirac, mais bien un échelon de 10 ms à amplitude 3,3 V.

dt = 0.0001;              % pas de temps (plus fin si on veut éviter un aliasing)
t  = 0 : dt : 0.03;       % on regarde 30 ms, par exemple

% Signal d’entrée "impulsion" de 10 ms (3.3 V)
% On l'assimile à un échelon qui retombe à 0 après 10 ms
T_imp = 0.01;             % durée = 10 ms
u_in = (t>=0 & t<=T_imp) * U_arduino; 

% -----------------
% a) Résidus pour H_global_Fostex
[RF1, PF1, KF1] = residue(global_num_Fostex, global_den_Fostex);
% Reconstruit la réponse impulsionnelle h_a1(t) correspondante au Dirac
h_a1 = zeros(size(t));
for ii = 1:length(PF1)
    % Somme R(i)*exp(P(i)*t)
    h_a1 = h_a1 + RF1(ii)*exp(PF1(ii)*t);
end
% Convolution h(t)*u_in(t) => p(t)
p_Fostex = conv(h_a1, u_in, 'same') * dt;  % pression Fostex

% b) Résidus pour H_global_SEAS
[RS1, PS1, KS1] = residue(global_num_SEAS, global_den_SEAS);
h_a2 = zeros(size(t));
for ii = 1:length(PS1)
    h_a2 = h_a2 + RS1(ii)*exp(PS1(ii)*t);
end
p_SEAS = conv(h_a2, u_in, 'same') * dt;    % pression SEAS

% c) Résidus pour H_elec_Fostex (courant i(t) / u(t))
[RF_e, PF_e, KF_e] = residue(e_num_Fostex, e_den_Fostex);
h_ie1 = zeros(size(t));
for ii = 1:length(PF_e)
    h_ie1 = h_ie1 + RF_e(ii)*exp(PF_e(ii)*t);
end
i_Fostex = conv(h_ie1, u_in, 'same') * dt; % courant Fostex

% d) Résidus pour H_elec_SEAS
[RS_e, PS_e, KS_e] = residue(e_num_SEAS, e_den_SEAS);
h_ie2 = zeros(size(t));
for ii = 1:length(PS_e)
    h_ie2 = h_ie2 + RS_e(ii)*exp(PS_e(ii)*t);
end
i_SEAS = conv(h_ie2, u_in, 'same') * dt;   % courant SEAS

% -----------------------------------------------------------
% Affichage des résultats (pression et courant vs. temps)
figure('Name', "Réponse impulsive (10ms) : pression & courant");
subplot(2,2,1);
plot(t*1000, p_Fostex, 'b', 'LineWidth',1.5); grid on;
xlabel('t (ms)'); ylabel('Pression (Pa)');
title('HP1 (Fostex) - Pression');

subplot(2,2,2);
plot(t*1000, p_SEAS, 'r', 'LineWidth',1.5); grid on;
xlabel('t (ms)'); ylabel('Pression (Pa)');
title('HP2 (SEAS) - Pression');

subplot(2,2,3);
plot(t*1000, i_Fostex*1000, 'b', 'LineWidth',1.5); grid on;
xlabel('t (ms)'); ylabel('Courant (mA)');
title('HP1 (Fostex) - Courant');

subplot(2,2,4);
plot(t*1000, i_SEAS*1000, 'r', 'LineWidth',1.5); grid on;
xlabel('t (ms)'); ylabel('Courant (mA)');
title('HP2 (SEAS) - Courant');

%% ------------------------------------------------------------------------
% 4) RÉPONSE À UN SIGNAL PWM (CALCUL PAR SÉRIE DE FOURIER)
% -------------------------------------------------------------------------
% On conserve l’approche initiale : décomposition du PWM en harmoniques
% puis multiplication par la fonction de transfert en régime sinusoïdal.

w0 = 2*pi*PWM_Freq;
t_pwm = linspace(0, PWM_Period, 200); % on échantillonne une période
y_pwm = (t_pwm < PWM_Period/2) .* U_arduino; % signal d’entrée sur [0, T]

% a) Calcul des coefficients X_k
N_harmoniques = 50;
k_vect = -N_harmoniques : 1 : N_harmoniques;
X_k = zeros(size(k_vect));

for ik = 1:length(k_vect)
    k_val = k_vect(ik);
    % Intégrale numérique : coefficient de Fourier
    X_k(ik) = (1/PWM_Period) * trapz(t_pwm, y_pwm .* exp(-1i * k_val * w0 * t_pwm));
end

% Affichage du spectre d’entrée
figure('name','Série de Fourier du signal PWM');
subplot(2,1,1);
stem(k_vect, abs(X_k), 'b'); grid on;
xlabel("Indice d’harmonique k"); ylabel("|X_k|"); title("Amplitude");
subplot(2,1,2);
stem(k_vect, angle(X_k), 'b'); grid on;
xlabel("Indice d’harmonique k"); ylabel("Phase(X_k)"); title("Phase");

% b) Fonction de transfert H_global pour chaque harmonique => HP1 et HP2
H_k1 = zeros(size(k_vect));
H_k2 = zeros(size(k_vect));
for ik = 1:length(k_vect)
    s_val = 1i * k_vect(ik) * w0;  % j*k*w0
    H_k1(ik) = evalfr(H_global_Fostex, s_val);
    H_k2(ik) = evalfr(H_global_SEAS,   s_val);
end

% c) Coefficients de sortie Y_k = X_k * H_k
Y_k1 = X_k .* H_k1;
Y_k2 = X_k .* H_k2;

% d) Reconstruction du signal de sortie dans une période
t_rec = linspace(0, PWM_Period, 500); % plus de points pour la reconstruction
y_in_recons = zeros(size(t_rec));
y_out1      = zeros(size(t_rec));
y_out2      = zeros(size(t_rec));

for ik = 1:length(k_vect)
    k_val = k_vect(ik);
    y_in_recons = y_in_recons + X_k(ik)*exp(1i*k_val*w0*t_rec);
    y_out1      = y_out1      + Y_k1(ik)*exp(1i*k_val*w0*t_rec);
    y_out2      = y_out2      + Y_k2(ik)*exp(1i*k_val*w0*t_rec);
end

% e) On répète sur quelques périodes pour visualiser
nPeriods = 3;
t_tot = [];
in_tot = [];
out1_tot = [];
out2_tot = [];

Tstep = t_rec(2)-t_rec(1);  % pas de la grille ci-dessus
for np = 0 : (nPeriods-1)
    t_shift = t_rec + np*PWM_Period;
    t_tot   = [t_tot,   t_shift];
    in_tot  = [in_tot,  y_in_recons];
    out1_tot= [out1_tot,y_out1];
    out2_tot= [out2_tot,y_out2];
end

% Affichage
figure('Name','Réponse au signal PWM');
plot(t_tot, real(in_tot), 'b', 'LineWidth',1.2); hold on; grid on;
plot(t_tot, real(out1_tot), 'r', 'LineWidth',1.2);
plot(t_tot, real(out2_tot), 'g', 'LineWidth',1.2);
legend("Entrée (PWM)", "Sortie HP1 (Fostex)", "Sortie HP2 (SEAS)");
xlabel('Temps (s)'); ylabel('Pression (Pa env.)');
title('Réponse au PWM');
