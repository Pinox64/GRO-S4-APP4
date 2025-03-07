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
U_arduino = 3.3;            % [V]  amplitude max
Imax_lim  = 0.05;           % [A]  (50 mA) limite courant

PWM_Freq = 100;             % [Hz] frequence du pwm
PWM_Period = 1/PWM_Freq;    % [S]  Période du pwm

%% 2) Fonction de transferts
% -------------------------------------------------------------------------

% Electrique + Mecanique (enlevé)
%H_elec_mec_Fostex = tf([0, HP1.Mm, HP1.Rm, HP1.Km], [(HP1.Mm*HP1.Le), (HP1.Mm*HP1.Re + HP1.Mm*Rs + HP1.Rm*HP1.Le), (HP1.Rm*HP1.Re + HP1.Rm*Rs + HP1.Km*HP1.Le + HP1.Bl*HP1.Bl), (HP1.Km*HP1.Re + HP1.Km*Rs)]);
%H_elec_mec_SEAS   = tf([0, HP2.Mm, HP2.Rm, HP2.Km], [(HP2.Mm*HP2.Le), (HP2.Mm*HP2.Re + HP2.Mm*Rs + HP2.Rm*HP2.Le), (HP2.Rm*HP2.Re + HP2.Rm*Rs + HP2.Km*HP2.Le + HP2.Bl*HP2.Bl), (HP2.Km*HP2.Re + HP2.Km*Rs)]);

s = tf('s');
% Mecanique
H_mec_Fostex =  minreal((HP1.Bl)/(HP1.Mm*s^2 + HP1.Rm*s + HP1.Km))
H_mec_SEAS   =  minreal((HP2.Bl)/(HP2.Mm*s^2 + HP2.Rm*s + HP2.Km))

% Electrique
H_elec_Fostex  = minreal(1/((HP1.Le + HP1.Bl * H_mec_Fostex)*s + HP1.Re + Rs))
H_elec_SEAS    = minreal(1/((HP2.Le + HP2.Bl * H_mec_SEAS)*s + HP2.Re + Rs))

[elec_num_coefs_Fostex, elec_den_coefs_Fostex] = tfdata(H_elec_Fostex, 'v');
[elec_num_coefs_SEAS, elec_den_coefs_SEAS] = tfdata(H_elec_SEAS, 'v');

% Accoustique
H_acc_Fostex = minreal(tf([(rho*HP1.Sm) / (2*pi*d), 0, 0],1, "InputDelay", d/c))
H_acc_SEAS   = minreal(tf([(rho*HP2.Sm) / (2*pi*d), 0, 0],1, "InputDelay", d/c))

% Global
H_global_Fostex = minreal(series(H_acc_Fostex, series(H_mec_Fostex,H_elec_Fostex)))
H_global_SEAS   = minreal(series(H_acc_SEAS, series(H_mec_SEAS, H_elec_SEAS)))

[global_num_coefs_Fostex, global_den_coefs_Fostex] = tfdata(H_global_Fostex, 'v');
[global_num_coefs_SEAS, global_den_coefs_SEAS] = tfdata(H_global_SEAS, 'v');

%% 3) Bode
%--------------------------------------------------------------------------
dt = 0.01;
t = 0:dt:10;
%ha_SEAS = impulse(minreal(H_global_SEAS),t);
%ha_Fostrex = impulse(minreal(H_global_Fostrex),t);

figure('Name',"Diagrammes de bode", 'Color', 'w')
subplot(2,2,1);
bode(H_elec_Fostex, 'b', H_elec_SEAS, '-r',{1,1e5});
grid on;
title("Elec");
subplot(2,2,2);
bode(H_mec_Fostex, 'b', H_mec_SEAS, '-r',{1,1e5});
grid on;
title("Mec");
subplot(2,2,3);
bode(H_global_Fostex, 'b', H_global_SEAS, '-r',{1,1e5});
grid on;
title("Global");
subplot(2,2,4);
bode(H_acc_Fostex, 'b', H_acc_SEAS, '-r',{1,1e5});
grid on;
title("Accoustique");
sgtitle('Diagrammes de bode', 'FontSize', 16, 'FontWeight', 'bold');

%% 3) Reponse à une impulsion
%--------------------------------------------------------------------------
Fe = 100e4;             % Frequence d'echantillonage (augmenté pour obtenir une approximation plus exacte se rapprochant plus du simulink)
dt = 1/Fe;              % Periode echantillonage
t = 0:dt:0.1;           % Echelle temporelle de l'evaluation de l'impulsion
tt = -0.1:dt:0.1;
Impultion_entree = (tt >= 0 & tt <= 0.01).*U_arduino;  % Signal d'entrée
idx = tt >= 0;
Impultion_entree_trunc = Impultion_entree(idx);

% Methode 1
[R1, P1, K1] = residue(H_global_Fostex.Numerator{1}, H_global_Fostex.Denominator{1});
[Re1, Pe1, Ke1] = residue(H_elec_Fostex.Numerator{1}, H_elec_Fostex.Denominator{1});
[R2, P2, K2] = residue(H_global_SEAS.Numerator{1}, H_global_SEAS.Denominator{1});
[Re2, Pe2, Ke2] = residue(H_elec_SEAS.Numerator{1}, H_elec_SEAS.Denominator{1});
h_a1=0;
h_e1=0;
h_a2=0;
h_e2=0;

for ii=1:length(P1)
    h_a1 = h_a1+R1(ii)*exp(P1(ii)*t);
end
for ii=1:length(Pe1)
    h_e1 = h_e1+Re1(ii)*exp(Pe1(ii)*t);
end
for ii=1:length(P2)
    h_a2 = h_a2+R2(ii)*exp(P2(ii)*t);
end
for ii=1:length(P2)
    h_e2 = h_e2+Re2(ii)*exp(Pe2(ii)*t);
end
%h_a1=impulse(H_global_Fostex,t);
%h_e1=impulse(H_elec_Fostex,t);
%h_a2=impulse(H_global_SEAS,t);
%h_e2=impulse(H_elec_SEAS,t);
% Calcul conv
yGlobalFostex = conv(h_a1, Impultion_entree, 'same')*dt;   % Convolution pour obtenir la reponse a l'impulsion de 10ms @ 3v3
yElecFostex = conv(h_e1, Impultion_entree, 'same')*dt;     % **% Division par 2 pour compenser l'intégration sur un vecteur symétrique doublant l'aire effective de l'impulsion.
yGlobalSEAS = conv(h_a2, Impultion_entree, 'same')*dt;   % Convolution pour obtenir la reponse a l'impulsion de 10ms @ 3v3
yElecSEAS = conv(h_e2, Impultion_entree, 'same')*dt;

%Remettre le delai qui est supprimé par la methode des residue
delay = d/c;  % delay in seconds
yGlobalFostex = interp1(t, yGlobalFostex, t - delay, 'linear', 0);
yGlobalSEAS = interp1(t, yGlobalSEAS, t - delay, 'linear', 0);

% Methode 2 (pour valider)
%yGlobalFostex = lsim(H_global_Fostex, Impultion_entree_trunc, t);
%yElecFostex = lsim(H_elec_Fostex, Impultion_entree_trunc, t);
%yGlobalSEAS = lsim(H_global_SEAS, Impultion_entree_trunc, t);
%yElecSEAS = lsim(H_elec_SEAS, Impultion_entree_trunc, t);

%-----affichage-----
% Création de la figure avec un fond blanc
figure('Name', "Réponse impulsionnelle", 'Color', 'w');

% Réponse impultionelle globale Fostex
subplot(2,2,1);
plot(t, yGlobalFostex, 'b-', 'LineWidth', 2);
grid on;
title('Réponse Globale Fostex', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Temps [s]', 'FontSize', 12);
ylabel('Pression [Pa]', 'FontSize', 12);
xlim([0 max(t)]);
ylim([min(yGlobalFostex)-0.1, max(yGlobalFostex)+0.1]);

% Réponse impultionelle en courant Fostex
subplot(2,2,2);
plot(t, yElecFostex, 'b-', 'LineWidth', 2);
grid on;
title('Réponse Courant Fostex', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Temps [s]', 'FontSize', 12);
ylabel('Courant [A]', 'FontSize', 12);
xlim([0 max(t)]);
ylim([min(yElecFostex)-0.1, max(yElecFostex)+0.1]);

% Réponse impultionelle globale SEAS
subplot(2,2,3);
plot(t, yGlobalSEAS, 'b-', 'LineWidth', 2);
grid on;
title('Réponse Globale SEAS', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Temps [s]', 'FontSize', 12);
ylabel('Pression [Pa]', 'FontSize', 12);
xlim([0 max(t)]);
ylim([min(yGlobalSEAS)-0.1, max(yGlobalSEAS)+0.1]);

% Réponse impultionelle en courant SEAS
subplot(2,2,4);
plot(t, yElecSEAS, 'b-', 'LineWidth', 2);
grid on;
title('Réponse Courant SEAS', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Temps [s]', 'FontSize', 12);
ylabel('Courant [A]', 'FontSize', 12);
xlim([0 max(t)]);
ylim([min(yElecSEAS)-0.1, max(yElecSEAS)+0.1]);

% Titre global de la figure
sgtitle('Simulation de la réponse impulsionnelle', 'FontSize', 16, 'FontWeight', 'bold');

%% 4) Reponse à un pwm
%--------------------------------------------------------------------------

% Définition du signal PWM
w0 = 2 * pi * PWM_Freq;                % Pulsation du signal PWM
t_pwm = linspace(0, PWM_Period, 100);   % Temps sur une période du PWM
y_pwm = (t_pwm < PWM_Period/2) * U_arduino;  % Signal carré (50% duty cycle)

%Calculer la série de fourrier du signal d'entrée
N_harmoniques = 10;                     % Nombre d'harmoniques à évaluer
X_k = zeros(2*N_harmoniques+1);         % Coefficients de fourrier
k = [-N_harmoniques:1:N_harmoniques];   % Indice des harmoniques

for ik = 1:length(X_k)
    X_k(ik) = 1/PWM_Period*trapz(t_pwm, y_pwm.*exp(-1i * k(ik)* w0 * t_pwm));
end

%Afficher les séries de fourrier du signal pwm
figure('name','fourrier');
subplot(2,1,1);
stem(k, abs(X_k));
xlabel("harmonique");
ylabel("abs(X_k)");
title("amplitude");
subplot(2,1,2);
stem(k, angle(X_k));
xlabel("harmonique");
ylabel("angle(X_k)");
title("Phase");

%Calculer la fct de transfert pour chaque harmonique
for ik = 1:length(k)
    H_k1(ik) = evalfr(H_global_Fostex, 1i*k(ik)*w0);
    H_k2(ik) = evalfr(H_global_SEAS, 1i*k(ik)*w0);
    He_k1(ik) = evalfr(H_elec_Fostex, 1i*k(ik)*w0);
    He_k2(ik) = evalfr(H_elec_SEAS, 1i*k(ik)*w0);
end

% Calculer les coefficients de fourrier de la sortie
for ik = 1:length(k)
    Y_k1(ik) = X_k(ik)*H_k1(ik);
    Y_k2(ik) = X_k(ik)*H_k2(ik);
    Ye_k1(ik) = X_k(ik)*He_k1(ik);
    Ye_k2(ik) = X_k(ik)*He_k2(ik);
end

% Reconstruction du signal de sortie à partir des coeffs de fourrier
yGlobalPWM_Fostex = zeros(size(t_pwm));
yGlobalPWM_SEAS = zeros(size(t_pwm));
yElecPWM_Fostex = zeros(size(t_pwm));
yElecPWM_SEAS = zeros(size(t_pwm));
x_t = zeros(size(t_pwm));
for ik = 1: length(k)
    x_t = x_t + X_k(ik)*exp(1i*k(ik)*w0*t_pwm);
    yGlobalPWM_Fostex = yGlobalPWM_Fostex + Y_k1(ik)*exp(1i*k(ik)*w0*t_pwm);
    yGlobalPWM_SEAS = yGlobalPWM_SEAS + Y_k2(ik)*exp(1i*k(ik)*w0*t_pwm);
    yElecPWM_Fostex = yElecPWM_Fostex + Ye_k1(ik)*exp(1i*k(ik)*w0*t_pwm);
    yElecPWM_SEAS =   yElecPWM_SEAS + Ye_k2(ik)*exp(1i*k(ik)*w0*t_pwm);
end

% Methode 2
% Simulation de la réponse pour Fostex et SEAS (réponses globale et en courant)
%yGlobalPWM_Fostex = lsim(H_global_Fostex, y_pwm, t_pwm);
%yElecPWM_Fostex   = lsim(H_elec_Fostex,   y_pwm, t_pwm);
%yGlobalPWM_SEAS   = lsim(H_global_SEAS,   y_pwm, t_pwm);
%yElecPWM_SEAS     = lsim(H_elec_SEAS,     y_pwm, t_pwm);

%-----affichage-----
% Création de la figure avec un fond blanc
figure('Name', "Réponse PWM", 'Color', 'w');
clf;

% Réponse globale Fostex
subplot(2,2,1);
yyaxis("left");
plot(t_pwm, yGlobalPWM_Fostex, 'b-', 'LineWidth', 2);
grid on;
title('Global Fostex', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Temps [s]', 'FontSize', 12);
ylabel('Amplitude', 'FontSize', 12);
xlim([0, max(t_pwm)]);
hold on
yyaxis("right");
plot(t_pwm, y_pwm, 'k');
ylabel("Tension d'entrée (V)")
ax = gca;
ax.YColor = 'k';  % Left axis color to black

% Réponse en courant Fostex
subplot(2,2,2);
yyaxis("left");
plot(t_pwm, yElecPWM_Fostex, 'b-', 'LineWidth', 2);
grid on;
title('Courant Fostex', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Temps [s]', 'FontSize', 12);
ylabel('Courant [A]', 'FontSize', 12);
xlim([0, max(t_pwm)]);
hold on

yyaxis("right");
plot(t_pwm, y_pwm, 'k');
ylabel("Tension d'entrée (V)")
ax = gca;
ax.YColor = 'k';  % Left axis color to black

% Réponse globale SEAS
subplot(2,2,3);
yyaxis("left");
plot(t_pwm, yGlobalPWM_SEAS, 'b-', 'LineWidth', 2);
grid on;
title('Global SEAS', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Temps [s]', 'FontSize', 12);
ylabel('Amplitude', 'FontSize', 12);
xlim([0, max(t_pwm)]);
hold on;
yyaxis("right");
plot(t_pwm, y_pwm, 'k');
ylabel("Tension d'entrée (V)")
ax = gca;
ax.YColor = 'k';  % Left axis color to black

% Réponse en courant SEAS
subplot(2,2,4);
yyaxis("left")
plot(t_pwm, yElecPWM_SEAS, 'b-', 'LineWidth', 2);
grid on;
title('Courant SEAS', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Temps [s]', 'FontSize', 12);
ylabel('Courant [A]', 'FontSize', 12);
xlim([0, max(t_pwm)]);
hold on
yyaxis("right");
plot(t_pwm, y_pwm, 'k');
ylabel("Tension d'entrée (V)")
ax = gca;
ax.YColor = 'k';  % Left axis color to black
% Titre global de la figure
sgtitle('Réponse à un signal PWM', 'FontSize', 16, 'FontWeight', 'bold');
