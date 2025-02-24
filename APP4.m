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
% Electrique + Mecanique
%H_elec_mec_Fostex = tf([0, HP1.Mm, HP1.Rm, HP1.Km], [(HP1.Mm*HP1.Le), (HP1.Mm*HP1.Re + HP1.Mm*Rs + HP1.Rm*HP1.Le), (HP1.Rm*HP1.Re + HP1.Rm*Rs + HP1.Km*HP1.Le + HP1.Bl*HP1.Bl), (HP1.Km*HP1.Re + HP1.Km*Rs)]);
%H_elec_mec_SEAS   = tf([0, HP2.Mm, HP2.Rm, HP2.Km], [(HP2.Mm*HP2.Le), (HP2.Mm*HP2.Re + HP2.Mm*Rs + HP2.Rm*HP2.Le), (HP2.Rm*HP2.Re + HP2.Rm*Rs + HP2.Km*HP2.Le + HP2.Bl*HP2.Bl), (HP2.Km*HP2.Re + HP2.Km*Rs)]);

s = tf('s');
% Mecanique
H_mec_Fostex =  (HP1.Bl)/(HP1.Mm*s^2 + HP1.Rm*s + HP1.Km)
H_mec_SEAS   =  (HP2.Bl)/(HP2.Mm*s^2 + HP2.Rm*s + HP2.Km)

% Electrique
H_elec_Fostex  = 1/((HP1.Le + HP1.Bl * H_mec_Fostex)*s + HP1.Re + Rs);
H_elec_SEAS    = 1/((HP2.Le + HP2.Bl * H_mec_Fostex)*s + HP2.Re + Rs);

% Accoustique
H_acc_Fostex = tf([(rho*HP1.Sm) / (2*pi()*d), 0, 0],[1], "InputDelay", d/c)
H_acc_SEAS   = tf([(rho*HP2.Sm) / (2*pi()*d), 0, 0],[1], "InputDelay", d/c)

% Globale
H_global_Fostex = H_elec_Fostex * H_acc_Fostex * H_mec_Fostex
H_global_SEAS   = H_elec_SEAS * H_acc_SEAS * H_mec_SEAS

dt = 0.01;
t = 0:dt:10;
%ha_SEAS = impulse(minreal(H_global_SEAS),t);
%ha_Fostrex = impulse(minreal(H_global_Fostrex),t);
figure('Name',"Diagrammes de bode")
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

%% Reponse à une impulsion
dt = 0.0001;                        % Pas de temps de l'evaluation de l'impulsion
t = 0:dt:0.1;                       % Echelle temporelle de l'evaluation de l'impulsion
tt = -0.1:dt:0.1;
Impultion_entree = (tt>= 0 & tt<=0.01).*U_arduino;  % Signal d'entrée
RI = impulse(H_global_Fostex,t);    % Réponse à l'impulsion de dirac
rep = conv(RI, Impultion_entree, 'same')*dt;   % Convolution pour obtenir la reponse a l'impulsion de 10ms @ 3v3

%-----affichage-----
% Find indices where t >= 0
idx = tt >= 0;
Impultion_entree_trunc = Impultion_entree(idx);
figure('Name',"Reponse à l'impulsion")
yyaxis left
plot(t, rep);
hold on;
grid on;
yyaxis right
plot(t, Impultion_entree_trunc);
title("Réponse du système à une impulsion de 10ms");

%% Reponse à un pwm

%Definition du signal
w0 = 2*pi*PWM_Freq;
t = linspace(0, PWM_Period, 100); %temps du signal pwm (une période)
y_pwm = (t < PWM_Period/2).*U_arduino;

%Calculer la série de fourrier du signal d'entrée
N_harmoniques = 50;                     % Nombre d'harmoniques à évaluer
X_k = zeros(2*N_harmoniques+1);         % Coefficients de fourrier
k = [-N_harmoniques:1:N_harmoniques];   % Indice des harmoniques

for ik = 1:length(X_k)
    X_k(ik) = 1/PWM_Period*trapz(t, y_pwm.*exp(-1i * k(ik)* w0 * t));
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
    H_k(ik) = evalfr(H_global_Fostex, j*k(ik)*w0);
end

