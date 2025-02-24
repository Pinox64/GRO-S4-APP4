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

%% 3) Reponse à une impulsion
%--------------------------------------------------------------------------
dt = 0.001;                        % Pas de temps de l'evaluation de l'impulsion
t = 0:dt:0.1;                       % Echelle temporelle de l'evaluation de l'impulsion
tt = -0.1:dt:0.1;
Impultion_entree = (tt>= 0 & tt<=0.01).*U_arduino;  % Signal d'entrée
%Methode2
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
% Calcul conv
rep1 = conv(h_a1, Impultion_entree, 'same')*dt;   % Convolution pour obtenir la reponse a l'impulsion de 10ms @ 3v3
i1 = conv(h_e1, Impultion_entree, 'same')*dt;
rep2 = conv(h_a2, Impultion_entree, 'same')*dt;   % Convolution pour obtenir la reponse a l'impulsion de 10ms @ 3v3
i2 = conv(h_e2, Impultion_entree, 'same')*dt;
%-----affichage-----
% Find indices where t >= 0
idx = tt >= 0;
Impultion_entree_trunc = Impultion_entree(idx);
figure('Name',"Reponse à l'impulsion")
clf;
subplot(2,1,1)
%yyaxis left
plot(t, rep1);
title("Haut-Parleur #1")
ylabel("Pression (pa)")
%hold on;
%yyaxis right
%ylabel("courant(A)");
%plot(t, i1);
%grid on;

subplot(2,1,2)
title("Haut-Parleur #2")
yyaxis left
plot(t, rep2);
hold on;
yyaxis right
ylabel("courant(A)");
plot(t, i2);
grid on;
yyaxis right
%title("Réponse du système à une impulsion de 10ms");

figure('name', 'test')
plot(t,rep1);


%% 4) Reponse à un pwm
%--------------------------------------------------------------------------

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
    H_k1(ik) = evalfr(H_global_Fostex, 1i*k(ik)*w0);
    H_k2(ik) = evalfr(H_global_SEAS, 1i*k(ik)*w0);
end

% Calculer les coefficients de fourrier de la sortie
for ik = 1:length(k)
    Y_k1(ik) = X_k(ik)*H_k1(ik);
    Y_k2(ik) = X_k(ik)*H_k2(ik);
end

% Reconstruction du signal de sortie à partir des coeffs de fourrier
y_t1 = zeros(size(t));
y_t2 = zeros(size(t));
x_t = zeros(size(t));
for ik = 1: length(k)
    x_t = x_t + X_k(ik)*exp(1i*k(ik)*w0*t);
    y_t1 = y_t1 + Y_k1(ik)*exp(1i*k(ik)*w0*t);
    y_t2 = y_t2 + Y_k2(ik)*exp(1i*k(ik)*w0*t);
end

% Additionne plusieurs périodes
nPeriods = 4;           % Number of periods to repeat
T = t(end) - t(1) + (t(2)-t(1));  % Compute period length (assuming uniform sampling)

t_new = [];
x_new = [];
y1_new = [];
y2_new = [];

for k = 0:nPeriods-1
    t_new = [t_new, t + k*T];
    x_new = [x_new, x_t];
    y1_new = [y1_new, y_t1];
    y2_new = [y2_new, y_t2];
end

% Affichage de la reponse au signal PWM
figure('Name', 'Reponse au pwm');
yyaxis right
p1 = plot(t_new, y1_new, '-r');
hold on;
p2 = plot(t_new, y2_new, '-g');
ylabel("Pression (Pa)");
xlabel("Temps (s)");
yyaxis left
p3 = plot(t_new, x_new, '-b');
ylabel("Tension (V)");
legend([p1,p2,p3], "y(t) (Hp1)", "y(t) (Hp2)","x(t)");


