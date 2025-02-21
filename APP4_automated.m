%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCRIPT EXEMPLE MATLAB - APP4 (GRO410)
%
% 1) Paramètres de deux hauts-parleurs (HP1, HP2)
% 2) Construction des fonctions de transfert:
%    - He(s) = I(s)/U(s)
%    - Ha(s) = P(s)/U(s)
% 3) Analyse fréquentielle (Bode)
% 4) Réponse à PWM (par superposition d'harmoniques)
% 5) Réponse impulsionnelle (convolution)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

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
U_arduino = 3.3;     % [V] amplitude max
Imax_lim  = 0.05;    % [A]  (50 mA) limite courant

%% 2) CONSTRUCTION FONCTIONS DE TRANSFERT
% -------------------------------------------------------------------------
% On va creer des fonctions pour generer He(s) et Ha(s)
% Input: structure HP => output: objets tf en Matlab

% On definit d'abord une fonction inline - tu peux la mettre en fin de script
makeTF = @(HP) buildTF(HP, Rs, rho, d, c);

[He1, Ha1] = makeTF(HP1);
[He2, Ha2] = makeTF(HP2);

disp('--- Haut-parleur #1 (Fostex) ---');
He1
Ha1

disp('--- Haut-parleur #2 (SEAS) ---');
He2
Ha2

%% 3) ANALYSE FREQUENTIELLE (BODE)
% -------------------------------------------------------------------------
% Bode du courant: He(s),  et Bode de la pression: Ha(s)
% Compare HP1 vs HP2
figure('Name','Bode He(s)');
bode(He1, 'b', He2, 'r--',{1,1e5}); grid on;
legend('HP1','HP2');
title('Fonction de transfert I(s)/U(s) - Bode');

figure('Name','Bode Ha(s)');
bode(Ha1, 'b', Ha2, 'r--',{1,1e5}); grid on;
legend('HP1','HP2');
title('Fonction de transfert P(s)/U(s) - Bode');

% Observations:
% - On verra la/les frequences de resonance
% - On peut comparer amplitude maxi du courant, etc.

%% 4) REPONSE A UN SIGNAL PWM (frequence 100 Hz, 50% duty)
% -------------------------------------------------------------------------
% Approche: decomposition en serie de Fourier -> on superpose n harmoniques
f0    = 100;             % [Hz]
omega0= 2*pi*f0;         % [rad/s]
nharm = 15;              % nombre d'harmoniques a considerer
Udc   = U_arduino/2;     % composante DC = ~1.65 V

% Coeffs amplitude pour l'onde carree 50%  (ex: a_k, b_k)
% Formule canonique: 
%   amplitude(k) = 4*U_arduino/(pi*(2k-1)) ?? 
%   => attention: l'onde carree +/- 1 => faire un offset 
% Ici, on peut generer numeriquement ou calculer analytique (ex. methodes).
% Pour la demonstration, je fais un petit calcul generique:
ampls = zeros(1,nharm);
phases= zeros(1,nharm);

% L'onde carree 50%:  pour k impairs, amplitude = 4*Vmax/(pi*k) sin ...
% Mais note: la moyenne (DC) = 1.65. 
% => On se concentre sur l'AC.  Variation: +/-1.65 autour de la moyenne.
for k=1:nharm
    % si k pair => amplitude = 0
    % si k impair => amplitude = 4*(U_arduino/2)/[pi*(2k-1)] => etc.
    if mod(k,2)==1
       % odd
       ampls(k) = 4*(U_arduino/2)/(pi*k);  % simple approximation
       phases(k)= 0;  % onde carree type cos...
    else
       ampls(k) = 0;
       phases(k)= 0;
    end
end

% => Dans la realite, verifie la correspondance a la bonne definition
%    offset DC = Udc

% On evalue la reponse harmonique pour HP1
Npts  = 1000;     % nb points temporels
T     = 1/f0;     % periode
tPWM  = linspace(0, 3*T, Npts);  % on trace sur 3 periodes, p.ex.
pHP1  = zeros(1,Npts);
iHP1  = zeros(1,Npts);

for tIdx=1:Npts
    tt = tPWM(tIdx);
    % la tension en tout point: u(t) = Udc + somme_{k=1->nharm} ...
    % On calcule la contribution en sortie p(t), i(t) en superposant.
    valP = 0;
    valI = 0;
    for k=1:nharm
        if ampls(k)==0, continue; end
        % amplitude = ampls(k)
        % freq = k*omega0
        % phase= phases(k)
        % => reponse p(t) = |Ha1(jk w0)| * ampls(k) cos(...) 
        [magHa,phaseHa] = bode(Ha1, k*omega0);
        [magHe,phaseHe] = bode(He1, k*omega0);
        magHa   = magHa(:);  % convert en array
        phaseHa = phaseHa(:);
        magHe   = magHe(:);
        phaseHe = phaseHe(:);

        % Contribution en pression
        amplitudeP = magHa*ampls(k);
        phiP       = phases(k) + phaseHa;
        valP       = valP + amplitudeP*cos(k*omega0*tt + phiP);

        % Contribution en courant
        amplitudeI = magHe*ampls(k);
        phiI       = phases(k) + phaseHe;
        valI       = valI + amplitudeI*cos(k*omega0*tt + phiI);
    end

    % + composante DC => dans l'excitation, la DC ne cree pas de courant AC 
    %   => si on veut inclure la DC, sur un HP reel, la DC agit sur x(t), 
    %   possible, mais on simplifie ici
    % A titre d'exemple, on ignore la contribution DC ds la pression 
    % (un offset DC ne cree pas de son). 
    % Idem pour le courant => un certain offset sur la bobine ?

    pHP1(tIdx) = valP; 
    iHP1(tIdx) = valI;
end

figure('Name','PWM HP1');
subplot(2,1,1)
plot(tPWM*1e3, pHP1, 'b','LineWidth',1.2);
xlabel('t [ms]'); ylabel('Pression [Pa]');
title('HP1 - reponse PWM (approx. harmoniques)');

subplot(2,1,2)
plot(tPWM*1e3, iHP1*1e3, 'r','LineWidth',1.2);
xlabel('t [ms]'); ylabel('Courant [mA]');
title('Courant HP1');

% On peut checker la valeur max
maxIHP1 = max(abs(iHP1));
disp(['Courant max HP1 ~ ', num2str(maxIHP1*1e3,'%.2f'), ' mA']);

%% 5) REPONSE IMPULSIONNELLE (convolution)
% -------------------------------------------------------------------------
% Approche 1: on recupere la reponse impulsionnelle h_a(t) = invLap{Ha(s)}
% Approche 2: on recupere la decomposition en fractions partielles
% Approche 3: on fait une simulation simulink
%
% Ici: illustration "pseudo-numerique" -> on echantillonne Ha(s) 
% via l'impulse() ou step() + numerics ?

% (A) Calcul numerique direct de la reponse impulsionnelle 
%    impulse(Ha1) donne p(t) suite a un Dirac en entree => mathematiquement
%    Ca sera dur a interpreter direct => on va echantillonner
tSim = linspace(0,0.1,2000);  % 0.1 s
[yImp1, tImp1] = impulse(Ha1, tSim);

% yImp1(t) = h_a(t), la reponse impulsionnelle
% L'impulsion qu'on applique reellement n'est pas un Dirac(t) mais 3.3V 
% pendant 10ms -> on fait la convolution discrete
dt = tImp1(2)-tImp1(1);

% Construire le signal d'entree
uImp = zeros(size(tImp1));
idxEnd = floor(0.010/dt);        % sur 10 ms
uImp(1:idxEnd) = U_arduino;      % 3.3 V

% Convolution discrete
pPercHP1 = conv(uImp, yImp1)*dt;  
% Le vecteur resultant aura une longueur = length(uImp)+length(yImp1)-1
tPerc = (0:length(pPercHP1)-1)*dt;

figure('Name','Impulsion HP1 (Convolution)');
plot(tPerc*1e3, pPercHP1, 'b','LineWidth',1.2);
xlabel('t [ms]');
ylabel('Pression [Pa]');
title('HP1 - Reponse a impulsion 3.3V/10ms (convolution)');

% Meme demarche pour le courant ?
% On recupere impulse(He1) -> i(t) suite a Dirac(t) = h_e(t),
% on convolution => i(t) pour l'impulsion 10ms
[yIE1, tIE1] = impulse(He1, tSim);
dtI = tIE1(2)-tIE1(1);
uImpI = zeros(size(tIE1));
idxEndI = floor(0.010/dtI);
uImpI(1:idxEndI) = U_arduino;

iPercHP1 = conv(uImpI, yIE1)*dtI;  
tIPerc   = (0:length(iPercHP1)-1)*dtI;

figure('Name','Impulsion HP1 - courant');
plot(tIPerc*1e3, iPercHP1*1e3,'r','LineWidth',1.2);
xlabel('t [ms]');
ylabel('Courant [mA]');
title('HP1 - Courant impulsion 3.3V/10ms (convolution)');

[maxIp, idxMaxIp] = max(abs(iPercHP1));
disp(['HP1 - Courant impulsion max ~ ', num2str(maxIp*1e3,'%.2f'), ' mA, vers t= ',num2str(tIPerc(idxMaxIp)*1e3),' ms']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FONCTION LOCALE: buildTF
% Construit He(s) et Ha(s) pour un HP donne, + distance d, c, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [He,Ha] = buildTF(HP, Rs, rho, dist, c)
    % HP contient: Re, Le, Bl, Mm, Rm, Km, Sm
    % Rs = resistance serie
    % dist, c, rho => acoustique

    s = tf('s');  % variable Laplace symbolique

    % Z elec(s) = Re + Rs + Le s
    Zelec = HP.Re + Rs + HP.Le*s;

    % Admittance mec = X(s)/I(s) => Hm(s)
    %   (Mm s^2 + Rm s + Km) X = Bl I
    %   => X/I = (Bl)/(Mm s^2 + Rm s + Km)
    Hm = (HP.Bl)/(HP.Mm*s^2 + HP.Rm*s + HP.Km);

    % => He(s) = I/U = 1 / [Zelec + Bl s * Hm]
    %    Sous forme => 1 / [Zelec + Bl s * (Bl/(Mm s^2 + Rm s + Km))]
    He = 1/(Zelec + HP.Bl*s * Hm);

    % Sous-systeme acoustique:
    %   p(t) = (rho Sm / 2 pi dist) * d^2 x/dt^2   *  delta(t - dist/c)
    % => P(s)/X(s) = (rho Sm / 2 pi dist) * s^2 * e^{- dist/c s}
    scaleAc = (rho * HP.Sm)/(2*pi*dist);
    Hacou = tf( [scaleAc 0 0], 1, 'InputDelay', dist/c );  % s^2 + delai

    % => Ha(s) = Hacou(s)* Hm(s)* He(s) = P(s)/U(s)
    Ha = series(series(Hacou, Hm), He); 
end
