clear, clc, close all
% Guida di Fairbairn modificata di primo tipo
% Analisi cinematica a campo intero

% Dichiarazione costanti:
a1 = 25;   %manovella
a3 = 120;  %glifo
a4 = 10;   %bielletta
a6 = 80;   %distanza  O1-O4 (telaio)
a7 = 40;   %distanza O1-guida corsoio (telaio)

% Angolo di manovella tra un'osservazione e l'altra all'interno del campo 0:360 deg:
delta = 1;

% Dichiarazione variabili indipendenti: posizione, velocità e accelerazione
% della manovella (membro motore)
th1 = deg2rad(0:delta:360);   %theta1 -> posizione angolare in rad
th1p = 1;                     %theta1-punto -> velocità angolare in rad/s
th1pp = 0;                    %theta1-puntopunto -> acceleraz angolare nulla

% Dichiarazione variabili dipendenti di PRIMO TENTATIVO:
th3 = deg2rad(80);   %posizione angolare glifo in rad
th4 = deg2rad(15);   %posizione angolare bielletta in rad
a2 = 25;             %posizione corsoio su glifo
a5 = 20;             %posizione corsoio finale rispetto all'asse y

% METODO DI NEWTON-RAPHSON + velocità + accelerazione
% Alla fine si avranno delle matrici con n righe (una per ogni
% osservazione) e tante colonne quanti sono le variabili dipendenti più il
% parametro libero th1
% è FONDAMENTALE l'ordine delle variabili dipendenti, di f e di J
P = zeros(length(th1),5);   %Posizioni
V = P;            %velocità
A = P;            %accelerazioni

for i=1:length(th1)
    
    % Calcolo iterativo per le posizioni
    while 1
        % Definizione equazioni di chiusura
        f = [a1*cos(th1(i))+a2*cos(th3)-a5;
             a1*sin(th1(i))+a2*sin(th3)-a7;
             a1*cos(th1(i))+a2*cos(th3)-a3*cos(th3)-a4*cos(th4);
             a1*sin(th1(i))+a2*sin(th3)-a3*sin(th3)-a4*sin(th4)+a6];

        % Definizione matrice Jacobiana
        J = [-a2*sin(th3), 0, cos(th3), -1;
             a2*cos(th3), 0, sin(th3), 0;
             -(a2-a3)*sin(th3), a4*sin(th4), cos(th3), 0;
             (a2-a3)*cos(th3), -a4*cos(th4), sin(th3), 0];

        % Calcolo iterativo
        % Calcolo delle correzioni
        dth = -inv(J)*f;
        % Aggiornamento delle variabili dipendenti
        th3 = th3 + dth(1);
        th4 = th4 + dth(2);
        a2 = a2 + dth(3);
        a5 = a5 + dth(4);
        % Calcolo errore massimo per la chiusura
        E = sum(abs(f));
        if E < 0.01 %somma delle 4 funzioni inferiore a un centesimo
           break
        end 
    end
    % Memorizzazione posizioni (angoli converititi in gradi)
    th = [rad2deg(th3), rad2deg(th4), a2, a5];
    P(i,:) = [rad2deg(th1(i)), th];
    
    % Calcolo delle velocità:
    % J da aggiornare:
    J = [-a2*sin(th3), 0, cos(th3), -1;
          a2*cos(th3), 0, sin(th3), 0;
         -(a2-a3)*sin(th3), a4*sin(th4), cos(th3), 0;
         (a2-a3)*cos(th3), -a4*cos(th4), sin(th3), 0];
    % Calcolo di B (vettore colonna):
    B = a1*[-sin(th1(i)); cos(th1(i)); -sin(th1(i)); cos(th1(i))];
    % Calcolo del vettore delle velocità dipendenti:
    thp = -inv(J)*B*th1p; % vettore colonna
    th3p = thp(1);
    th4p = thp(2);
    a2p = thp(3);
    a5p = thp(4);    
    V(i,:) = [rad2deg(th1(i)), thp'];
    
    % Calcolo delle accelerazioni
    % Calcolo della derivata di J
    Jp = [-a2p*sin(th3)-a2*th3p*cos(th3), 0, -th3p*sin(th3), 0;
          a2p*cos(th3)-a2*th3p*sin(th3), 0, th3p*cos(th3), 0;
          -a2p*sin(th3)-(a2-a3)*th3p*cos(th3), a4*th4p*cos(th4), -th3p*sin(th3), 0;
          a2p*cos(th3)-(a2-a3)*th3p*sin(th3), a4*th4p*sin(th4), th3p*cos(th3), 0];
    Bp = [-a1*th1p*cos(th1(i)); -a1*th1p*sin(th1(i)); -a1*th1p*cos(th1(i)); -a1*th1p*sin(th1(i))];
    thpp = -inv(J)*(B*th1pp + Bp*th1p + Jp*thp); % vettore colonna
    A(i,:) = [rad2deg(th1(i)), thpp'];
end