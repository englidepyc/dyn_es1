clear, clc, close all
% 55:44 fine di ascoltare;
% linea 33 non è vettore colonna (?) imoporta poco

% Dichiarazione costanti:
a1 = 25;   %manovella
a3 = 120;  %glifo
a4 = 10;   %bielletta
a6 = 80;   %distanza  O1-O4 (telaio)
a7 = 40;   %distanza O1-guida corsoio (telaio)

%caricamento posizioni ed estrarre variabili colonna
load Posizioni_Itipo.mat
th1d = P(:,1); %deg
th3d = P(:,2); %deg
th4d = P(:,3); %deg
a2 = P(:,4); %mm
a5 = P(:,5); %mm

%caricamento velocità ed estrarre variabili colonna
th1p = 1; %rad/s
load Velocità_Itipo.mat
th3p = V(:,2);
th4p = V(:,3);
a2p = V(:,4);
a5p = V(:,5);
th2p = th3p;

%forza esterna nota (10N - sempre resistente)
F_e5 = -sign(a5p)*10;
RR = zeros(length(F_e5),5); % matrice con tutte le forze (sol) per angolo di manovella

for i=1:length(F_e5);
%definizione F_e
F_e=[F_e5(i); 0; 0; 0; 0]; %vettore colonnna (no?)

%definizione di C_F
C_F = zeros(5,5);
C_F(1,1)=1;

%definizione di C_R
C_R=[-1, 0, 0, 0, 0;
    0, -1, 1, 0, 0;
    1, 0, 0, 1, 1;
    0, 1, 0, tand(th4d(i)), tand(th3d(i)+90);
    0, 0, 0, a3*(sind(th3d(i))-tand(th4d(i))*cosd(th3d(i))), a2(i)*(sind(th3d(i))-tand(th3d(i)+90)*cosd(th3d(i)));];

%calcolo di R
R = -inv(C_R)*C_F*F_e;
RR(i,:) = R';

end


F_53_x = RR(:,1);
F_53_y = RR(:,2);
F_T5 = RR(:,3);
F_43_x = RR(:,4);
F_43_y = F_43_x.*tand(th4d);
F_23_x = RR(:,5);
F_23_y = F_23_x.*tand(th3d+90);


F_53 = sqrt(F_53_y.^2 + F_53_x.^2);
th_F53 = atan2d(F_53_y,F_53_x);

F_43 = sqrt(F_43_y.^2 + F_43_x.^2);
th_F43 = atan2d(F_43_y,F_43_x);

F_23 = sqrt(F_23_y.^2 + F_23_x.^2);
th_F23 = atan2d(F_23_y,F_23_x);


F_21_x = -F_23_x;
F_21_y = -F_23_y;

M_e1 = a1*F_21_x.*sind(th1d) - a1*F_21_y.*cosd(th1d);


figure(30)
plot(th1d,M_e1),grid on
xlim([0 360])
xlabel('\theta_{1} [deg]')
ylabel('M_{e_1 } [Nmm]')

figure(40)
plot(th1d,M_e1),grid on
xlim([0 360])
xlabel('\theta_{1} [deg]')
ylabel('M_{e_1 } [Nmm]')