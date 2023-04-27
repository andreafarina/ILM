% STUDIO DELLA POLARIZZAZIONE CON LE MATRICI DI JONES
% convenzione cos(kz-wt + phi) e guardare l'onda che ci arriva
% agli occhi.
% Le matrici di Jones sono delle matrici 2x2 che permettono di trasformare
% le componenti complesse Ex, Ey dei campi. Sono comode per gestire la luce
% polarizzata.
% 
% Author:       Andrea Farina
% Institution:  CNR - IFN
% email:        andrea.farina@polimi.it 
% March 2021; Last revision: 17-March-2022

close all;clearvars;
%% DEFINIZIONE DELLE MATRICI DI JONES
% Polarizzatore con asse orizzontale
P_H = [1,0;0 0];

% Polarizzatore con asse verticale
P_V = [0 0; 0 1];
%% rotazione di una matrice J
% matrice rotazione rispetto all'orizzontale
Jt = @(theta)[cos(theta) sin(theta);-sin(theta) cos(theta)];

% matrice lamina ritardatrice con asse fast parallelo all'asse X 
% (quindi ritardo la componenti y)
LR = @(phi)[1 0;0 exp(1i*phi)];

% ritardatore circolare Faraday
CR = @(delta)[cos(delta) sin(delta);-sin(delta) cos(delta)];

%% COMPONENTI OTTICI
% polarizzatore con asse ottico orientato di un angolo alfa rispetto a X
P= @(alfa)Jt(-alfa)*P_H*Jt(alfa);

% lamina a l/4 con asse ottico orientato di un angolo alfa rispetto a X
L4 = @(alfa)Jt(-alfa)*LR(pi/2)*Jt(alfa);

% lamina a l/2 con asse ottico orientato di un angolo alfa rispetto a X
L2 = @(alfa)Jt(-alfa)*LR(pi)*Jt(alfa);

% lamina con sfasamento generico eta con asse ottico orientato di un angolo
% alfa rispetto a X
Leta = @(eta,alfa)Jt(-alfa)*LR(eta)*Jt(alfa);

%% PROPAGAZIONE
% Le matrici di Jones sono inserite in ordine inverso alla propagazione
%% Inserire qui il vettore [Ex;Ey] complesso del campo di partenza
% Es 1
E = [sqrt(41/2); sqrt(41/2)*exp(1i*1.3495)];
% Es 2 
E=[1;sqrt(2)*exp(-1i*pi/4)];
%E = [1;exp(1i*pi/2)];
%E = [1;0.5*exp(1i*3*pi/4)];%*exp(1i*pi/2)];
%E = [5/sqrt(2);5/sqrt(2)*exp(1i*0.2838)];
%E = [sqrt(97/2);sqrt(97/2)*exp(1i*0.8364)];
%E=E./norm(E);

%% Inserire qui la matrice di propagazione
% comporre in ordine inverso le funzioni implicite sopra nella sezione
% "Componenti Ottici"
%M=P_H;L2(pi/6);
M = 1;%L4(pi/2);%L2(pi/3);%L2(0);%L4(0);%1;%L4(pi/3);%
%% ========================================================================
% campo in uscita
Ef=M*E;
%% PLOT
N = 100;
key = 6.3/N;
zed = 1:N;
prop=exp(1i*key*zed);
figure(1),
% plot campo in ingresso (osservare la rotazione!)
subplot(1,2,1),axis square,title('Campo in ingresso'),
xlabel('E_x');ylabel('E_y'),
xlim([-norm(E) norm(E)]);ylim([-norm(E) norm(E)]);
hold;drawnow;
for k=1:numel(prop)
    plot(real(E(1)*prop(k)),real(E(2)*prop(k)),'g*'),grid;
    drawnow;%pause(0.001);
end
yline(0,'LineWidth',1.5),xline(0,'LineWidth',1.5),grid
subplot(1,2,2),axis square,title('Campo in uscita'),
xlabel('E_x');ylabel('E_y'),
xlim([-norm(Ef) norm(Ef)]);ylim([-norm(Ef) norm(Ef)]);
hold;drawnow;
for k=1:numel(prop)
    plot(real(Ef(1)*prop(k)),real(Ef(2)*prop(k)),'g*'),grid;
    drawnow;%pause(0.001);
end
yline(0,'LineWidth',1.5),xline(0,'LineWidth',1.5),grid

%% PARAMETRI DI STOKES
S(1)=(Ef'*Ef);
S(2)=(abs(Ef(1))).^2-(abs(Ef(2))).^2;
S(3)=2*real(Ef(1)*conj(Ef(2)));
S(4)=-2*imag(Ef(1)*conj(Ef(2)));
S(abs(S)<=eps)=0;
disp('------ Parametri di Stokes -------');
disp(['I = ' num2str(S(1))]);
disp(['Q = ' num2str(S(2))]);
disp(['U = ' num2str(S(3))]);
disp(['V = ' num2str(S(4))]);

%% PARAMETRI DELL'ELLISSE
% assi dell'ellisse
a=sqrt(1/2*(S(1)+sqrt(S(2)^2+S(3)^2)));
b=sqrt(1/2*(S(1)-sqrt(S(2)^2+S(3)^2)));

%cos2psi=S(2)./sqrt(S(1)^2-S(4)^2);
%ax_ang=acos(cos2psi)/2;

tan2psi=S(3)/S(2);
ax_ang=atan(tan2psi)/2;

%sin2psi = S(3)./sqrt(S(2)^2+S(3)^2);
%ax_ang = asin(sin2psi)/2;

ax_ang=ax_ang/pi*180;
if abs(Ef(1))<abs(Ef(2))
    ax_ang = ax_ang + 90;
end

disp('----- Parametri dell''ellisse -----');
disp(['Psi =' num2str(ax_ang)]);
disp(['a = ' num2str(a)]);
disp(['b = ' num2str(b)]);
%% ==== campo in uscita
disp('------- Campo in uscita ---------');
disp(['Ex = ',num2str(Ef(1))]);
disp(['Ey = ',num2str(Ef(2))]);
%% utilizzo del toolbox PhaseArray
   figure,polellip(Ef)
 %  [TAU,EPSILON,AR,RS] = polellip(Ef);
%% 03/03/2021 per il corso ILM facciamo il plot al variare dell'angolo di
%% un polarizzatore. Misura di intensità
theta = linspace(-pi/4,pi,360);
for i = 1:numel(theta)
    aa = P(theta(i))*Ef;
    I(i) = norm(aa)^2;%abs(aa(1))^2+abs(aa(2))^2;
end
thetad = theta*180/pi;
figure1 = figure;
plot(thetad,I),xlim([-45 180]),xlabel('\theta','FontSize',16),
ylabel('I_{riv}(\theta)','FontSize',16),grid
yline(norm(Ef)^2/2,'LineWidth',1.5)
xline(0,'LineWidth',1);
ka = 1.1;
kb = 0.8;
ylimb = b^2*kb;
ylima = a^2*ka;
ylim([ylimb ylima])
%% disegnamo le linee e punti interessanti
yline(a^2,'r','LineWidth',1);
yline(b^2,'b','LineWidth',1);
xline(ax_ang,'g','LineWidth',1);
hold on
%plot(45,0.5*(S(1) + S(2)),'*')
%% crea labels sul grafico
Dy = ylima-ylimb; %range y del grafico
ya = (a^2-ylimb)/Dy - 0.1;
yb = (b^2-ylimb)/Dy + 0.1;
ypsi = (mean(I)-ylimb)/Dy;
annotation(figure1,'textbox',...
    [0.84 ya 0.04 0.05],'String',{'a^2'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(figure1,'textbox',...
    [0.84 yb 0.04 0.05],...
    'String','b^2',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');
if ~isnan(ax_ang)
    annotation(figure1,'textbox',...
        [(ax_ang+48)/225 ypsi 0.04 0.05],...
        'String','\Psi',...
        'FontSize',16,...
        'FitBoxToText','off',...
        'EdgeColor','none');
end
