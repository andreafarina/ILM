% STUDIO DELLA POLARIZZAZIONE CON LE MATRICI DI MUELLER
% convenzione cos(kz-wt + phi) e guardare l'onda che ci arriva
% agli occhi.
% Le matrici di Mueller sono delle matrici 4x4 che permettono di trasformare
% i parametri di Stokes I,Q,U,V. Sono comode per gestire la luce
% polarizzata e parzialmente polarizzata.
% convenzioni: parallela - ortogonale
% 
% Author:       Andrea Farina
% Institution:  CNR - IFN
% email:        andrea.farina@polimi.it 
% March 2021; Last revision: 17-March-2022

close all;clearvars;
%% DEFINIZIONE DELLE MATRICI DI MUELLER
% Polarizzatore con asse parallelo alla direzione //
P_H = 0.5 *[1,1,0,0;
            1,1,0,0;
            0,0,0,0;
            0,0,0,0];

% rotatore angolo theta rispetto alla direzione //
Mt = @(theta)[1,0,0,0; 
              0,cos(2*theta), sin(2*theta),0;
              0,-sin(2*theta),cos(2*theta),0;
              0,0,0,1];

% matrice lamina ritardatrice con asse slow //  
% (quindi ritardo la componente parallela e l'asse fast e' nella direzione ortogonale)
LR = @(phi)[1,0,0,0;
            0,1,0,0;
            0,0,cos(phi),sin(phi);
            0,0,-sin(phi),cos(phi)];


%% COMPONENTI OTTICI
% polarizzatore con asse ottico orientato di un angolo alfa rispetto a //
P = @(alfa)Mt(-alfa)*P_H*Mt(alfa);

% lamina a l/4 con asse slow orientato di un angolo alfa rispetto a //
L4 = @(alfa)Mt(-alfa)*LR(pi/2)*Mt(alfa);

% lamina a l/2 con asse slow orientato di un angolo alfa rispetto a //
L2 = @(alfa)Mt(-alfa)*LR(pi)*Mt(alfa);

% lamina con sfasamento generico eta con asse slow orientato di un angolo
% alfa rispetto a //
Leta = @(eta,alfa)Mt(-alfa)*LR(eta)*Mt(alfa);

%% PROPAGAZIONE
% Le matrici di Mueller sono inserite in ordine inverso alla propagazione
%% Inserire qui il vettore dei parametri di Stokes del campo di partenza
% Es 1
%S0 = [41; 0; 9; 40];

% Es 2
S0 = [3; -1; 2; -2];

% luce non polarizzata
S0 = [1;0;0;0];

% luce polarizzata //
%S0 = [1;1;0;0];

% luce polarizzata circ dx
%S0 = [1;0;0;1];
%% Inserire qui la matrice di propagazione
% comporre in ordine inverso le funzioni implicite sopra nella sezione
% "Componenti Ottici"
M1 = P(-pi/4);
M2 = Leta(pi/3,pi/6);
M = M2*M1;
%% ========================================================================
% campo in uscita
S = M*S0;
S(abs(S)<eps) = 0;
disp('------ Parametri di Stokes -------');
disp(['I = ' num2str(S(1))]);
disp(['Q = ' num2str(S(2))]);
disp(['U = ' num2str(S(3))]);
disp(['V = ' num2str(S(4))]);

%% Dai parametri di Stokes all'ellisse
Ef(1) = sqrt(0.5*(S(1) + S(2)));
Ef(2) = sqrt(0.5*(S(1) - S(2)));
% delta necessita uno studio caso per caso per i 4 quadranti
if S(3)==0
    delta = pi/2*sign(S(4)); %+/- pi/2
end

if S(4)==0
    delta = pi*sign(S(3)); %+/- pi
end

if (S(3)>0 && S(4)>0)   % primo quadrante (no problem)
    delta = asin(0.5*S(4)./(Ef(1)*Ef(2)));
end

if (S(3)<0 && S(4)>0)   % secondo quadrante
    delta = acos(0.5*S(3)./(Ef(1)*Ef(2))); % uso acos tra [0,pi]
end

if (S(3)<0 && S(4)<0)    % terzo quadrante
    delta = -acos(0.5*S(3)./(Ef(1)*Ef(2)));
end

if (S(3)>0 && S(4)<0)   % quarto quadrante
    delta = asin(0.5*S(4)./(Ef(1)*Ef(2)));
end
    

Ef(2) = Ef(2) * exp(1i*delta);
Ef = Ef';

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
  figure,polellip(Ef),axis image
  [TAU,EPSILON,AR,RS] = polellip(Ef);
%% per il corso ILM facciamo il plot al variare dell'angolo di
%% un polarizzatore. Misura di intensità
theta = linspace(-pi/4,pi,360);
for i = 1:numel(theta)
    aa = P(theta(i))*S;
    I(i) = aa(1);   % L'intensità è il primo elemento del vettore di Stokes
end
thetad = theta*180/pi;
figure1 = figure; plot(thetad,I),xlim([-45 180]),xlabel('\theta','FontSize',16),
ylabel('I_{riv}(\theta)','FontSize',16),grid
yline(norm(Ef)^2/2,'LineWidth',1.5)
xline(0,'LineWidth',1);
ylim([b^2*(0.8) a^2*(1.1)])
%% disegnamo le linee e punti interessanti
yline(a^2,'r','LineWidth',1);
yline(b^2,'b','LineWidth',1);
xline(ax_ang,'g','LineWidth',1);
hold on
%plot(45,0.5*(S(1) + S(2)),'*')
%% crea labels sul grafico
annotation(figure1,'textbox',...
    [0.84 0.85 0.04 0.05],'String',{'a^2'},...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(figure1,'textbox',...
    [0.84 0.14 0.04 0.05],...
    'String','b^2',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(figure1,'textbox',...
    [(ax_ang+48)/225 0.44 0.04 0.05],...
    'String','\Psi',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');
