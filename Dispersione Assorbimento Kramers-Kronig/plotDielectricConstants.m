% Plot della dell'indice di rifrazione e del coefficiente d'assorbimento
% di alcuni materiali usando il modello di Lorentz-Debye

% Apri LD.m per l'elenco di materiali e modelli disponibili.
% 
% Author:       Andrea Farina
% Institution:  CNR - IFN
% email:        andrea.farina@polimi.it 
% March 2021; Last revision: 17-March-2022

close all;
clearvars;

c = 3e8; %m/s
nu = logspace(-1,9,100000);     % GHz
lambda = c./nu;                 % nm
%% look at LD.m for options
[er, eim, n] = LD(lambda*1e-9,'H2O','LD');

% figure,subplot(2,2,1),loglog(nu,er),
% subplot(2,2,2),loglog(nu,eim)
% subplot(2,2,3),loglog(nu,real(n))
% subplot(2,2,4),loglog(nu,imag(n))

%% select visible area
lambda_VIS = [450,750]; %nm
nu_VIS = c./lambda_VIS;
figure(2),subplot(2,1,1),loglog(nu,er,'LineWidth',2),hold on
xline(nu_VIS(1)),xline(nu_VIS(2)),xlabel('Freq (GHz)'),
ylabel('\epsilon_r''')

subplot(2,1,2),loglog(nu,eim,'LineWidth',2),hold on
xline(nu_VIS(1)),xline(nu_VIS(2)),xlabel('Freq (GHz)'),
ylabel('\epsilon_r''''')

%% refractive index and absorption coefficient
figure(3),subplot(2,1,1),loglog(nu,real(n),'LineWidth',2),hold on
xline(nu_VIS(1)),xline(nu_VIS(2)),xlabel('Freq (GHz)'),
ylabel('refractive index')

mua = 2*2*pi*nu*1e9./c.*imag(n);
subplot(2,1,2),loglog(nu,mua,'LineWidth',2),hold on
xline(nu_VIS(1)),xline(nu_VIS(2)),xlabel('Freq (GHz)'),
ylabel('absorption')
% %return;
% %% use toolbox Communication
% material = 'wood';
% nu = nu*1e9;
% [~,~,epsilon] = arrayfun(@(y)buildingMaterialPermittivity(material,y),nu);
% figure, subplot(2,1,1),loglog(nu,real(epsilon)),
% subplot(2,1,2),loglog(nu,imag(epsilon))