% Plot della suscettività e verfica della trasformata di Hilbert.
%
% Author:       Andrea Farina
% Institution:  CNR - IFN
% email:        andrea.farina@polimi.it 
% March 2021; Last revision: 17-March-2022
%% 
close all
clearvars;

wp = 100;       % MHz
w0 = 1e3;   
gamma = 200;    % us   

x = linspace(w0/100,w0*100,10000009);
chi = wp.^2./(w0^2 - wp^2/3 - x.^2 - 1i*gamma*x);

chi_r = real(chi);
chi_im = imag(chi);

figure('Name','Suscettivita');
subplot(1,2,1),plot(x,chi_r,'LineWidth',2),xlim([w0/2,w0*2]),
    xlabel('\omega','FontSize',18),ylabel('Re\{\chi\}','FontSize',18),grid
subplot(1,2,2),plot(x,chi_im,'LineWidth',2),xlim([w0/2,w0*2]),xlabel('\omega','FontSize',18),
ylabel('Im\{\chi\}','FontSize',18),
grid

%% test with hilbert transform
fK1 = hilbert(chi_r);
figure('Name','Trasformate di Hilbert');
subplot(1,2,1),plot(x,imag(fK1),'LineWidth',2),xlim([w0/2,w0*2]),xlabel('\omega','FontSize',18),
ylabel('Hilb\{Re\{\chi\}\}','FontSize',18),grid,
fK2 = hilbert(chi_im);
subplot(1,2,2),plot(x,imag(-fK2),'LineWidth',2),xlim([w0/2,w0*2]),xlabel('\omega','FontSize',18),
ylabel('Hilb\{Im\{\chi\}\}','FontSize',18),grid,

%% test Kramers-Kronig
% wrange = linspace(w0-50,w0+50,100);
% 
% fK1 = zeros(size(wrange));
% fK2 = zeros(size(wrange));
% for i = 1:numel(wrange)
% w = wrange(i);
% %chi_val = wp.^2./(w0^2 - wp^2/3 - w.^2 - 1i*gamma*w)
% maskleft = x<w;
% maskright = x>w;
% xleft = x(x<w);
% xright = x(x>w);
% %figure(2),plot(x,chi_im.*x./(w^2-x.^2))
% fK1(i) = -2/pi*(trapz(xleft,chi_im(maskleft).*xleft./(w^2-xleft.^2)) +...
%     trapz(xright,chi_im(maskright).*xright./(w^2-xright.^2)));
% fK2(i) = 2/pi.*w.*(trapz(xleft,chi_r(maskleft)./(w^2-xleft.^2)) +...
%     trapz(xright,chi_r(maskright)./(w^2-xright.^2)));
% end
% figure(3),plot(wrange,fK1),
% figure(4),plot(wrange,fK2)
