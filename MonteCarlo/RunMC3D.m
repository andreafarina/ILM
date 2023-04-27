%% Script for run and read an MC_3D simulation
% the script shows:
% - the statistics of step-length, cos(theta), phi
% - trace of the output plane
% - pathlength histogram
close all;
PHASE = 'HG';   % HG: Heyniey-Greenstein, RAY: Rayleigh phase functions
PLOT = 1;       % when PLOT = 1, set N = 1 to see a single photon trajectory
N = 1;
mus = 10;
g = 0.5;
musp = (1-g)*mus;
thick = 5;
%% run the simulation
[out,~,SStep,CCost,PPhi] = MC_3D(mus,g,thick,N,PHASE,PLOT);
%% statistics of the step size
Nbins = 100;
figure,subplot(2,2,1),
h = histogram(SStep,Nbins); hold on
% theoretical function calculation
% midpoint between the bins edges
dz = diff(h.BinEdges(1:2));
z = mean(h.BinEdges(1:2)) + (0:Nbins-1)*dz;
% exponential distribution
f = mus*exp(-mus*z)*dz*sum(h.BinCounts);
plot(z,f,'LineWidth',2),title(['<s> = ',num2str(mean(SStep))])
xlabel('$s$','FontSize',16,'interpreter','latex'),legend('simulation','theoretical')
%% statistics of  cos(theta) azimuthal angle
subplot(2,2,2)
h = histogram(CCost,Nbins); hold on
% theoretical function calculation
% midpoint between the bins edges
dcost = diff(h.BinEdges(1:2));
cost = mean(h.BinEdges(1:2)) + (0:Nbins-1)*dcost;
switch lower(PHASE)
    case 'hg'
        % HG function
        f = (1-g^2)./(2.*(1 + g^2 -2*g*cost).^(1.5))*dcost*sum(h.BinCounts);
    case 'ray'
        % Rayleigh
        f = 2*pi*3/(16*pi).*(1 + cost.^2)*dcost*sum(h.BinCounts);
end
plot(cost,f,'LineWidth',2),title(['g = <cos\theta> = ',num2str(mean(CCost))],...
    'interpreter','tex'),
xlabel('$cos\theta$','FontSize',16,'interpreter','latex'),legend('simulation','theoretical')
    
%% statistics of phi zenithal angle (uniform)
subplot(2,2,3)
h = histogram(PPhi,Nbins); hold on
% theoretical function calculation
% midpoint between the bins edges
dphi = diff(h.BinEdges(1:2));
phi = mean(h.BinEdges(1:2)) + (0:Nbins-1)*dphi;
% uniform
f = 1/(2*pi)*dphi*sum(h.BinCounts)*ones(Nbins,1);
plot(phi,f,'LineWidth',2),title(['<\phi> = ',num2str(mean(PPhi))]),
xlabel('$\phi$','FontSize',16,'interpreter','latex'),legend('simulation','theoretical')
%% plot of the pathlength
Nbins = 100;
subplot(2,2,4),h = histogram(out(:,3),Nbins);
title(['<pathlength> = ', num2str(mean(out(:,3)))]),
xlabel('$\ell$','FontSize',16,'interpreter','latex');
%% plot of the output plane
Nbins_xy = 50;
F = 1;
% xB = linspace(-musp*thick*F,musp*thick*F,Nbins_xy);
% yB = linspace(-musp*thick*F,musp*thick*F,Nbins_xy);
xB = linspace(-10*F,10*F,Nbins_xy);
yB = linspace(-10*F,10*F,Nbins_xy);

[N,x,y] = histcounts2(out(:,1),out(:,2),xB,yB);%Nbins_xy);
dx = x(2) - x(1);
dy = y(2) - y(1);
xx = mean(x(1:2)) + (0:Nbins_xy-1)*dx;
yy = mean(y(1:2)) + (0:Nbins_xy-1)*dy;
figure,imagesc(yy,xx,N),xlabel('y'),ylabel('x'),
% xlim([-min(abs(yy(1)),abs(yy(end))) min(abs(yy(1)),abs(yy(end)))]);
% ylim([-min(abs(xx(1)),abs(xx(end))) min(abs(xx(1)),abs(xx(end)))]);
xlim([xB(1) xB(end)]);
ylim([yB(1) yB(end)]);

