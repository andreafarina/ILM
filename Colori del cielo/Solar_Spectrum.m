% This scripts shows a typical solar spectrum and simulates the color
% perception fo the sky during daylight and sunset.
% The solar spectrum is multiplied:
% 1) by 1 - exp(-opt_thickness * Rayleigh scattering amplitude) for daylight
%    and cloudy.
% 2) by exp(-mu_s * opt_thick) in case of sunset (line of sight towards the
% sun)
% the resulting spectrum is converted to the XYZ tristimuls space thanks to
% the colour matching funcion. Finally the XYZ colour is converted to RGB
% and plotted as a uniform image.
% PRESS ANY KEY AFTER EACH FIGURE

% Author:       Andrea Farina
% Institution:  CNR - IFN
% email:        andrea.farina@polimi.it 
% March 2021; Last revision: 17-March-2022

clearvars;
close all
%% TYPE
choice = menu('Choice:','daylight','sunset','cloudy');
switch choice
    case 1
        type = 'daylight';
    case 2
        type = 'sunset';
        answer = inputdlg('Enter optical thickness (typ. value 10)');
        opt_thick = str2double(answer{1});
    case 3
        type = 'cloudy';
        answer = inputdlg('Enter optical thickness (typ. value 10)');
        opt_thick = str2double(answer{1});
end
%% load solar spectrum
a = dlmread('spectral_solar.txt');
lambda = a(:,3);    % nm
y = a(:,7);         % W/m^2/nm
LMAX = 1050;
lambda = lambda(1:LMAX);
y = y(1:LMAX);

%% visible range
lambda_VIS = [400,700];

%% plot
FONTSIZE = 14;
f = figure('Name','Solar Spectrum');%('Position',get(0,'ScreenSize'));
p1 = plot(lambda,y);
xlabel('Wavelength (nm)');
ylabel('arb units');
f.CurrentAxes.FontSize = FONTSIZE;
p1.LineStyle = '-';
p1.LineWidth = 3;
p1.Color = 'k';

%% PLOT RAYLEIGH
LAMBDA0 = 300;
switch lower(type)
    case 'daylight'
        Srayleigh = 1 - exp(-0.4*(LAMBDA0./lambda).^4*max(y)); % for 0.4 see
        % Bohren - Atmospheric Optics
    case 'sunset'
        Srayleigh = exp(-opt_thick*(LAMBDA0./lambda).^4*max(y));
    case 'cloudy'
        Srayleigh = 1- exp(-opt_thick*(LAMBDA0./lambda).^4*max(y)); 
end

hold on
plot(lambda,Srayleigh,'-b','LineWidth',2),ylim([min(y) max(y)])

%% load photopic curve
a = dlmread('ssvl2e_1.txt');
lambda_p = a(:,1);
y_p = a(:,2)./max(a(:,2))*max(y);
plot(lambda_p,y_p,'-g','LineWidth',2),ylim([min(y) max(y)]),
xline(lambda_VIS(1),'--p','LineWidth',2),xline(lambda_VIS(2),'--r','LineWidth',2);
legend('solar spectrum',['Rayleigh (',type,')'],'Photopic','violet','red')
pause;

%% calculate the diffused spectrum
y_diff = Srayleigh.*y;
f = figure('Name','Diffused spectrum');
    p1 = plot(lambda,y_diff);
    xlabel('Wavelength (nm)');
    ylabel('arb units');
    f.CurrentAxes.FontSize = FONTSIZE;
    p1.LineStyle = '-';
    p1.LineWidth = 3;
    p1.Color = 'k';
    xline(lambda_VIS(1),'--p','LineWidth',2),xline(lambda_VIS(2),'--r','LineWidth',2);
legend('Diffused spectrum','violet','red'),hold on
pause;
%% calculate the perceived colour
% requires: Get_xyz() Matlab AddOn for the colour matching-functions
% color-matching function
[X,Y,Z,~,~] = Get_xyz();
lambda_XYZ = X(:,1);

% interpolate the diffused spectrum over the initial lambda axis
y_diff_interp = interp1(lambda,y_diff,lambda_XYZ);

% matching function matrix
XYZ = [X(:,2),Y(:,2),Z(:,2)]';

%% plot color-matching 
area(lambda_XYZ,X(:,2),'FaceColor','r','FaceAlpha',0.3,'EdgeColor','r','LineWidth',2)
area(lambda_XYZ,Y(:,2),'FaceColor','g','FaceAlpha',0.3,'EdgeColor','g','LineWidth',2)
area(lambda_XYZ,Z(:,2),'FaceColor','b','FaceAlpha',0.3,'EdgeColor','b','LineWidth',2)
legend('Diffused spectrum','violet','red','X','Y','Z'),hold on
pause;

%% scalar product of the XYZ functions and the diffused spectrum
S_xyz = XYZ*y_diff_interp;
S_xyz = S_xyz./norm(S_xyz);

%% create RGB flat image
Npixels = 1024;
RGB = xyz2rgb(S_xyz');
RGB_image = reshape(RGB'*ones(1,Npixels^2),[3,Npixels,Npixels]);
figure('Name','The sky!'),imshow(permute(RGB_image,[2,3,1])),title('The sky!','FontSize',16)
