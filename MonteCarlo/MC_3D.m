function [out,launched,SStep,CCost,PPhi] =  MC_3D(mus,g,thick,N,TYPESCA,PLOT)
%% Simple MonteCarlo for light transport
% No refractive-index mismatch
% zero absorption
% Heyney-Greenstein or Raylegh scattering functions
% mus:      scattering coefficient
% g:        anisotropy factor g = <cos(theta)>
% thick:    thickness of the medium
% N:        number of photons to receive
% TYPESCA:  'HG' o 'RAY'
% PLOT:     1 o 0 for trajectory plot (use with only 1/2 photons)
% -------------------------------------------------------------------
% out:      coordinates (x,y) of the output photons and pathlength
% launched: launched photons numbers
% SStep:    steplengths of interactions
% CCost:    cosine of the azimuthal angles
% PPhi:     Zenithal angles
%% define scattering function
switch lower(TYPESCA)
    case 'hg'
        Spintheta = @()SpinthetaHG;
    case 'ray'
        Spintheta = @()SpinthetaRAY;
end
%% set maximum length for trajectories
LMAX = 1000*1/(mus*(1-g));

%% initialize some variables
out = zeros(N,3); %x,y,pathlength
det = 0;
launched = 0;
i_int = 1;

while det<N
%% launch a photon at p0 with direction mu0
launched = launched + 1;
photon.weight = 1;
photon.mu = [0;0;1];
photon.P = [0;0;0];
photon.path = 0;

if PLOT
    initPlot(photon.P);
end

while photon.weight
    s = stepsize();
    SStep(i_int) = s;
    i_int = i_int + 1;
    hit = HitBoundary(photon.P,photon.mu,s);
    if hit==1               % the photon exits from the detection side
        s = HopToBoundary(photon.P,photon.mu);
        photon.P = Hop(photon.P,photon.mu,s); % move the photon
        photon.path = photon.path + s;
        photon.weight = 0;         
        det = det + 1;
        % save output xy coordinates
        out(det,1:2) = photon.P(1:2);
        out(det,3) = photon.path;
    elseif hit==0           % the photon is inside
        P_new = Hop(photon.P,photon.mu,s); % move the photon
        photon.path = photon.path + s;
        if photon.path > LMAX
            photon.weight = 0;
        end
    elseif hit == 2 % the photon exits from the source side
        %clf;
        photon.weight = 0;
    end
    if PLOT
        PlotTrajectories(photon.P,P_new);
    end
    if photon.weight == 1
        photon.P = P_new;
        costheta = Spintheta(); % sample azimuthal angle
        phi = Spinphi();        % sample zenithal angle
        photon.mu = Spin(photon.mu,costheta,phi); % calculate new direction
    end
end
end
%% launched photons
disp(['Launched photons: ',num2str(launched)]);
disp(['Total transmittance = ',num2str(det./launched)]);
%% discard zeros elements
SStep(SStep==0) = [];
CCost(CCost==0) = [];
PPhi(PPhi==0) = [];
%% =======================================================================
function s = stepsize()
% sample the new step-length
    s = -log(rand())./mus;
end

function cost = SpinthetaHG()
% sample the new azimuthal direction
% Heyney-Greenstein phase function
    if g == 0
        cost = 2*rand() - 1;
    else        
        cost = 1./(2*g)*(1+g^2-((1-g^2)./(1-g+2*g*rand()))^2);
    end
    CCost(i_int) = cost;
end

function cost = SpinthetaRAY()
% sample the new azimuthal direction
% Rayleigh phase function
p = [1/8;0;3/8;1/2-rand()];
c = roots(p);
cost = c(abs(imag(c))<eps);
CCost(i_int) = cost;
end

function phi = Spinphi()
% sample the new zenithal direction
    phi = 2*pi*rand();
    PPhi(i_int) = phi;
end

function mup = Spin(mu,cost,phi)
% calculate the new direction
mup = zeros(3,1);
if abs(mu(3))>0.99999
    mup(1) = sqrt(1-cost^2)*cos(phi);
    mup(2) = sqrt(1-cost^2)*sin(phi);
    mup(3) = sign(mu(3))*cost;
else
    A = sqrt(1-cost^2)./sqrt(1-mu(3)^2);
    mup(1) = A.*(mu(1)*mu(3)*cos(phi)-mu(2)*sin(phi)) + mu(1)*cost;
    mup(2) = A.*(mu(2)*mu(3)*cos(phi)+mu(1)*sin(phi)) + mu(2)*cost;
    mup(3) = -sqrt(1-cost^2)*cos(phi)*sqrt(1-mu(3)^2)+mu(3)*cost;
end
end

function p = Hop(p_old,mu,s)
% move the photon
 p = p_old + mu*s;   
end

function hit = HitBoundary(p_old,mu,s)
% check if the photon hits the boundary
    hit = 0; % photon inside
    p = Hop(p_old,mu,s);
    if p(3) > thick
        hit = 1; % photon detected
    elseif p(3) < 0
        hit = 2; % photon exits and not detected
    end
end

function s = HopToBoundary(p,mu)
% calcualte a reduced step length to detection boundary
    s = (thick - p(3))./mu(3);        
end
%% visualization functions
function initPlot(P)
    L = 10;
    figure(1),subplot(2,2,1),
    plot3(P(1),P(2),P(3),'o'),
    xlim([-L,L]),ylim([-L,L]),zlim([0 thick]),
    xlabel('x'),ylabel('y'),zlabel('z'),hold on
    
    subplot(2,2,2)
    plot(P(1),P(3),'o'),
    xlim([-L,L]),ylim([0,thick])
    xlabel('x'),ylabel('z'),hold on
    
    subplot(2,2,3)
    plot(P(1),P(2),'o'),
    xlim([-L,L]),ylim([-L,L])
    xlabel('x'),ylabel('y'),hold on
end

function PlotTrajectories(P1,P2)
    subplot(2,2,1),plot3([P1(1);P2(1)],[P1(2);P2(2)],[P1(3);P2(3)],'.-'),
   grid on
   subplot(2,2,2),plot([P1(1);P2(1)],[P1(3);P2(3)],'.-'),
   grid on
   subplot(2,2,3),plot([P1(1);P2(1)],[P1(2);P2(2)],'.-'),
   grid on
 drawnow;
end

end