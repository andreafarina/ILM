% ====== pi calculation by a Monte Carlo itegration ==========
%
%   pi=4*A_cirlce/A_square
%
%
%
clear all
close all
c=0;
d=0;
%i=0;
N=50000; %MAX 50000

hold on
X=0:0.01:2;
circ_up=sqrt(1-(X-1).^2)+1;
circ_down=-sqrt(1-(X-1).^2)+1;
subplot(1,2,1),hold on
plot(X,circ_up,'r',X,circ_down,'r'),xlabel('X','FontSize',16),...
    ylabel('Y','FontSize',16),grid,axis square

pi_est=zeros(1,N);
for i=1:N
    P=2*rand([1 2]);    %generic point in the unit square
    
    if (P(1)-1)^2+(P(2)-1)^2<1
        c=c+1;
        In(c,:)=P;
        %plot(P(1),P(2),'.r'),grid
    else
        d=d+1;
        Out(d,:)=P;
        %plot(P(1),P(2),'.k'),grid
        
    end
    
    pi_est(i)=4*c/i;
end

plot(In(:,1),In(:,2),'.r','MarkerSize',1);

plot(Out(:,1),Out(:,2),'.k','MarkerSize',1),grid
subplot(1,2,2),hold on
plot((1:N),repmat(pi,1,N),'r');
subplot(1,2,2)

plot((1:N),pi_est,'-'),grid,ylim([3. 3.3]),xlabel('Num trials','FontSize',16),...
    ylabel('\pi','FontSize',16)
%pi;
%pi_est(N);
ll=legend(['\pi=' num2str(pi)],['\pi_{ext}=' num2str(pi_est(N))]);
set(ll,'FontSize',16);