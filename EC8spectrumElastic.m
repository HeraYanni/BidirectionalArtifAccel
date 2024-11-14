function [ Se,Td,tc ] = EC8spectrumElastic(ag,soil,Td,damp)
%function [ Se,t1 ] = EC8spectrumElastic( ag,soil,t1,damp)
%
%   damp = 5 (default)
%
%Example:
%EC8 responce spectrum
%tb=0.15;
%tc=0.50;
%td=2.0;
%S = 1.2;
%pga=0.40;
%damp=5;
%figDir=cd;
%figure(4);
%hold on; box on; grid on;
%xlabel('period T (sec)','Fontsize',22);
%ylabel('Spectral acceleration (g)','Fontsize',22);
%set(gca,'Fontsize',16)
%tb=0.2;
%tc=0.8;
%td=2.0;
%S = 1.35;
%pga = 0.3
%[ t1,psa ] = EC8spectrum( damp,pga,tb,tc,td,S );
%plot(t1,psa,'k-','LineWidth',3)
%fname='ec8spectrum';
%printScript

if size(ag,2)>size(ag,1)
   ag = ag'; 
end

if strcmp(soil,'A')
    tb = 0.15;
    tc = 0.40;
    td = 2.5;
    S = 1.0;    
elseif strcmp(soil,'B')
    tb=0.15;
    tc=0.50;
    td=2.50;
    S = 1.2;
elseif strcmp(soil,'C')
    tb=0.20;
    tc=0.60;
    td=2.50;
    S = 1.15;
elseif strcmp(soil,'D')
    tb=0.20;
    tc=0.80;
    td =2.5;
    S = 1.35;
elseif strcmp(soil,'E')
    tb=0.15;
    tc=0.50;
    td=2.50;
    S = 1.4;
end

if ~exist('damp','var')
    damp = 5;
end


if damp > 1
    eta=sqrt(10./(5+damp));
else
    eta=sqrt(0.1/(0.05+damp));
end

% if ~exist('Td','var')
%     Td=linspace(0,4,1000);
% end

i1 = find(Td<=tb);
i2 = find(Td>tb & Td <=tc);
i3 = find(Td>=tc & Td <=td);
i4 = find(Td>td);

Se(i1) = ag*S*(1 +(Td(i1)/tb)*(eta*2.5 - 1));
Se(i2) = ag*S*2.5*eta;
Se(i3) = (ag*S*2.5*eta).*(tc./Td(i3));
Se(i4) = (ag*S*2.5*eta).*(tc*td)./(Td(i4).^2);

end