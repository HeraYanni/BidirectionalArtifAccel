clc; clearvars; close all;
set(groot,'defaultFigureColor','w')
set(groot,'defaultLineLineWidth',1)
g=9.81;

% Tool for generating bidirectional target spectrum-compatible 
% fully-nonstationary artificial seismic ground motions.

% Further details are provided in the following document:
%
% Yanni H., Fragiadakis M., and Mitseas I.P. 
% "Wavelet-based stochastic model for the generation of fully
% non-stationary bidirectional seismic accelerograms". 
% Earthquake Engineering and Structural Dynamics. In review.

% Version 1.0 created by Hera Yanni, first release: 14th of November, 2024. 

% Copyright (c) 2024
% Hera Yanni
% Structural Engineer NTUA, MSc in ADERS
% Ph.D. Candidate, Laboratory for Earthquake Engineering NTUA
% email: heragian@mail.ntua.gr, heragian@gmail.com 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Models used:

% GMM to obtain Sa(Ti,zeta) and sigma_ln
% The median and logarithmic standard deviation values are obtained from 
% the BSSA14 NGA-West2 model, based on the following:
%
% Boore DM, Stewart JP, Seyhan E, Atkinson GM. 
% NGA-West2 Equations for Predicting PGA, PGV, and 5% Damped PSA for 
% Shallow Crustal Earthquakes. Earthquake Spectra. 2014;30:1057-1085.
% https://doi.org/10.1193/070113EQS184M.

% The correlation between pairs of periods are obtained from the following:
%
% Baker JW, Jayaram N. Correlation of Spectral Acceleration Values from 
% NGA Ground Motion Models. Earthquake Spectra. 2008;24(1):299-317.
% https://doi.org/10.1193/1.2857544

% The generation of the artificial ground motions uses the 
% Yanni H., Fragiadakis M., and Mitseas I.P. 
% "Wavelet-based stochastic model for the generation of fully
% non-stationary bidirectional seismic accelerograms". 
% Earthquake Engineering and Structural Dynamics.

% The code can be modified to use other spectrum-compatible ground motion 
% generation methods.

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For the generation of the spectra Sax, Say
% A GMM for the geometric mean target response spectrum Sa(Ti,zeta) and
% the standard deviation sigma_ln (Eq.28)
% Ntrials = number of trials for the MCS procedure
% weights for the scoring procedure (Eq.31):
% w1 = weight for the Sa(Ti,zeta) vs. GMxy = sqrt(Sax.*Say) error
% w2 = weight for the Sa(Ti,zeta) vs. Sax deviation
% w3 = weight for the Sa(Ti,zeta) vs. Say deviation
% rhoSaxSay_lim = spectral correlation minimum limit, 0.95
% MSE = mean squared error limit for the Sa*(Ti,zeta) vs. 
% GMxy = sqrt(Sax.*Say) error

%% For the generation of the spectra Sax, Say
% f_0 = lowest frequency (Hz) of the frequency range of the generated signals 
% fmax = highest frequency (Hz) of the frequency range of the generated signals
% rec_aG_L = the seed record in time and acceleration columns format for
% the X component
% rec_aG_T = the seed record in time and acceleration columns format for
% the X component (may be the same with rec_aG_L)
% N_iter = number of corrective iterations

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start of User inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters for the geometric mean target response spectrum Sa*(Ti,zeta)
Mw = 6.5; % Mw = Moment Magnitude
Rjb = 10; % Rjb = Joyner-Boore distance (km)
% Fault_Type    = 0 for unspecified fault
%               = 1 for strike-slip fault
%               = 2 for normal fault
%               = 3 for reverse fault
% normal fault
Fault_Type = 2; 

Vs30 = 600; % Vs30 = shear wave velocity averaged over top 30 m in m/s
% region        = 0 for global (incl. Taiwan)
%               = 1 for California
%               = 2 for Japan
%               = 3 for China or Turkey
%               = 4 for Italy
region = 4;
% z1            = Basin depth (km); depth from the groundsurface to the
%                   1km/s shear-wave horizon.
%               = 999 if unknown
z1 = 999;
zeta = 0.05; % damping ratio

% number of trials for the MCS procedure
Ntrials = 1000;

% weights
w1 = 1;
w2 = 0.07;
w3 = 0.07;

% spectral correlation minimum limit
rhoSaxSay_lim = 0.95;

% target error of 
MSE = 5; % value in percentage

%% Parameters for ground motion generation

% Load the seed records
% L component
 rec_aG_L =  load('Kozani_1995_L.dat'); 
% T component
 rec_aG_T =  load('Kozani_1995_L.dat');
% Always convert the acceleration values to m/s^2 !!!
 ut_L = rec_aG_L(:,2).*g;     % convert L component to m/s2
 ut_T = rec_aG_T(:,2).*g;     % convert T component to m/s2

 t = rec_aG_L(:,1); % time

% Define the desired frequency range of the ground motions
f_0 = 1/(2*pi); % lowest frequency (Hz), can't be less than 0.36/(2*pi)
fmax = 125/(2*pi); % highest frequency (Hz)

% Define the desired number of corrective iterations
N_iter = 3; % close spectrum matching is desired in this case

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of user inputs. 
% Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generation of the target spectra Sax , Say

% Discretization for periods
N = 2000; 
omega = linspace(2*pi*f_0,2*pi/0.01,N)';
dom = mean(diff(omega));
T = flip(2*pi./omega);

% Define th geometric mean target response spectrum Sa*(Ti,zeta)
% PGA
[PGA sigma_ag] = gmpe_bssa_2014(Mw, 0, Rjb, Fault_Type, region, z1, Vs30);
ag = PGA*g;

% Ground Motion Model
for i=1:N
    [Sa_mean(i) sigma_ln(i)] = gmpe_bssa_2014(Mw, T(i), Rjb, Fault_Type, region, z1, Vs30);
end
Sa_mean = Sa_mean.*g;

% Calculate the correlation matrix (Eq. 29)
for i=1:N
    for j=1:N
        Ti = T(i);
        Tj = T(j);
        rho(i,j) = baker_jayaram_correlation(Ti, Tj).*(0.79-0.023.*sqrt(Ti.*Tj));
    end

end

% Generate the spectra (Eqs. 28, 30-31)
[Sa_x, Sa_y] = Gen_Sax_Say(T,Sa_mean,sigma_ln,rho,w1,w2,w3,Ntrials,rhoSaxSay_lim,MSE);

% Geometric mean of the generated spectra
GM_xy = sqrt(Sa_x.*Sa_y);
 

%% Produce the bidirectional spectrum compatible accelerograms

Sax = flip(Sa_x'); % m/s2
Say = flip(Sa_y'); % m/s2

dt = t(3)-t(2); % time step
Fs = 1/dt;   % sampling frequency
Fn = Fs/2;   % Nyquist frequency

omega_max = 2*pi*Fn; % maximum frequency
omega_m = 2*pi*fmax; % rad/sec

omega_u = min(2*pi*fmax, omega_max); % rad/sec
omega_0 = 2*pi*f_0; % rad/sec

% Cut-off frequency upper limit check
if omega_m > omega_max
    msgbox('Maximum frequency must be less than the Nyquist frequency')
      return
end

% Cut-off frequency lower limit check
if omega_0 < 0.36 
    msgbox('Minimum frequency must be larger than 0.36 rad/s')
      return
end

% Find the omega values within the desired frequency range
indx_artif = find(omega>omega_0 & omega<(omega_u+dom));
indx_artif = [1 
              indx_artif];
omega_artif = omega(indx_artif);

% The two components
nrecs=2;

% Temporal correlation check: rho_temp < 0.30 according to 
% NIST GCR 11-917-15
rho_temp = 10;

while rho_temp>0.30

for i=1:nrecs

ag_comp =0;
Sa1 = 0;
Sa_ = 0;
    if i==1
        ag_comp = Sa_x(1);
        Sa1 = flip(Sax(indx_artif));
        Sa_ = Sa_x;
        ut = ut_L;
    elseif i==2
        ag_comp = Sa_y(1);
        Sa1 = flip(Say(indx_artif));
        Sa_ = Sa_y;
        ut = ut_T;
    end

% Returns [accelerogram, response spectrum, cwt filterbank for the plots]
[accel(i,:),PSa(i,:),fb_amor]=Wav_gen_single_accel(ut,t,Sa1,Sa_,omega_artif,ag_comp,N_iter,T',zeta);

end

r = corrcoef(accel(1,:),accel(2,:));
rho_temp = abs(r(1,2))
 end

% Geometric mean of the response spectra of the two generated components 
GM_gen_xy = sqrt(PSa(1,:).*PSa(2,:));

% MSE between the geometric mean of the response spectra of the two 
% generated components vs the target spectrum
disp('MSE between Sa(Ti,zeta) and  GM_gen_xy = sqrt(Sax_gen.*Say_gen) is error_gen (%) = ')
error_gen = 100*sum(abs(GM_gen_xy-Sa_mean).^2)/N

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Write the results 
% tempx=accel(1,:);
% accel_x = [t, tempx'];
% writematrix(accel_x,'accel_x.txt','Delimiter',' ') 
% 
% tempy=accel(2,:);
% accel_y = [t, tempy'];
% writematrix(accel_y,'accel_y.txt','Delimiter',' ') 

%% PLOTS

% Plot the x component
figure
subplot(3,3,[1,4])
hold on; grid on; box on;
[fft1,FAm1] = FFTp(t(3)-t(2),accel(1,:)); 
plot(FAm1./g,2*pi*fft1,'b') 
ylim([0,omega_u])
xlim([0,max(FAm1./g)]);
set(gca,'XDir','reverse','FontSize',21,'FontName','Times New Roman')
xlabel('$|F(\omega)|$ [g]','FontSize',24,'interpreter','latex');
ylabel('$\omega$ [rad/s]','FontSize',24,'interpreter','latex');

subplot(3,3,[2,3,5,6])
hold on; grid on; box on;  
ylim([omega_0,omega_u])
xlim([0,t(end)]);
[wt_f1,f1]=cwt_L2(accel(1,:),'Filterbank',fb_amor);
pcolor(t,2*pi.*f1,abs(wt_f1));
shading interp
colormap('jet')
set(gca,'FontSize',21,'FontName','Times New Roman')
ylabel('$\omega$ [rad/s]','FontSize',24,'interpreter','latex');
xlabel('Time [s]','FontSize',24,'interpreter','latex');
s.EdgeColor = 'none';
hcb = colorbar('location','EastOutside');

subplot(3,3,7)
hold on; grid on; box on;
set(gca,'FontSize',21,'FontName','Times New Roman')
plot(T,Sa_x./g,'r-.','LineWidth',3);
plot(T,PSa(1,:)./g,'b','LineWidth',2);
xlabel('Period [s]','FontSize',24,'interpreter','latex');
ylabel('Spec. accel. [g]','FontSize',24,'interpreter','latex');
legend('${S^*_a}(T_i,\zeta)$','${S_a}(T_i,\zeta)$','FontSize',18,'interpreter','latex','Location','northeast')
xlim([0, 4])
    
subplot(3,3,[8,9])
title('\alpha_x(t) component');
hold on; grid on; box on;
plot(t,accel(1,:)./g,'b')
set(gca,'FontSize',21,'FontName','Times New Roman')
yline(0,'k')
xlim([0, t(end)])
ylim([-Sa_x(1)./g, Sa_x(1)./g])
xlabel('Time [s]','FontSize',24,'interpreter','latex');
ylabel('$\alpha_x(t)$ [g]','FontSize',24,'interpreter','latex');
yline(Sa_x(1)./g,'r--','LineWidth',2);
yline(-Sa_x(1)./g,'r--','LineWidth',2);


% Plot the y component
figure
subplot(3,3,[1,4])
hold on; grid on; box on;
[fft1,FAm1] = FFTp(t(3)-t(2),accel(2,:)); 
plot(FAm1./g,2*pi*fft1,'b') 
ylim([0,omega_u])
xlim([0,max(FAm1./g)]);
set(gca,'XDir','reverse','FontSize',21,'FontName','Times New Roman')
xlabel('$|F(\omega)|$ [g]','FontSize',24,'interpreter','latex');
ylabel('$\omega$ [rad/s]','FontSize',24,'interpreter','latex');

subplot(3,3,[2,3,5,6])
hold on; grid on; box on;  
ylim([omega_0,omega_u])
xlim([0,t(end)]);
[wt_f2,f2]=cwt_L2(accel(2,:),'Filterbank',fb_amor);
pcolor(t,2*pi.*f2,abs(wt_f2));
shading interp
colormap('jet')
set(gca,'FontSize',21,'FontName','Times New Roman')
ylabel('$\omega$ [rad/s]','FontSize',24,'interpreter','latex');
xlabel('Time [s]','FontSize',24,'interpreter','latex');
s.EdgeColor = 'none';
hcb = colorbar('location','EastOutside');

subplot(3,3,7)
hold on; grid on; box on;
set(gca,'FontSize',21,'FontName','Times New Roman')
plot(T,Sa_y./g,'r-.','LineWidth',3);
plot(T,PSa(2,:)./g,'b','LineWidth',2);
xlabel('Period [s]','FontSize',24,'interpreter','latex');
ylabel('Spec. accel. [g]','FontSize',24,'interpreter','latex');
legend('${S^*_a}(T_i,\zeta)$','${S_a}(T_i,\zeta)$','FontSize',18,'interpreter','latex','Location','northeast')
xlim([0, 4])
    
subplot(3,3,[8,9])
title('\alpha_y(t) component');
hold on; grid on; box on;
plot(t,accel(2,:)./g,'b')
set(gca,'FontSize',21,'FontName','Times New Roman')
yline(0,'k')
xlim([0, t(end)])
ylim([-Sa_y(1)./g, Sa_y(1)./g])
xlabel('Time [s]','FontSize',24,'interpreter','latex');
ylabel('$\alpha_y(t)$ [g]','FontSize',24,'interpreter','latex');
yline(Sa_y(1)./g,'r--','LineWidth',2);
yline(-Sa_y(1)./g,'r--','LineWidth',2);


% Plot the target spectrum vs the geometric mean of the produced spectra vs
% the geometric mean of the generated ground motions' spectra
figure()
title('Geometric means comparison');
hold on; grid on; box on;
set(gca,'FontSize',19,'FontName','Times New Roman')
plot(T,Sa_mean./g,'r','LineWidth',3);
plot(T,GM_xy./g,'b--','LineWidth',3);
plot(T,GM_gen_xy./g,'g-.','LineWidth',3);
xlim([0 , 4])
ylim([0 0.9])
xlabel('Period [s]','FontSize',22,'interpreter','latex');
ylabel('Spectral acceleration [g]','FontSize',22,'interpreter','latex');
legend('${S^*_a}(T_i,\zeta)$','$S_{a,GM}^*(T_i,\zeta)$','$S_{a,GM}(T_i,\zeta)$','FontSize',18,'interpreter','latex','Location','southwest')
set(gca, 'XScale','log')
set(gca, 'YScale','log')
   % print('SpectraGen.png','-dpng')

