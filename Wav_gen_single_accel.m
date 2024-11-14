function[accel,PSa,fb_amor]=Wav_gen_single_accel(ut_G1,t,Sa_mean,Sa_,omega,PGA,N_iter,T1,zeta)

N = length(omega);
dom = omega(3)-omega(2);

dt = t(3)-t(2); % time step
L = t(end);

omega_0 = min(omega);
omega_u = max(omega);

%% Stationary accelrogram    

Ts = L; % accelerogram duration
% generate the stationary accelerogram
[a_G,~] = Stat_Accel_Gen(zeta,N,omega,Sa_mean,Ts,dom,t,Sa_,PGA,T1);

% CWT of the produced stationary artificial accelerogram
% Create filterbank
% Morlet wavelet = 'amor', Morse wavelet = 'morse'
f0=omega_0/(2*pi);
fu=omega_u/(2*pi);
fb_amor = cwtfilterbank_L2('SignalLength',numel(a_G),'SamplingFrequency',1/dt,'FrequencyLimits',[f0 fu],'VoicesPerOctave',48,'Wavelet','amor');
%fb_amor = cwtfilterbank_L2('SignalLength',numel(a_G),'SamplingFrequency',1/dt,'FrequencyLimits',[f0 fu],'VoicesPerOctave',48,'Wavelet','morse')

[wt_aGartif,f,~]=cwt_L2(a_G,'Filterbank',fb_amor);

[~,~,scls]=scales(fb_amor); % get the scales

%% Seed record

% CWT of the seed record
[wt_rec,~,~]=cwt_L2(ut_G1,'Filterbank',fb_amor);

%% Produce the new accelerogram

% parameters of Eq. 17
max_pga_rec=max(abs(ut_G1),[],'all');
max_pga_artif=max(abs(a_G),[],'all');

lambda = max_pga_artif./max_pga_rec; % Eq. 17

% parameter of Eq. 18
mod_aG = abs(wt_aGartif);

% parameters of Eq. 16
mod_rec = abs(wt_rec);
w_max=max(mod_rec,[],'all'); 

wt_norm = (mod_rec./w_max); % Eq. 16

% Produce the new accelerogram's CWT coefficients
for ik=1:length(f)
    E_aG(ik) = trapz(t,mod_aG(ik,:)); % Eq.18
    E_rec(ik) = trapz(t,mod_rec(ik,:)); % Eq 19

    WT(ik,:) = (E_aG(ik)./(lambda*E_rec(ik))).*wt_norm(ik,:).*wt_aGartif(ik,:); % Eq. 20
end

% Inverse CWT 
ia_G = 0;
ia_G = icwt_L2(WT,scls,'amor',f,[f0+(dom/2*pi) fu]); % the new artificial accelerogram
%ia_G = icwt_L2(WT,scls,'morse',f,[f0+(dom/2*pi) fu]);

[ut_G,PSa_aG]=FT_iter([t'; ia_G],N_iter,[T1;Sa_],T1,zeta,PGA);

% Baseline correction: 
% L: linear, Q: quadratic, C: cubic
[ut_G] = baseline_correction(ut_G','C');

accel = ut_G;
PSa=PSa_aG(N_iter,:);


end
