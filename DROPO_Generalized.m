clear all

% bright-bright

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-----general definition-----%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global c eps_0;  
c =  2.99792458E8;           % velocity of light in vacuum, m/s
eps_0 = 8.854E-12;           % dielectric constant in vacuum, F/m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-----------------Sim parameters---------%%%%%%%%%%%%%%%%%%%%%%%%
Nw = 1200;             % slicing number for time
Nzfi = 31;             % slicing number for PPLN1

%%%% OPO Configuration %%%%%%%%
lamp0 = 2100e-9;                  % wavelength of pump laser for OPO, m
lams0 = 2658.2e-9;                 % wavelength of signal laser for OPO, m
lami0 = 1/(1/lamp0-1/lams0);     % wavelength of idler laser for OPO, m
%%%% Nonlinearity and geometry... %%%%
Ac = pi*(3e-6)^2;                % m^2, mode area in PPLN1
L = 3E-3;                       % length of PPLN1, m
deff1 = 27.6E-12;                  % nonlinear coefficient for PPLN, m/V
%%%% Coupling and losses
alphap= pi/160/L;   alphas = alphap;    alphai = alphap; %loss 1/m
% Coupling Transmission set too 0 for non resonant condition
Rs = 1-pi/160;
Rp = 1-pi/160;
Ri = 1-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---------Pumping scheme--------%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ip0 = 500; Ip00 = 2*Ip0; % 00 for inital pulse seeding
Is0 = 0;             Is00 =2*Ip0;
Ii0 = 0;             Ii00 = 0;
noise = 1e-12; % Noise seed for non-resonant conditions
tseed = 1.0E-12;

% Detuning and phase-mismatch
detunep = -0.022784571428571*alphas*0+-0.149124671281650;
detunes = lamp0/lams0*detunep;
detunei = lamp0/lami0*detunep;

Nrt =2000;                     % roundtrip number
dk = 1*pi/L; % Phase mismatch

name = strcat('dk_',num2str(dk./pi*L),'pi___','Ip0_',num2str(Ip0));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%----------GVM and GDD of PPLN--------%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%---------GV-------%%%%%%%%% %%%%%%%%%---------GVM Compensation -------%%%%%%%%%
beta1_p = 7389e-12 + 0E-12;                              
beta1_s = 7389e-12 + 0E-12;
beta1_i = 7389e-12 + 0E-12;         % Compensation given in units of s for GVM and s^2 for GVD
beta_offset1 = beta1_p-beta1_s;      beta1_comp_p = (beta1_p-beta1_s).*L;
beta_offset2 = beta1_s-beta1_s;      beta1_comp_s = (beta1_p-beta1_s).*L;
beta_offset3 = beta1_i-beta1_s;      beta1_comp_i = (beta1_p-beta1_s).*L;
%%%%%%%---------GVD-------%%%%%%%%% %%%%%%%---------GVD Compensation-------%%%%%%%%%
beta2_p = 281e-27;                   beta2_comp_p = (0e-27).*L;
beta2_s = 213e-27;                   beta2_comp_s = (0e-27).*L;
beta2_i = - 1096e-27;                beta2_comp_i = (0e-27).*L;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

w_p = 2*pi*c/lamp0; 
w_s = 2*pi*c/lams0; 
w_i = 2*pi*c/lami0;

n_p = beta1_p*c;
n_s = beta1_s*c;
n_i = beta1_i*c;

kappa_p = sqrt(2)*w_p*deff1/sqrt(n_p*n_s*n_i*c^3*eps_0*Ac); 
kappa_s = sqrt(2)*w_s*deff1/sqrt(n_p*n_s*n_i*c^3*eps_0*Ac);
kappa_i = sqrt(2)*w_i*deff1/sqrt(n_p*n_s*n_i*c^3*eps_0*Ac);


%     kappa_p = 156.6925;
%     kappa_s = 65.587;
%     kappa_i = 91.1055;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tww= beta1_s*L*1;      % time window, related to the roundtrip time
t = linspace(-tww,tww,Nw);    
dt = mean(diff(t));                                                     
w = 2*pi*linspace(-1/2/dt,1/2/dt,Nw); % frequency 
dw = mean(diff(w));                                                                                     
z = linspace(0,L,Nzfi);  % PPLN1                                
dz = mean(diff(z));      % PPLN1                                                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-----------FFT transfer -----------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D_p = exp(-alphap/2.*(dz/2)).*exp(1j*(beta_offset1.*w + beta2_p/2.*w.^2)*dz/2); 
D_s = exp(-alphas/2.*(dz/2)).*exp(1j*(beta_offset2.*w + beta2_s/2.*w.^2)*dz/2);
D_i = exp(-alphai/2.*(dz/2)).*exp(1j*(beta_offset3.*w + beta2_i/2.*w.^2)*dz/2); 
% Dispersion Compensation
Dcomp_p = exp(-alphap/2.*(dz/2)).*exp(1j*(beta1_comp_p.*w + beta2_comp_p/2.*w.^2)*dz/2); 
Dcomp_s = exp(-alphas/2.*(dz/2)).*exp(1j*(beta1_comp_s.*w + beta2_comp_s/2.*w.^2)*dz/2);
Dcomp_i = exp(-alphai/2.*(dz/2)).*exp(1j*(beta1_comp_i.*w + beta2_comp_i/2.*w.^2)*dz/2); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% CW Pumping variables
Ap0 = sqrt(Ip0)*ones(1,Nw); 
As0 = sqrt(Is0)*ones(1,Nw); 
Ai0 = sqrt(Ii0)*ones(1,Nw); 

% Pulse Pump Variables
Ap00 = sqrt(Ip0)*ones(1,Nw) + sqrt(Ip00)*exp(-2*sqrt(2)*(t/tseed).^2);
As00 = sqrt(Is0)*ones(1,Nw) + sqrt(Is00)*exp(-2*sqrt(2)*(t/tseed).^2);  
Ai00 = sqrt(Ii0)*ones(1,Nw) + sqrt(Ii00)*exp(-2*sqrt(2)*(t/tseed).^2);  

% Initialize
Ap = noise*exp(2*pi*1j*rand(size(Ap0)));
As = noise*exp(2*pi*1j*rand(size(As0)));
Ai = noise*exp(2*pi*1j*rand(size(Ai0)));





%%%%%% Analytic guide %%%%%%
kappa_e = (beta2_s*kappa_i + beta2_i*kappa_s)./(kappa_s+kappa_i);
% Ap2 = sqrt((alphas.^2+detunes.^2)*(alphai.^2+detunei.^2))./(kappa_s*kappa_i*L.^2*sinc(dk*L/2).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LP = zeros(Nrt,Nw);
LS = zeros(Nrt,Nw);
LI = zeros(Nrt,Nw);
AZp = zeros(Nzfi,Nw);
AZs = zeros(Nzfi,Nw);
AZi = zeros(Nzfi,Nw);

wbd = 2*pi*9E12;
BD = fftshift(exp(-(w/wbd).^8));

% fftshift before stepping
D_p = fftshift(D_p);
D_s = fftshift(D_s);
D_i = fftshift(D_i);

Dcomp_p = fftshift(Dcomp_p);
Dcomp_s = fftshift(Dcomp_s);
Dcomp_i = fftshift(Dcomp_i);

for indrt = 1:Nrt
    Nrt - indrt

    if indrt<50
        Ap = (1-sqrt(Rp)).*Ap00 + sqrt(Rp)*Ap*exp(-1j*detunep) + noise*exp(2*pi*1j*rand(size(Ap)));
        As = (1-sqrt(Rs)).*As00 + sqrt(Rs)*As*exp(-1j*detunes) + noise*exp(2*pi*1j*rand(size(As)));
        Ai = (1-sqrt(Ri)).*Ai00 + sqrt(Ri)*Ai*exp(-1j*detunei) + noise*exp(2*pi*1j*rand(size(Ai)));
    else
        Ap = (1-sqrt(Rp)).*Ap0 + sqrt(Rp)*Ap*exp(-1j*detunep) + noise*exp(2*pi*1j*rand(size(Ap)));
        As = (1-sqrt(Rs)).*As0 + sqrt(Rs)*As*exp(-1j*detunes) + noise*exp(2*pi*1j*rand(size(As)));
        Ai = (1-sqrt(Ri)).*Ai0 + sqrt(Ri)*Ai*exp(-1j*detunei) + noise*exp(2*pi*1j*rand(size(Ai)));
    end

    % Propagation (1st half), split-step Fourier method and
    % DispersionCompensation
    sAp = BD.*(ifft(ifftshift(Ap)));
    sAs = BD.*(ifft(ifftshift(As)));
    sAi = BD.*(ifft(ifftshift(Ai)));
    Ap = fftshift(fft((D_p.*sAp.*Dcomp_p)));
    As = fftshift(fft((D_s.*sAs.*Dcomp_s)));
    Ai = fftshift(fft((D_i.*sAi.*Dcomp_i)));
    jj = 0;
    
    for indz = z

        jj = jj + 1;


        % nonlinear step using Runga-Kutta 4th order  
        [Ap, As, Ai] = OPOsol(Ap,As,Ai,kappa_p,kappa_s,kappa_i,dk,indz,dz);
        % Propagation (1st half), split-step Fourier method 
        sAp = BD.*(ifft(ifftshift(Ap)));
        sAs = BD.*(ifft(ifftshift(As)));
        sAi = BD.*(ifft(ifftshift(Ai)));
     % Propagation (2st half), split-step Fourier method
     if indz == z(end)
         Ap = fftshift(fft((D_p.^1.*sAp)));
         As = fftshift(fft((D_s.^1.*sAs)));
         Ai = fftshift(fft((D_i.^1.*sAi)));
     else
         Ap = fftshift(fft((D_p.^2.*sAp)));
         As = fftshift(fft((D_s.^2.*sAs)));
         Ai = fftshift(fft((D_i.^2.*sAi)));
     end

         AZp(jj,:) = Ap;
         AZs(jj,:) = As;
         AZi(jj,:) = Ai;

    end

    
    LP(indrt,:) = Ap;  
    LS(indrt,:) = As; 
    LI(indrt,:) = Ai;          
end
 

run('run_plots.m')
f1 = figure(1);
f2 = figure(2);
f6 = figure(6);
f6.Position = [2 463 1776 420];
savefig(f6,fullfile('Sim_Results',strcat(name,'_Time_Evolution.fig')))
savefig(f1,fullfile('Sim_Results',strcat(name,'_PeakPower_Evolution.fig')))
savefig(f2,fullfile('Sim_Results',strcat(name,'_OutputLzEvolution.fig')))
