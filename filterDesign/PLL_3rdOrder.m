clc
clear
% PLL circuit block parameters
ref_freq =2.875e+9; % reference frequency
Icp =150e-6;
Kpd = Icp/(2*pi);
Kvco =50e6;
N =8;
C1 =5.4e-12;
C2 =54.5e-12;
C3 =54.5e-12;
R2 =34.38e+3;
R3 =34.38e+3;
order_lpf =3; % order of the LPF
dsm_order =2; % order of the DSM
bit_in =12; % input bit number of DSM.
% This is a 1-bit - output DSM
N_freq_multi =6; % frequency multiplication ratio after VCO
f1_vco =23e+9; % VCO starting frequency
f2_vco = f1_vco + FMCW_bandwidth / N_freq_multi ;% VCO ending frequency
N0= f1_vco / fref ; % starting division ratio =8
N1= f2_vco / fref ; % ending division ratio =8.4638



ref_noise_on =1; % noise of reference 1= on 0= off


fref = ref_freq ;
wref =2* pi* fref ;
PFD_cycle =1/ref_freq ;

% Choose to use divide -by -8/10 prescaler .
% Following designs are based on such prescaler .
alpha =(N1 -N0)/ T_chirp ; % slope of division ratio
q_N0 = ceil ((N0 -9+0.0001) *2^( bit_in -1) );
q_N1 = ceil ((N1 -9) *2^( bit_in -1));
N_step =q_N1 - q_N0 +1; % number of frequency steps
T_step = T_chirp / N_step ;
N_PFD_cycle = T_chirp / N_step / PFD_cycle ;% number of PFD cycles in one chirp step
freq_resolution = fref /2^ bit_in ;% frequency step / resolution
% simulation parameters for phase noise
fstart = 10e3;
fstep = 1e3;
fend = 100e6;
foff = fstart : fstep : fend ;% offset frequency from the carrier
s = 1i*2* pi* foff ;
Kboltzmann = 1.38e-23;
T_room = 300;

% FMCW chirp parameters
FMCW_bandwidth =8e+9;
chirp_slope =200*1e+6/1e-6;
T_chirp = FMCW_bandwidth / chirp_slope ;
n_ref_floor =10^( -140/10) ;
fc_ref =2.8184e+6;
phi_ref = sqrt ( n_ref_floor *( fc_ref ./ foff +1) );


phi_e_mean =0.363; % from linearity calculation
phi_e_mean_degree = phi_e_mean /pi *180;
alpha0 =6.6125/100;
ncurr0 =3.5e-12/ sqrt ( alpha0 *2);
fc0_ncurr =8.3111e+6/( alpha0 *2) ;
if abs ( phi_e_mean )<pi
duty_cycle = abs ( phi_e_mean )/2/ pi;
elseif abs ( phi_e_mean )>pi && abs ( phi_e_mean ) <2* pi
duty_cycle =(2* pi - abs ( phi_e_mean ))/2/ pi;
end
n_floor = ncurr0 * sqrt ( alpha0 *2+ duty_cycle );
fc_pfdcp = fc0_ncurr *( alpha0 *2+ duty_cycle );
ncurr_pfdcp = n_floor * sqrt ( fc_pfdcp ./ foff +1) ;% noise current
% spectrum , in A/ sqrt (Hz)

A0=C1+C2+C3;
A1=C2*R2 *( C1+C3)+C3*R3 *( C1+C2);
A2=C1*C2*C3*R2*R3;
A3 =0;
Zlpf =(1+ s*R2*C2)./(s .*( A0+s*A1+s .^2* A2+s .^3* A3));

ncurr_filter = sqrt (4* Kboltzmann * T_room .* real (1./ Zlpf ));
% above : noise current of the filter
% VCO
PN_1M = -105; % phase noise @ 1 MHz [ dBc ]
fc_vco = 10e6; % 1/f^3 cut frequency
fout = 23e9; % output frequency
f1M =1e6;
n_floorVCO = 10^( PN_1M /10) *f1M ^2/(1+ fc_vco / f1M );
phi_vco = sqrt ( n_floorVCO *(1./ foff .^2+ fc_vco ./ foff .^3) );
% DSM
fs= fout /N;
nv_dsm = sqrt((2* pi) ^2/12/ fs *(2* sin (pi* foff ./ fs)) .^(2*( dsm_order-1) ));



%% Transfer functions

Gfwd = Icp * Kvco * Zlpf ./s;
Grev =1/N;
T= Gfwd * Grev ;
Gnref =N* Gfwd ./(N+ Gfwd );
% PFD and CP
Gicp =2* pi/ Icp * Gnref;
% LPF
Glpf = (2* pi/ Icp ).* Gnref ;
% VCO
Gvco =1./(1+ T);
% divider
Gdiv =- Gnref ;

