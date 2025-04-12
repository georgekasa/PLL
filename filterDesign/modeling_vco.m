clc
%clear
F0 = 640e6;
T0 = 1.0/F0;
L1_db = -142;%is for white noise
L1 = 10^(L1_db/10.0);
L2_db = -119;%20 db per decade
L2 = 10^(L2_db/10.0);
df2 = 10e6;
pn_1khz = -35;%1khz, 30 db per decade
sigma_flicker = (1e3/F0)*sqrt(T0)*sqrt(10^(pn_1khz/10.0));
sigma_t_NF = sqrt(F0*L1)/(2*pi*F0);
sigma_t_NF20 =(df2/F0) *sqrt(L2/F0);
cycles = 100e3;

% Filter parameters (corner frequencies in Hz)
f_corners = [1e3, 10e3, 100e3, 1e6, 10e6];
gains = [1, 0.3162, 0.1, 0.03162, 0.01];  % Corresponding gains
a = 2*pi.*f_corners./F0;
% Initialize output
flicker_noise = zeros(1, N);
designed_filters = cell(1, 5);
y1 = zeros(1,N);
y2 = y1;
y3 = y2;
y4 = y3;
y5 = y4; 
tv1 = 0;
tv2 = 0;
tv_test = 0;
CKV = zeros(1,cycles);
jitter_flicker = 0 ;
for i = 1:cycles
    if (i > 1)
        flicker_gaussian_noise = randn(1) * sigma_flicker;
        y1(i) = (1-a(1))*y1(i-1) + a(1)*gains(1)*flicker_gaussian_noise;
        y2(i) = (1-a(2))*y2(i-1) + a(2)*gains(2)*flicker_gaussian_noise;
        y3(i) = (1-a(3))*y3(i-1) + a(3)*gains(3)*flicker_gaussian_noise;
        y4(i) = (1-a(4))*y4(i-1) + a(4)*gains(4)*flicker_gaussian_noise;
        y5(i) = (1-a(5))*y5(i-1) + a(5)*gains(5)*flicker_gaussian_noise;
        jitter_flicker = y1(i) + y2(i) + y3(i) + y4(i) + y5(i) + jitter_flicker;
    end
    tv1 = i*T0 + randn(1)*sigma_t_NF;
   % dco_wander = sigma_t_NF20*randn(1);
    tv2 = tv2 + sigma_t_NF20*randn(1);
    CKV(i) = tv1 + tv2 + 0*jitter_flicker;

end

tdev = CKV - T0.*(1:cycles);
mean_tdev = mean(tdev);
std_tdev = std(tdev);
std_phase_deg = std_tdev / T0 * 360;

period = CKV(2:cycles) - CKV(1:cycles-1);
dT = 1 * (period - T0);
jitter_pp = (max(dT)-min(dT));
phase = 2 * pi * tdev / T0;

[Y,f] = pwelch(phase,blackmanharris(length(phase)/20),[],[],F0,'onesided');
YdB = 10 * log10(Y);
PN_SSB = YdB - 3;

semilogx(f, PN_SSB, 'b-','LineWidth',3)
%histogram(tdev, 'Normalization', 'count', 'EdgeColor', 'none', 'FaceColor', [0 0 0.5]);

% Labels and title
xlabel('Freq', 'FontSize', 12);
ylabel('PN ssb (dBc)', 'FontSize', 12);
title('Phase noise', 'FontSize', 14);
wannabe_cycle = sqrt((1.0/cycles)*sum((period - mean(period)).^2))