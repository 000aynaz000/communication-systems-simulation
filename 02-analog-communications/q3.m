%% Part 1 - FM

clear; close all; clc;

% 1.1) Signal Generation
Fs = 200e3;                
dt = 1/Fs;
N  = 2^16;                 
t  = ((0:N-1) - N/2) * dt;

m = 2000 * tri_fun(2000 * t);  

% 1.2) Time- and frequency-domain plots
M  = fftshift(fft(m)) * dt;
f  = ((0:N-1) - N/2) * (Fs/N);     % Frequency axis (Hz)
MagM = abs(M);

% Time-domain plot over a suitable interval
figure('Name','1.2 Time Domain');
t_ms = t * 1e3;
mask_t = (t >= -1.5e-3) & (t <= 1.5e-3);
plot(t_ms(mask_t), m(mask_t), 'LineWidth', 1.5);
grid on;
xlabel('Time (ms)');
ylabel('m(t)');
title('m(t) = 2000 \cdot tri(2000t)');

% Frequency-domain magnitude plot and mark the first zero at f = 2 kHz
figure('Name','1.2 Frequency Domain');
mask_f = (f >= -10e3) & (f <= 10e3);
plot(f(mask_f)/1e3, MagM(mask_f), 'LineWidth', 1.5); hold on;
xline(2,  '--', 'First zero (2 kHz)', 'LabelHorizontalAlignment','left');
xline(-2, '--', 'First zero (-2 kHz)', 'LabelHorizontalAlignment','right');
grid on;
xlabel('Frequency (kHz)');
ylabel('|M(f)|');
title('Magnitude Spectrum |M(f)| (FFT-based)');

% 1.3) Effective bandwidth Beff (90% energy)
PSD = abs(M).^2;    
df  = Fs / N;

posMask = (f >= 0);
fpos = f(posMask);
Ppos = PSD(posMask);

E_total = (Ppos(1) + 2*sum(Ppos(2:end))) * df;

E_cum = Ppos(1)*df + 2*cumsum(Ppos(2:end))*df;    
ratio = E_cum / E_total;

idx = find(ratio >= 0.90, 1, 'first');  
Beff = fpos(idx+1);                     

fprintf('Estimated Beff (90%% energy) ≈ %.2f Hz (%.4f kHz)\n', Beff, Beff/1e3);

figure('Name','1.3 Effective Bandwidth');
plot(fpos(2:end)/1e3, ratio, 'LineWidth', 1.5); hold on;
yline(0.90, '--', '90% energy');
xline(Beff/1e3, '--', sprintf('Beff = %.3f kHz', Beff/1e3));
grid on;
xlabel('B (kHz)');
ylabel('Energy ratio in [-B, B]');
title('Effective Bandwidth B_{eff} for 90% Energy');

% 1.4) Lowpass filtering with cutoff approximately Beff
Wn = Beff / (Fs/2);         % Normalized cutoff (0..1)
Wn = min(max(Wn, 1e-6), 0.999999);  % Safety clamp

filtOrder = 400;           
b = fir1(filtOrder, Wn, 'low', hamming(filtOrder+1), 'scale');

mLP = filtfilt(b, 1, m);

MLP = fftshift(fft(mLP)) * dt;
MagMLP = abs(MLP);

figure('Name','1.4 Time Domain: original vs filtered');
plot(t_ms(mask_t), m(mask_t), 'LineWidth', 1.2); hold on;
plot(t_ms(mask_t), mLP(mask_t), 'LineWidth', 1.5);
grid on;
xlabel('Time (ms)');
ylabel('Amplitude');
legend('m(t)', 'm_{LP}(t)', 'Location','best');
title(sprintf('Lowpass Filtering with cutoff \\approx Beff = %.3f kHz', Beff/1e3));

figure('Name','1.4 Spectrum: original vs filtered');
plot(f(mask_f)/1e3, MagM(mask_f),   'LineWidth', 1.2); hold on;
plot(f(mask_f)/1e3, MagMLP(mask_f), 'LineWidth', 1.5);
xline(Beff/1e3, '--', 'Beff');
xline(-Beff/1e3, '--');
grid on;
xlabel('Frequency (kHz)');
ylabel('Magnitude');
legend('|M(f)|', '|M_{LP}(f)|', 'Location','best');
title('Spectrum Before/After Lowpass Filtering');

% function y = tri_fun(x)
%     ax = abs(x);
%     y = zeros(size(x));
%     idx = (ax < 1);
%     y(idx) = 1 - ax(idx);
% end



%% Part 2 – Direct FM Modulation (2.1 – 2.4)

% 2.1) Direct FM modulation using mLP(t)
fc_RF   = 10e6;      
N       = length(t);

t   = t(:);
mLP = mLP(:);

% 2.2) Selecting the frequency deviation (Δf)
fm_max = Beff;                
df_target = 10 * fm_max;       
df = min(max(df_target, 20e3), 40e3);  

beta = df / fm_max;

fprintf('Chosen fm_max = %.2f Hz (%.4f kHz)\n', fm_max, fm_max/1e3);
fprintf('Chosen Δf     = %.2f Hz (%.2f kHz)\n', df, df/1e3);
fprintf('Resulting beta= %.2f\n', beta);

x = mLP / max(abs(mLP));

dt = 1/Fs;
phi = 2*pi*df * cumsum(x) * dt;     
sBB = exp(1j * phi);                

fc_plot = 50e3;  
sFM = real( sBB .* exp(1j*2*pi*fc_plot*t) );

% 2.3) Time- and frequency-domain plots of FM signal
figure('Name','2.3 FM Time Domain (IF visualization)');
mask_t = (t >= -200e-6) & (t <= 200e-6);

plot(t(mask_t)*1e6, sFM(mask_t), 'LineWidth', 1.2);
grid on;
xlabel('Time (\mus)');
ylabel('s_{FM}(t) (real, IF)');
title(sprintf('FM signal at IF carrier f_c=%.1f kHz (RF carrier is 10 MHz)', fc_plot/1e3));

S   = fftshift(fft(sFM)) * dt;
f   = ((0:N-1)' - N/2) * (Fs/N);   
MagS = abs(S);


figure('Name','2.3 FM Magnitude Spectrum (around IF)');
span = (df + fm_max) * 2.5;        
mask_f = (f >= (fc_plot - span)) & (f <= (fc_plot + span));
plot(f(mask_f)/1e3, MagS(mask_f), 'LineWidth', 1.2); grid on;
xlabel('Frequency (kHz)');
ylabel('|S_{FM}(f)|');
title(sprintf('FM Spectrum around IF (f_c=%.1f kHz), Δf=%.1f kHz, fm,max≈%.2f kHz', ...
    fc_plot/1e3, df/1e3, fm_max/1e3));

% 2.4) Effective FM bandwidth (90% energy) and comparison to Carson’s rule

% Remove DC component from baseband FM signal
sBB_nodc = sBB - mean(sBB);

S = fftshift(fft(sBB_nodc)) * dt;
f = ((0:N-1)' - N/2) * (Fs/N);

PSD = abs(S).^2;
df_bin = Fs / N;

[f, idx] = sort(f);
PSD = PSD(idx);

posMask = (f >= 0);
fpos = f(posMask);
Ppos = PSD(posMask);

% Total energy
E_total = (Ppos(1) + 2*sum(Ppos(2:end))) * df_bin;

% Cumulative energy
E_cum = Ppos(1)*df_bin + 2*cumsum(Ppos(2:end))*df_bin;
ratio = E_cum / E_total;

% 90% bandwidth
idx90 = find(ratio >= 0.90, 1, 'first');
Beff_FM_half = fpos(idx90);
Beff_FM_full = 2 * Beff_FM_half;

% Carson bandwidth
BC = 2 * (df + fm_max);

fprintf('\n FM Bandwidth Results\n');
fprintf('Effective FM BW (90%% energy): half-width = %.2f kHz\n', Beff_FM_half/1e3);
fprintf('Effective FM BW (90%% energy): full BW   = %.2f kHz\n', Beff_FM_full/1e3);
fprintf('Carson BW approximation (full): BC       = %.2f kHz\n', BC/1e3);



%% Part 3 – FM Modulation with Armstrong Method (Indirect FM)

% 3.1) Narrowband FM at 500 kHz
Fs = 5e6;                 
dt = 1/Fs;
N  = 2^17;               
t  = ((0:N-1) - N/2) * dt;    
t  = t(:);

m = 400 * tri_fun(400 * t);

x = m / max(abs(m));

fm_max = 400;             

beta0 = 0.25;            
df0   = beta0 * fm_max;  

fprintf('--- Part 3 Setup ---\n');
fprintf('fm_max ≈ %.2f Hz\n', fm_max);
fprintf('Chosen narrowband beta0 = %.3f\n', beta0);
fprintf('Initial deviation Δf0   = %.2f Hz\n\n', df0);

phi0 = 2*pi*df0 * cumsum(x) * dt;

fbase_RF = 500e3;         

fc1_plot = 10e3;          

K = fbase_RF / fc1_plot; 

fprintf('Plot-carrier fc1_plot = %.2f kHz maps to RF fbase = %.2f kHz (scale K=%.1f)\n', ...
    fc1_plot/1e3, fbase_RF/1e3, K);

s1 = exp(1j * (2*pi*fc1_plot*t + phi0));   

% 3.2) Frequency multiplication
Nmult = 200;              
s_mult = s1 .^ Nmult;

fc2_plot = Nmult * fc1_plot;       
fc2_RF   = Nmult * fbase_RF;       

df_final = Nmult * df0;            
beta_final = Nmult * beta0;

fprintf('\n--- 3.2 Frequency Multiplication ---\n');
fprintf('Multiplier Nmult = %d\n', Nmult);
fprintf('Carrier after mult (plot): %.2f kHz\n', fc2_plot/1e3);
fprintf('Carrier after mult (RF):   %.2f MHz\n', fc2_RF/1e6);
fprintf('Deviation after mult Δf:   %.2f kHz\n', df_final/1e3);
fprintf('Beta after mult:           %.2f\n', beta_final);

% 3.3) Frequency translation to final carrier 

fc_final_RF   = 10e6;      
fc_final_plot = 200e3;     

fLO_plot = fc2_plot - fc_final_plot;

s_indirect_cx = s_mult .* exp(-1j * 2*pi*fLO_plot*t);

s_indirect = real(s_indirect_cx);

fprintf('\n--- 3.3 Frequency Translation ---\n');
fprintf('Final carrier (plot): %.2f kHz\n', fc_final_plot/1e3);
fprintf('Final carrier (RF):   %.2f MHz\n', fc_final_RF/1e6);
fprintf('Mixing LO (plot):     %.2f kHz\n', fLO_plot/1e3);
fprintf('In RF analogy: LO ≈ %.2f MHz (e.g., 100-10=90 MHz)\n', (fc2_RF - fc_final_RF)/1e6);

% 3.4) Time- and frequency-domain plots of indirect FM
figure('Name','3.4 Indirect FM Time Domain');
mask_t = (t >= -200e-6) & (t <= 200e-6);
plot(t(mask_t)*1e6, s_indirect(mask_t), 'LineWidth', 1.2);
grid on;
xlabel('Time (\mus)');
ylabel('s_{indirect}(t) (real)');
title(sprintf('Indirect FM (Armstrong) at plot carrier %.1f kHz (maps to RF 10 MHz)', fc_final_plot/1e3));

S  = fftshift(fft(s_indirect)) * dt;
f  = ((0:N-1)' - N/2) * (Fs/N);
MagS = abs(S);

BC_arm = 2*(df_final + fm_max);     
span = 2.5*(BC_arm/2);              
mask_f = (f >= (fc_final_plot - span)) & (f <= (fc_final_plot + span));

figure('Name','3.4 Indirect FM Spectrum');
plot(f(mask_f)/1e3, MagS(mask_f), 'LineWidth', 1.2); grid on;
xlabel('Frequency (kHz)');
ylabel('|S_{indirect}(f)|');
title(sprintf('Indirect FM Spectrum (plot carrier %.1f kHz), \\Delta f=%.1f kHz, fm,max≈%.2f kHz', ...
    fc_final_plot/1e3, df_final/1e3, fm_max/1e3));

% 3.4) Effective bandwidth (90% energy) for indirect FM
[Beff_half, Beff_full] = eff_bw_90(S, f, fc_final_plot, Fs);

fprintf('\n--- 3.4 Effective BW (Indirect FM) ---\n');
fprintf('Carson BW (full) BC_arm      = %.2f kHz\n', BC_arm/1e3);
fprintf('Effective BW (90%% energy) full= %.2f kHz\n', Beff_full/1e3);

figure('Name','3.4 Energy Ratio vs B (Indirect FM)');
plot(Beff_half.axis_khz, Beff_half.ratio, 'LineWidth', 1.2); hold on; grid on;
yline(0.90, '--', '90% energy');
xline(Beff_half.B90_khz, '--', sprintf('B_{90} = %.2f kHz', Beff_half.B90_khz));
xline((BC_arm/2)/1e3, '--', sprintf('B_C/2 = %.2f kHz', (BC_arm/2)/1e3));
xlabel('B (kHz) (half-width around carrier)');
ylabel('Energy ratio in [f_c-B, f_c+B]');
title('Indirect FM: Effective BW (90% energy) vs Carson half-width');

% 3.5) Comparison with direct FM (Part 2)
fprintf('\n--- 3.5 Comparison with Direct FM (Part 2) ---\n');

fprintf('INDIRECT (Armstrong): Δf = %.2f kHz, Carson BW = %.2f kHz, Eff(90%%) BW = %.2f kHz\n', ...
    df_final/1e3, BC_arm/1e3, Beff_full/1e3);

if exist('df_direct','var') && exist('fm_max_direct','var') && exist('Beff_FM_full_direct','var') && exist('BC_direct','var')
    fprintf('DIRECT (Part 2):      Δf = %.2f kHz, Carson BW = %.2f kHz, Eff(90%%) BW = %.2f kHz\n', ...
        df_direct/1e3, BC_direct/1e3, Beff_FM_full_direct/1e3);

    fprintf('\nSpectral shape note:\n');
    fprintf('- Ideal FM spectral shape depends mainly on beta and message spectrum.\n');
    fprintf('- If direct and indirect end with similar Δf and fm,max, spectra should be very similar.\n');
else
    fprintf(['Direct-FM comparison variables not found.\n' ...
             'To compare numerically, store Part 2 outputs like:\n' ...
             '  df_direct, fm_max_direct, BC_direct, Beff_FM_full_direct\n']);
end

function y = tri_fun(x)
    ax = abs(x);
    y = zeros(size(x));
    idx = (ax < 1);
    y(idx) = 1 - ax(idx);
end

function [outHalf, Beff_full] = eff_bw_90(S, f, fc, Fs)
    PSD = abs(S).^2;
    df_bin = Fs / length(f);

    f_rel = f - fc;

    [f_rel, ord] = sort(f_rel);
    PSD = PSD(ord);

    pos = (f_rel >= 0);
    fpos = f_rel(pos);
    Ppos = PSD(pos);

    E_total = (Ppos(1) + 2*sum(Ppos(2:end))) * df_bin;
    E_cum   = Ppos(1)*df_bin + 2*cumsum(Ppos(2:end))*df_bin;
    ratio   = E_cum / E_total;

    idx = find(ratio >= 0.90, 1, 'first');
    B90 = fpos(idx+1);                
    Beff_full = 2 * B90;              

    outHalf.B90_khz  = B90/1e3;
    outHalf.axis_khz = fpos(2:end)/1e3;
    outHalf.ratio    = ratio;
end


%% Part 4 – FM and PM through AWGN and SNR Analysis 

Fs = 5e6;                
dt = 1/Fs;
N  = 2^17;
t  = ((0:N-1) - N/2) * dt;  t = t(:);

m = 400 * tri_fun(400*t);
x = m / max(abs(m));       

fc_plot = 200e3;           

fm_max = 400;             
betaFM_target = 10;
dfFM = min(max(betaFM_target * fm_max, 20e3), 40e3); 

betaPM = 1.0;              

phiFM = 2*pi*dfFM * cumsum(x) * dt;
sFM = cos(2*pi*fc_plot*t + phiFM);

sPM = cos(2*pi*fc_plot*t + betaPM*x);

% 4.1) AWGN channel model with 50 dB loss and N0/2 = 0.5e-12
Loss_dB = 50;
a = 10^(-Loss_dB/20);       

N0_over2 = 0.5e-12;        
sigma2 = N0_over2 * Fs;     
sigma = sqrt(sigma2);

nFM = sigma * randn(N,1);
nPM = sigma * randn(N,1);

rFM = a*sFM + nFM;
rPM = a*sPM + nPM;

% 4.2) Time- and frequency-domain inspection
mask_t = (t >= -200e-6) & (t <= 200e-6);

figure('Name','4.2 Time: FM clean vs noisy');
plot(t(mask_t)*1e6, sFM(mask_t), 'LineWidth', 1.0); hold on;
plot(t(mask_t)*1e6, rFM(mask_t), 'LineWidth', 1.0);
grid on; xlabel('Time (us)'); ylabel('Amplitude');
legend('sFM(t)', 'rFM(t)', 'Location','best'); 
title('FM: clean vs received (50 dB loss + AWGN)');

figure('Name','4.2 Time: PM clean vs noisy');
plot(t(mask_t)*1e6, sPM(mask_t), 'LineWidth', 1.0); hold on;
plot(t(mask_t)*1e6, rPM(mask_t), 'LineWidth', 1.0);
grid on; xlabel('Time (us)'); ylabel('Amplitude');
legend('sPM(t)', 'rPM(t)', 'Location','best');
title('PM: clean vs received (50 dB loss + AWGN)');

[faxis, Mag_sFM] = mag_spectrum(sFM, Fs);
[~,     Mag_rFM] = mag_spectrum(rFM, Fs);
[~,     Mag_sPM] = mag_spectrum(sPM, Fs);
[~,     Mag_rPM] = mag_spectrum(rPM, Fs);

span = 250e3; 
mask_f = (faxis >= (fc_plot - span)) & (faxis <= (fc_plot + span));

figure('Name','4.2 Spectrum: FM clean vs noisy');
plot(faxis(mask_f)/1e3, Mag_sFM(mask_f), 'LineWidth', 1.0); hold on;
plot(faxis(mask_f)/1e3, Mag_rFM(mask_f), 'LineWidth', 1.0);
grid on; xlabel('Frequency (kHz)'); ylabel('Magnitude');
legend('SFM(f)', 'RFM(f)', 'Location','best'); 
title('FM Spectrum: clean vs received');

figure('Name','4.2 Spectrum: PM clean vs noisy');
plot(faxis(mask_f)/1e3, Mag_sPM(mask_f), 'LineWidth', 1.0); hold on;
plot(faxis(mask_f)/1e3, Mag_rPM(mask_f), 'LineWidth', 1.0);
grid on; xlabel('Frequency (kHz)'); ylabel('Magnitude');
legend('SPM(f)', 'RPM(f)', 'Location','best'); 
title('PM Spectrum: clean vs received');

% 4.3) Noise reduction using simple filtering
BC_FM = 2*(dfFM + fm_max);              
dfPM_equiv = betaPM * fm_max;           
BC_PM = 2*(dfPM_equiv + fm_max);        

Bhalf_FM = 0.6 * BC_FM;                 
Bhalf_PM = 0.6 * BC_PM;

ord = 600; 

bpFM = fir1(ord, [(fc_plot - Bhalf_FM)/(Fs/2), (fc_plot + Bhalf_FM)/(Fs/2)], 'bandpass', hamming(ord+1));
bpPM = fir1(ord, [(fc_plot - Bhalf_PM)/(Fs/2), (fc_plot + Bhalf_PM)/(Fs/2)], 'bandpass', hamming(ord+1));

rFM_tilde = filtfilt(bpFM, 1, rFM);
rPM_tilde = filtfilt(bpPM, 1, rPM);

figure('Name','4.3 Time: FM noisy vs filtered');
plot(t(mask_t)*1e6, rFM(mask_t), 'LineWidth', 1.0); hold on;
plot(t(mask_t)*1e6, rFM_tilde(mask_t), 'LineWidth', 1.0);
grid on; xlabel('Time (us)'); ylabel('Amplitude');
legend('rFM(t)', 'rFM filtered', 'Location','best'); 
title('FM: received vs filtered');

figure('Name','4.3 Time: PM noisy vs filtered');
plot(t(mask_t)*1e6, rPM(mask_t), 'LineWidth', 1.0); hold on;
plot(t(mask_t)*1e6, rPM_tilde(mask_t), 'LineWidth', 1.0);
grid on; xlabel('Time (us)'); ylabel('Amplitude');
legend('rPM(t)', 'rPM filtered', 'Location','best'); 
title('PM: received vs filtered');

[~, Mag_rFMf] = mag_spectrum(rFM_tilde, Fs);
[~, Mag_rPMf] = mag_spectrum(rPM_tilde, Fs);

figure('Name','4.3 Spectrum: FM noisy vs filtered');
plot(faxis(mask_f)/1e3, Mag_rFM(mask_f), 'LineWidth', 1.0); hold on;
plot(faxis(mask_f)/1e3, Mag_rFMf(mask_f), 'LineWidth', 1.0);
grid on; xlabel('Frequency (kHz)'); ylabel('Magnitude');
legend('RFM(f) unfiltered', 'RFM(f) filtered', 'Location','best'); 
title('FM Spectrum: received vs filtered');

figure('Name','4.3 Spectrum: PM noisy vs filtered');
plot(faxis(mask_f)/1e3, Mag_rPM(mask_f), 'LineWidth', 1.0); hold on;
plot(faxis(mask_f)/1e3, Mag_rPMf(mask_f), 'LineWidth', 1.0);
grid on; xlabel('Frequency (kHz)'); ylabel('Magnitude');
legend('RPM(f) unfiltered', 'RPM(f) filtered', 'Location','best'); 
title('PM Spectrum: received vs filtered');

% 4.4) Output SNR comparison (before and after AWGN + loss, and after filtering)
PsFM = mean(sFM.^2);
PsPM = mean(sPM.^2);

Pn = mean(nFM.^2);

SNRin_FM_dB = 10*log10(PsFM / Pn);
SNRin_PM_dB = 10*log10(PsPM / Pn);

Psout_FM = mean((a*sFM).^2);
Psout_PM = mean((a*sPM).^2);
SNRout_FM_dB = 10*log10(Psout_FM / Pn);
SNRout_PM_dB = 10*log10(Psout_PM / Pn);

nFM_f = filtfilt(bpFM, 1, nFM);
nPM_f = filtfilt(bpPM, 1, nPM);

sFM_loss_f = filtfilt(bpFM, 1, a*sFM);
sPM_loss_f = filtfilt(bpPM, 1, a*sPM);

Psout_FM_f = mean(sFM_loss_f.^2);
Psout_PM_f = mean(sPM_loss_f.^2);
Pn_FM_f = mean(nFM_f.^2);
Pn_PM_f = mean(nPM_f.^2);

SNRout_FM_f_dB = 10*log10(Psout_FM_f / Pn_FM_f);
SNRout_PM_f_dB = 10*log10(Psout_PM_f / Pn_PM_f);

fprintf('\n================ SNR RESULTS ================\n');
fprintf('Noise PSD: N0/2 = %.3e W/Hz, Fs = %.3e Hz => sigma^2 = %.3e\n', N0_over2, Fs, sigma2);
fprintf('Loss: %d dB => amplitude a = %.3e, power scale a^2 = %.3e\n\n', Loss_dB, a, a^2);

fprintf('Reference SNR (no loss, with same AWGN):\n');
fprintf('  SNRin_FM = %.2f dB\n', SNRin_FM_dB);
fprintf('  SNRin_PM = %.2f dB\n\n', SNRin_PM_dB);

fprintf('After loss + AWGN (unfiltered):\n');
fprintf('  SNRout_FM = %.2f dB\n', SNRout_FM_dB);
fprintf('  SNRout_PM = %.2f dB\n\n', SNRout_PM_dB);

fprintf('After loss + AWGN + bandpass filtering:\n');
fprintf('  SNRout_FM_filtered = %.2f dB\n', SNRout_FM_f_dB);
fprintf('  SNRout_PM_filtered = %.2f dB\n', SNRout_PM_f_dB);
fprintf('============================================\n');

function [f, Mag] = mag_spectrum(x, Fs)
    x = x(:);
    N = length(x);
    dt = 1/Fs;
    X = fftshift(fft(x)) * dt;
    f = ((0:N-1)' - N/2) * (Fs/N);
    Mag = abs(X);
end



%% Part 5 – FM Demodulation and Error Comparison

Fs = 5e6;
dt = 1/Fs;
N  = 2^17;
t  = ((0:N-1) - N/2) * dt;  t = t(:);

m = 400 * tri_fun(400*t);               
m = m(:);

fc_plot = 200e3;
fm_max = 400;
betaFM_target = 10;
dfFM = min(max(betaFM_target * fm_max, 20e3), 40e3);

x = m / max(abs(m));
phiFM = 2*pi*dfFM * cumsum(x) * dt;
sFM = cos(2*pi*fc_plot*t + phiFM);

N0_over2 = 0.5e-12;
sigma = sqrt(N0_over2 * Fs);
Loss_dB = 50; a = 10^(-Loss_dB/20);
rFM = a*sFM + sigma*randn(N,1);

BC_FM = 2*(dfFM + fm_max);
Bhalf_FM = 0.6 * BC_FM;
ord = 600;
bpFM = fir1(ord, [(fc_plot - Bhalf_FM)/(Fs/2), (fc_plot + Bhalf_FM)/(Fs/2)], 'bandpass', hamming(ord+1));
rFM_tilde = filtfilt(bpFM, 1, rFM);

% 5.1) FM demodulation using fmdemod
r_in = rFM_tilde / max(abs(rFM_tilde) + 1e-12);
mhat1 = fmdemod(r_in, fc_plot, Fs, dfFM);
mhat1 = align_gain(mhat1, m);

figure('Name','5.1 fmdemod output vs message');
plot(t(mask_t)*1e3, m(mask_t), 'LineWidth', 1.2); hold on;
plot(t(mask_t)*1e3, mhat1(mask_t), 'LineWidth', 1.2);
grid on; xlabel('Time (ms)'); ylabel('Amplitude');
legend('m(t)', 'mhat1 (fmdemod)', 'Location','best'); 
title('FM Demodulation using fmdemod');

% 5.2) FM demodulation using a discriminator:
z = hilbert(r_in);               
dz = [0; diff(z)] * Fs;           
env = abs(dz);

env_dc_removed = env - mean(env);

fcut = 2.0e3;                     
lpOrd = 400;
bLP = fir1(lpOrd, fcut/(Fs/2), 'low', hamming(lpOrd+1), 'scale');

mhat2 = filtfilt(bLP, 1, env_dc_removed);

mhat2 = align_gain(mhat2, m);

figure('Name','5.2 Discriminator output vs message');
plot(t(mask_t)*1e3, m(mask_t), 'LineWidth', 1.2); hold on;
plot(t(mask_t)*1e3, mhat2(mask_t), 'LineWidth', 1.2);
grid on; xlabel('Time (ms)'); ylabel('Amplitude');
legend('m(t)', 'mhat2 (discriminator)', 'Location','best'); 
title('FM Demodulation using Differentiator + Envelope + LPF');

E1 = sum((m - mhat1).^2) / sum(m.^2);
E2 = sum((m - mhat2).^2) / sum(m.^2);

fprintf('\n--- Part 5 NMSE Results ---\n');
fprintf('E_norm for mhat1 (fmdemod):       %.6e\n', E1);
fprintf('E_norm for mhat2 (discriminator): %.6e\n', E2);

figure('Name','5.3 Error signals');
plot(t(mask_t)*1e3, (m(mask_t)-mhat1(mask_t)), 'LineWidth', 1.1); hold on;
plot(t(mask_t)*1e3, (m(mask_t)-mhat2(mask_t)), 'LineWidth', 1.1);
grid on; xlabel('Time (ms)'); ylabel('Error');
legend('m - mhat1', 'm - mhat2', 'Location','best'); 
title('Demodulation Error Signals');

function y_aligned = align_gain(y, ref)
    y0 = y - mean(y);
    r0 = ref - mean(ref);
    g = (r0' * y0) / (y0' * y0 + 1e-12);
    y_aligned = g * y0;
end


%% Part 6 – Audio FDM

doRecord = true;           
recSeconds = 3.5;          
Fs_rec = 16000;            

FsFDM = 200e3;            

DeltaF_default = 8e3;      

guard = 3e3;              

N0_over2 = 0.5e-12;        
useChannelLoss = false;    
Loss_dB = 50;              

plotSpanKHz = 60;         

% 6.1) Record three voice signals
if doRecord
    fprintf('--- 6.1 Recording ---\n');
    fprintf('Record 1: say your FULL NAME.\n');
    m1 = record_voice('FULL NAME', Fs_rec, recSeconds);

    fprintf('Record 2: say your STUDENT ID.\n');
    m2 = record_voice('STUDENT ID', Fs_rec, recSeconds);

    fprintf('Record 3: say your FIELD OF STUDY + COURSE NAME.\n');
    m3 = record_voice('FIELD + COURSE', Fs_rec, recSeconds);
else
    [m1, Fs1] = audioread('voice1_name.wav');  m1 = mean(m1,2);
    [m2, Fs2] = audioread('voice2_id.wav');    m2 = mean(m2,2);
    [m3, Fs3] = audioread('voice3_field.wav'); m3 = mean(m3,2);
    assert(Fs1==Fs2 && Fs2==Fs3, 'All WAV files must have the same Fs.');
    Fs_rec = Fs1;
end

m1 = preprocess_audio(m1);
m2 = preprocess_audio(m2);
m3 = preprocess_audio(m3);

Lmin = min([length(m1), length(m2), length(m3)]);
m1 = m1(1:Lmin); m2 = m2(1:Lmin); m3 = m3(1:Lmin);

t_rec = (0:Lmin-1)'/Fs_rec;

% 6.2) Effective bandwidth Beff,i (90% energy)
[Beff1, f1, M1mag] = eff_bw_90_spectrum(m1, Fs_rec);
[Beff2, f2, M2mag] = eff_bw_90_spectrum(m2, Fs_rec);
[Beff3, f3, M3mag] = eff_bw_90_spectrum(m3, Fs_rec);

T = table((1:3)', [Beff1; Beff2; Beff3], [Beff1; Beff2; Beff3]/1e3, ...
    'VariableNames', {'Channel','Beff_Hz','Beff_kHz'});
disp('--- 6.2 Effective bandwidth (90% energy) ---');
disp(T);

figure('Name','6.2 Voice spectra (baseband)');
subplot(3,1,1); plot(f1/1e3, M1mag); grid on; xline(Beff1/1e3,'--'); xline(-Beff1/1e3,'--');
title('Voice 1 spectrum |M_1(f)|'); xlabel('kHz'); ylabel('Mag');
subplot(3,1,2); plot(f2/1e3, M2mag); grid on; xline(Beff2/1e3,'--'); xline(-Beff2/1e3,'--');
title('Voice 2 spectrum |M_2(f)|'); xlabel('kHz'); ylabel('Mag');
subplot(3,1,3); plot(f3/1e3, M3mag); grid on; xline(Beff3/1e3,'--'); xline(-Beff3/1e3,'--');
title('Voice 3 spectrum |M_3(f)|'); xlabel('kHz'); ylabel('Mag');

% 6.3) Bandwidth limiting via lowpass filtering
m1LP = lowpass_fir_zero_phase(m1, Fs_rec, 1.05*Beff1);
m2LP = lowpass_fir_zero_phase(m2, Fs_rec, 1.05*Beff2);
m3LP = lowpass_fir_zero_phase(m3, Fs_rec, 1.05*Beff3);

figure('Name','6.3 Time comparison (original vs LP)');
mask = (t_rec <= 0.08);
subplot(3,1,1); plot(t_rec(mask), m1(mask)); hold on; plot(t_rec(mask), m1LP(mask));
grid on; title('Voice 1'); legend('m_1','m_{1,LP}');
subplot(3,1,2); plot(t_rec(mask), m2(mask)); hold on; plot(t_rec(mask), m2LP(mask));
grid on; title('Voice 2'); legend('m_2','m_{2,LP}');
subplot(3,1,3); plot(t_rec(mask), m3(mask)); hold on; plot(t_rec(mask), m3LP(mask));
grid on; title('Voice 3'); legend('m_3','m_{3,LP}'); xlabel('Time (s)');

m1r = resample(m1LP, FsFDM, Fs_rec);
m2r = resample(m2LP, FsFDM, Fs_rec);
m3r = resample(m3LP, FsFDM, Fs_rec);

LminF = min([length(m1r), length(m2r), length(m3r)]);
m1r = m1r(1:LminF); m2r = m2r(1:LminF); m3r = m3r(1:LminF);
tF = (0:LminF-1)'/FsFDM;

fm1 = Beff1; fm2 = Beff2; fm3 = Beff3;

% 6.4) FM modulation and audio FDM construction
Df1 = DeltaF_default;
Df2 = DeltaF_default;
Df3 = DeltaF_default;

BW1 = 2*(Df1 + fm1);
BW2 = 2*(Df2 + fm2);
BW3 = 2*(Df3 + fm3);

fc1 = 25e3;
fc2 = fc1 + (BW1/2) + guard + (BW2/2);
fc3 = fc2 + (BW2/2) + guard + (BW3/2);

nyq = FsFDM/2;
margin = 5e3;
if (fc3 + BW3/2) > (nyq - margin)
    scale = (nyq - margin - fc3) / (BW3/2);
    scale = max(scale, 0.2);
    Df1 = Df1*scale; Df2 = Df2*scale; Df3 = Df3*scale;
    BW1 = 2*(Df1 + fm1); BW2 = 2*(Df2 + fm2); BW3 = 2*(Df3 + fm3);
    fc2 = fc1 + (BW1/2) + guard + (BW2/2);
    fc3 = fc2 + (BW2/2) + guard + (BW3/2);
end

fprintf('\n--- 6.4 Carrier plan (Hz) ---\n');
fprintf('FsFDM = %.1f kHz, Nyquist = %.1f kHz\n', FsFDM/1e3, nyq/1e3);
fprintf('Ch1: fc=%.1f kHz, Δf=%.1f kHz, CarsonBW≈%.1f kHz\n', fc1/1e3, Df1/1e3, BW1/1e3);
fprintf('Ch2: fc=%.1f kHz, Δf=%.1f kHz, CarsonBW≈%.1f kHz\n', fc2/1e3, Df2/1e3, BW2/1e3);
fprintf('Ch3: fc=%.1f kHz, Δf=%.1f kHz, CarsonBW≈%.1f kHz\n', fc3/1e3, Df3/1e3, BW3/1e3);

s1 = fm_mod_real(m1r, FsFDM, fc1, Df1);
s2 = fm_mod_real(m2r, FsFDM, fc2, Df2);
s3 = fm_mod_real(m3r, FsFDM, fc3, Df3);

sFDM = s1 + s2 + s3;

% 6.5) AWGN channel for FDM signal
if useChannelLoss
    a = 10^(-Loss_dB/20);
else
    a = 1.0;
end

sigma2 = N0_over2 * FsFDM;
sigma = sqrt(sigma2);
n = sigma * randn(length(sFDM),1);

r = a*sFDM + n;

figure('Name','6.5 FDM waveform');
mask = (tF <= 1.0e-3);
plot(tF(mask)*1e3, sFDM(mask)); hold on;
plot(tF(mask)*1e3, r(mask));
grid on; xlabel('Time (ms)'); ylabel('Amplitude');
legend('s_{FDM}(t)','r(t)','Location','best');
title('FDM signal before/after AWGN');

[fF, MagSFDM] = mag_spectrum(sFDM, FsFDM);
[~,  MagR]    = mag_spectrum(r, FsFDM);

figure('Name','6.5 FDM spectrum');
plot(fF/1e3, MagSFDM, 'LineWidth', 1.0); hold on;
plot(fF/1e3, MagR, 'LineWidth', 1.0);
grid on; xlabel('Frequency (kHz)'); ylabel('Magnitude');
legend('|S_{FDM}(f)|','|R(f)|','Location','best');
title('FDM spectrum (clean vs received)');

% 6.6) Separation of FM channels (bandpass filters)
Bhalf1 = 0.6*BW1; 
Bhalf2 = 0.6*BW2;
Bhalf3 = 0.6*BW3;

ordBP = 900;      

bp1 = fir1(ordBP, [(fc1-Bhalf1)/(FsFDM/2), (fc1+Bhalf1)/(FsFDM/2)], 'bandpass', hamming(ordBP+1));
bp2 = fir1(ordBP, [(fc2-Bhalf2)/(FsFDM/2), (fc2+Bhalf2)/(FsFDM/2)], 'bandpass', hamming(ordBP+1));
bp3 = fir1(ordBP, [(fc3-Bhalf3)/(FsFDM/2), (fc3+Bhalf3)/(FsFDM/2)], 'bandpass', hamming(ordBP+1));

s1_tilde = filtfilt(bp1, 1, r);
s2_tilde = filtfilt(bp2, 1, r);
s3_tilde = filtfilt(bp3, 1, r);

% 6.7) FM demodulation of separated channels
m1hatF = fm_demod_auto(s1_tilde, FsFDM, fc1, Df1);
m2hatF = fm_demod_auto(s2_tilde, FsFDM, fc2, Df2);
m3hatF = fm_demod_auto(s3_tilde, FsFDM, fc3, Df3);

m1hatF = lowpass_fir_zero_phase(m1hatF, FsFDM, 1.05*Beff1);
m2hatF = lowpass_fir_zero_phase(m2hatF, FsFDM, 1.05*Beff2);
m3hatF = lowpass_fir_zero_phase(m3hatF, FsFDM, 1.05*Beff3);

m1hat = resample(m1hatF, Fs_rec, FsFDM);
m2hat = resample(m2hatF, Fs_rec, FsFDM);
m3hat = resample(m3hatF, Fs_rec, FsFDM);

LminR = min([length(m1LP), length(m1hat)]);
m1LPc = m1LP(1:LminR); m1hat = m1hat(1:LminR);

LminR = min([length(m2LP), length(m2hat)]);
m2LPc = m2LP(1:LminR); m2hat = m2hat(1:LminR);

LminR = min([length(m3LP), length(m3hat)]);
m3LPc = m3LP(1:LminR); m3hat = m3hat(1:LminR);

[m1hatA, ~] = align_delay_gain(m1hat, m1LPc);
[m2hatA, ~] = align_delay_gain(m2hat, m2LPc);
[m3hatA, ~] = align_delay_gain(m3hat, m3LPc);

t_show = (0:length(m1LPc)-1)'/Fs_rec;
mask = (t_show <= 0.08);

figure('Name','6.7 Recovered audio (time-domain)');
subplot(3,1,1); plot(t_show(mask), m1LPc(mask)); hold on; plot(t_show(mask), m1hatA(mask));
grid on; title('Channel 1'); legend('m_{1,LP}','\hat{m}_1');
subplot(3,1,2); plot(t_show(mask), m2LPc(mask)); hold on; plot(t_show(mask), m2hatA(mask));
grid on; title('Channel 2'); legend('m_{2,LP}','\hat{m}_2');
subplot(3,1,3); plot(t_show(mask), m3LPc(mask)); hold on; plot(t_show(mask), m3hatA(mask));
grid on; title('Channel 3'); legend('m_{3,LP}','\hat{m}_3'); xlabel('Time (s)');

% 6.8) Normalized error metrics and discussion helpers
E1 = sum((m1LPc - m1hatA).^2) / (sum(m1LPc.^2) + 1e-12);
E2 = sum((m2LPc - m2hatA).^2) / (sum(m2LPc.^2) + 1e-12);
E3 = sum((m3LPc - m3hatA).^2) / (sum(m3LPc.^2) + 1e-12);

TE = table((1:3)', [E1;E2;E3], 'VariableNames', {'Channel','Enorm'});
disp('--- 6.8 Normalized errors ---');
disp(TE);

% soundsc(m1LPc, Fs_rec); pause(length(m1LPc)/Fs_rec + 0.5);
% soundsc(m1hatA, Fs_rec);

function x = record_voice(label, Fs, seconds)
    recObj = audiorecorder(Fs, 16, 1);
    input(sprintf('Press ENTER to start recording: %s (%.1f s)\n', label, seconds), 's');
    recordblocking(recObj, seconds);
    x = getaudiodata(recObj);
    fprintf('Recorded "%s"\n\n', label);
end

function y = preprocess_audio(x)
    x = x(:);
    x = x - mean(x);
    mx = max(abs(x)) + 1e-12;
    y = x / mx;
end

function [Beff, f, Mag] = eff_bw_90_spectrum(x, Fs)
    x = x(:);
    N = 2^nextpow2(length(x));
    xz = [x; zeros(N-length(x),1)];
    dt = 1/Fs;

    X  = fftshift(fft(xz)) * dt;
    f  = ((0:N-1)' - N/2) * (Fs/N);
    Mag = abs(X);

    PSD = abs(X).^2;
    df = Fs/N;

    pos = (f >= 0);
    fpos = f(pos);
    Ppos = PSD(pos);

    Etotal = (Ppos(1) + 2*sum(Ppos(2:end))) * df;
    Ecum = Ppos(1)*df + 2*cumsum(Ppos(2:end))*df;

    ratio = Ecum / (Etotal + 1e-18);
    idx = find(ratio >= 0.90, 1, 'first');
    Beff = fpos(idx);  
end

function y = lowpass_fir_zero_phase(x, Fs, fcut)
    x = x(:);
    fcut = min(max(fcut, 10), 0.49*(Fs/2)); 
    ord = 400;
    b = fir1(ord, fcut/(Fs/2), 'low', hamming(ord+1), 'scale');
    y = filtfilt(b, 1, x);
end

function s = fm_mod_real(m, Fs, fc, dF)
    m = m(:);
    mx = max(abs(m)) + 1e-12;
    x = m / mx;

    dt = 1/Fs;
    phi = 2*pi*dF * cumsum(x) * dt;
    t = (0:length(m)-1)'/Fs;
    s = cos(2*pi*fc*t + phi);
end

function mhat = fm_demod_auto(s, Fs, fc, dF)
    s = s(:);
    s = s / (max(abs(s)) + 1e-12);

    if exist('fmdemod', 'file') == 2
        mhat = fmdemod(s, fc, Fs, dF);
    else
        z = hilbert(s);
        ph = unwrap(angle(z));
        inst_f = [0; diff(ph)] * Fs/(2*pi);    
        mhat = (inst_f - fc) / (dF + 1e-12);   
        mhat = mhat - mean(mhat);
    end
end

function [y_aligned, lag] = align_delay_gain(y, ref)
    y = y(:); ref = ref(:);

    maxLag = min(round(0.02 * length(ref)), 2000);
    [c, lags] = xcorr(ref, y, maxLag, 'coeff');
    [~, idx] = max(c);
    lag = lags(idx);

    if lag > 0
        yshift = [zeros(lag,1); y(1:end-lag)];
    elseif lag < 0
        lag2 = -lag;
        yshift = [y(lag2+1:end); zeros(lag2,1)];
    else
        yshift = y;
    end

    y0 = yshift - mean(yshift);
    r0 = ref - mean(ref);
    g = (r0' * y0) / (y0' * y0 + 1e-12);

    y_aligned = g * y0;
end

