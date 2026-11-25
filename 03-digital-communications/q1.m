%% q1
Rb = 1e3;           % Bit rate 
Tb = 1/Rb;          % Bit duration
fs = 20*Rb;         % Sampling frequency 
Ts = 1/fs;          % Sampling period
L = 20;             % Samples per symbol
N_total_bits = 1e5; % Total bits for spectrum estimation
num_bits_plot = 10; % Number of bits to plot for waveform visualization

% Generate a large set of random binary data for accurate PSD estimation
data_bits = randi([0 1], 1, N_total_bits);
data_symbols = 2 * data_bits - 1; % Map to BPSK symbols (+1 or -1)

g_rect = 1; 
n_pulse = 0:L-1;
g_half_sine = sin(pi * n_pulse / L); 

total_samples = N_total_bits * L;

% Rectangular Signal (NRZ)
rect_upsampled_L = repelem(data_symbols, L);
rect_signal_t = rect_upsampled_L;

% Half-Sine Signal
half_sine_signal_t = zeros(1, total_samples);
for k = 1:N_total_bits
    start = (k-1)*L + 1;
    end_idx = k*L;
    half_sine_signal_t(start:end_idx) = data_symbols(k) * g_half_sine;
end

%1
figure(1);
t_plot = 0:(num_bits_plot*L - 1);
t_plot = t_plot * Ts; % Time vector in seconds

subplot(2,1,1);
plot(t_plot, rect_signal_t(1:(num_bits_plot*L)), 'b', 'LineWidth', 1.5);
title('1. Rectangular (NRZ) BPSK Signal (10 bits)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
xlim([0 10*Tb]);

subplot(2,1,2);
plot(t_plot, half_sine_signal_t(1:(num_bits_plot*L)), 'r', 'LineWidth', 1.5);
title('1. Half-Sine BPSK Signal (10 bits)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
xlim([0 10*Tb]);
sgtitle('Time-Domain Waveforms');

% 2
NFFT = 4096; % FFT size
overlap = floor(NFFT/2); % 50% overlap
window = hamming(NFFT);

% Calculate PSD for Rectangular (NRZ)
[Pxx_rect, F_rect] = pwelch(rect_signal_t, window, overlap, NFFT, fs, 'psd');

% Calculate PSD for Half-Sine
[Pxx_half, F_half] = pwelch(half_sine_signal_t, window, overlap, NFFT, fs, 'psd');

figure(2);
Pxx_rect_dB = 10*log10(Pxx_rect);
Pxx_half_dB = 10*log10(Pxx_half);

plot_limit_idx = find(F_rect > 4*Rb, 1, 'first');

subplot(2,1,1);
plot(F_rect(1:plot_limit_idx), Pxx_rect_dB(1:plot_limit_idx), 'b', 'LineWidth', 1.5);
title('PSD of Rectangular Pulse');
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
grid on;
xlim([0 4*Rb]);

subplot(2,1,2);
plot(F_half(1:plot_limit_idx), Pxx_half_dB(1:plot_limit_idx), 'r', 'LineWidth', 1.5);
title('PSD of Half-Sine Pulse');
xlabel('Frequency (Hz)');
ylabel('PSD (dB/Hz)');
grid on;
xlim([0 4*Rb]);
sgtitle('Power Spectral Density Estimates');

% 3
figure(3);

plot(F_rect(1:plot_limit_idx), Pxx_rect_dB(1:plot_limit_idx), 'b', 'LineWidth', 1.5);
hold on;
plot(F_half(1:plot_limit_idx), Pxx_half_dB(1:plot_limit_idx), 'r', 'LineWidth', 1.5);
hold off;

title('Comparative Spectrum Analysis');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('Rectangular (NRZ)', 'Half-Sine', 'Location', 'northeast');
grid on;
xlim([0 4*Rb]);

%% q2
% Eb/N0 values to test 
EbN0_dB_test = [100 0];
EbN0_lin_test = 10.^(EbN0_dB_test/10);

num_samples_plot = num_bits_plot * L;

% Preallocate noise-corrupted signals
rx_rect = cell(1, length(EbN0_dB_test));
rx_half = cell(1, length(EbN0_dB_test));

for i = 1:length(EbN0_lin_test)
    
    EbN0 = EbN0_lin_test(i);
    
    % Discrete-time noise variance
    sigma2 = 1 / (2 * L * EbN0);
    sigma = sqrt(sigma2);
    
    % Generate AWGN
    noise_rect = sigma * randn(1, total_samples);
    noise_half = sigma * randn(1, total_samples);
    
    % Add noise
    rx_rect{i} = rect_signal_t + noise_rect;
    rx_half{i} = half_sine_signal_t + noise_half;
end
figure(4);

% Eb/N0 = 100 dB 
subplot(2,2,1);
plot(t_plot, rect_signal_t(1:num_samples_plot), 'k', 'LineWidth', 1.2);
hold on;
plot(t_plot, rx_rect{1}(1:num_samples_plot), 'r--');
hold off;
title('Rectangular BPSK: Eb/N0 = 100 dB');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Noise-Free', 'With Noise');
grid on;

subplot(2,2,2);
plot(t_plot, half_sine_signal_t(1:num_samples_plot), 'k', 'LineWidth', 1.2);
hold on;
plot(t_plot, rx_half{1}(1:num_samples_plot), 'r--');
hold off;
title('Half-Sine BPSK: Eb/N0 = 100 dB');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Noise-Free', 'With Noise');
grid on;

% Eb/N0 = 0 dB 
subplot(2,2,3);
plot(t_plot, rect_signal_t(1:num_samples_plot), 'k', 'LineWidth', 1.2);
hold on;
plot(t_plot, rx_rect{2}(1:num_samples_plot), 'r--');
hold off;
title('Rectangular BPSK: Eb/N0 = 0 dB');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Noise-Free', 'With Noise');
grid on;

subplot(2,2,4);
plot(t_plot, half_sine_signal_t(1:num_samples_plot), 'k', 'LineWidth', 1.2);
hold on;
plot(t_plot, rx_half{2}(1:num_samples_plot), 'r--');
hold off;
title('Half-Sine BPSK: Eb/N0 = 0 dB');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Noise-Free', 'With Noise');
grid on;

sgtitle(' Noise Sensitivity Test');

%% 3

% Matched filter for half-sine pulse
h_mf = fliplr(g_half_sine);

EbN0_vis_dB = 10;
EbN0_vis = 10^(EbN0_vis_dB/10);

% Noise variance
sigma2_vis = 1 / (2 * L * EbN0_vis);
sigma_vis = sqrt(sigma2_vis);

noise_vis = sigma_vis * randn(1, total_samples);
rx_half_vis = half_sine_signal_t + noise_vis;

mf_output = conv(rx_half_vis, h_mf, 'same');

% Optimal sampling instants 
sample_indices = L:L:(num_bits_plot*L);
sample_times = sample_indices * Ts;
sample_values = mf_output(sample_indices);

% Plot MF output
figure(5);
t_mf = (0:(num_bits_plot*L-1)) * Ts;

plot(t_mf, mf_output(1:num_bits_plot*L), 'b', 'LineWidth', 1.3);
hold on;
stem(sample_times, sample_values, 'r', 'filled');
hold off;

xlabel('Time (s)');
ylabel('Amplitude');
title('Matched Filter Output (Half-Sine Pulse)');
legend('MF Output', 'Optimal Sampling Points');
grid on;

%3.2
EbN0_dB = 0:2:16;
EbN0_lin = 10.^(EbN0_dB/10);
SQNR_dB = zeros(size(EbN0_dB));

% Noise-free matched filter output
mf_clean = conv(half_sine_signal_t, h_mf, 'same');
clean_samples = mf_clean(L:L:end);

for i = 1:length(EbN0_lin)
   
    sigma2 = 1 / (2 * L * EbN0_lin(i));
    sigma = sqrt(sigma2);
    
    noise = sigma * randn(1, total_samples);
    rx_noisy = half_sine_signal_t + noise;
    
    mf_noisy = conv(rx_noisy, h_mf, 'same');
    noisy_samples = mf_noisy(L:L:end);
 
    signal_power = mean(clean_samples.^2);
    noise_power = mean((noisy_samples - clean_samples).^2);
    
    SQNR_dB(i) = 10*log10(signal_power / noise_power);
end

%SQNR curve
figure(6);
plot(EbN0_dB, SQNR_dB, 'o-', 'LineWidth', 1.5);
xlabel('E_b/N_0 (dB)');
ylabel('SQNR (dB)');
title('SQNR vs E_b/N_0');
grid on;

% 3.3
p = polyfit(EbN0_dB, SQNR_dB, 1);
SQNR_fit = polyval(p, EbN0_dB);

% Plot fitted line
figure(7);
plot(EbN0_dB, SQNR_dB, 'bo', 'LineWidth', 1.5);
hold on;
plot(EbN0_dB, SQNR_fit, 'r--', 'LineWidth', 1.5);
hold off;

xlabel('E_b/N_0 (dB)');
ylabel('SQNR (dB)');
title(['Linearity of SQNR vs E_b/N_0, Slope = ', num2str(p(1), '%.2f'), ' dB/dB']);
legend('Measured SQNR', 'Linear Fit');
grid on;

%% 4
EbN0_dB = 0:1:20;                 
EbN0_lin = 10.^(EbN0_dB/10);

Ne_target = 200;                  % Minimum bit errors
Nmax = 1e6;                       % Maximum transmitted bits

BER_BPSK  = zeros(size(EbN0_dB));
BER_QPSK  = zeros(size(EbN0_dB));
BER_QAM16 = zeros(size(EbN0_dB));

%  BPSK 
for i = 1:length(EbN0_dB)
    errors = 0; bits = 0;
    while errors < Ne_target && bits < Nmax
        Nb = 1e4;
        b = randi([0 1], Nb, 1);
        s = 2*b - 1;                           % BPSK symbols (Es = 1)

        sigma = sqrt(1/(2*EbN0_lin(i)));
        r = s + sigma*randn(size(s));

        b_hat = r > 0;
        errors = errors + sum(b ~= b_hat);
        bits = bits + Nb;
    end
    BER_BPSK(i) = errors / bits;
end

%  QPSK 
for i = 1:length(EbN0_dB)
    errors = 0; bits = 0;
    while errors < Ne_target && bits < Nmax
        Nb = 2e4;
        b = randi([0 1], Nb, 1);

        bI = b(1:2:end);
        bQ = b(2:2:end);

        s = (2*bI-1 + 1j*(2*bQ-1))/sqrt(2);    % Unit energy QPSK

        k = 2;                                % bits/symbol
        sigma = sqrt(1/(2*k*EbN0_lin(i)));
        r = s + sigma*(randn(size(s)) + 1j*randn(size(s)));

        bI_hat = real(r) > 0;
        bQ_hat = imag(r) > 0;
        b_hat = reshape([bI_hat.'; bQ_hat.'], [], 1);

        errors = errors + sum(b ~= b_hat);
        bits = bits + Nb;
    end
    BER_QPSK(i) = errors / bits;
end

% 16-QAM 
M = 16;
k = log2(M);

for i = 1:length(EbN0_dB)
    errors = 0; bits = 0;
    while errors < Ne_target && bits < Nmax
        Nb = 4e4;
        b = randi([0 1], Nb, 1);

        symbols = bi2de(reshape(b, [], k), 'left-msb');
        s = qammod(symbols, M, 'gray', 'UnitAveragePower', true);

        sigma = sqrt(1/(2*k*EbN0_lin(i)));
        r = s + sigma*(randn(size(s)) + 1j*randn(size(s)));

        symbols_hat = qamdemod(r, M, 'gray', 'UnitAveragePower', true);
        b_hat = reshape(de2bi(symbols_hat, k, 'left-msb').', [], 1);

        errors = errors + sum(b ~= b_hat);
        bits = bits + Nb;
    end
    BER_QAM16(i) = errors / bits;
end

%  Theoretical BER 
BER_BPSK_th  = qfunc(sqrt(2*EbN0_lin));
BER_QPSK_th  = BER_BPSK_th;
BER_QAM16_th = (3/8)*erfc(sqrt((4/10)*EbN0_lin));

%  Plot BER Curves 
figure;
semilogy(EbN0_dB, BER_BPSK,  'bo-', 'LineWidth', 1.5); hold on;
semilogy(EbN0_dB, BER_QPSK,  'rs-', 'LineWidth', 1.5);
semilogy(EbN0_dB, BER_QAM16, 'kd-', 'LineWidth', 1.5);

semilogy(EbN0_dB, BER_BPSK_th,  'b--', 'LineWidth', 1.5);
semilogy(EbN0_dB, BER_QPSK_th,  'r--', 'LineWidth', 1.5);
semilogy(EbN0_dB, BER_QAM16_th, 'k--', 'LineWidth', 1.5);

grid on;
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER Performance Comparison');
legend('BPSK (Sim)', 'QPSK (Sim)', '16-QAM (Sim)', ...
       'BPSK (Theory)', 'QPSK (Theory)', '16-QAM (Theory)', ...
       'Location', 'southwest');

% MSE Computation 
MSE_BPSK  = mean((BER_BPSK  - BER_BPSK_th ).^2);
MSE_QPSK  = mean((BER_QPSK  - BER_QPSK_th ).^2);
MSE_QAM16 = mean((BER_QAM16 - BER_QAM16_th).^2);

fprintf('MSE BPSK  = %.3e\n', MSE_BPSK);
fprintf('MSE QPSK  = %.3e\n', MSE_QPSK);
fprintf('MSE 16QAM = %.3e\n', MSE_QAM16);

%  Power Penalty @ BER = 1e-3 
targetBER = 1e-3;

Eb_BPSK  = EbN0_dB(find(BER_BPSK  <= targetBER, 1, 'first'));
Eb_QPSK  = EbN0_dB(find(BER_QPSK  <= targetBER, 1, 'first'));
%Eb_QAM16 = EbN0_dB(find(BER_QAM16 <= targetBER, 1, 'first'));
Eb_QAM16 = interp1(BER_QAM16_th, EbN0_dB, targetBER);


power_penalty = Eb_QAM16 - Eb_BPSK;

fprintf('\nEb/N0 required for BER = 10^-3:\n');
fprintf('BPSK   = %.2f dB\n', Eb_BPSK);
fprintf('QPSK   = %.2f dB\n', Eb_QPSK);
fprintf('16-QAM = %.2f dB\n', Eb_QAM16);
fprintf('Power penalty (16-QAM vs BPSK) = %.2f dB\n', power_penalty);

%% 5

Rb = 1e3;                
Tb = 1/Rb;                 
L  = 20;                  
fs = L*Rb;                  
Ts = 1/fs;
Nb = 2e4;                   

n = 0:L-1;
p = sin(pi*(n+1)/L);
p = p / norm(p);           
h = fliplr(p);              

b = randi([0 1], Nb, 1);
bI = b(1:2:end);
bQ = b(2:2:end);
symbols = (2*bI-1 + 1j*(2*bQ-1))/sqrt(2);

x_up = upsample(symbols, L);
tx = conv(x_up, p);

EbN0_dB_phase = 10;
EbN0 = 10^(EbN0_dB_phase/10);
sigma = sqrt(1/(2*EbN0));

rx = tx + sigma*(randn(size(tx)) + 1j*randn(size(tx)));

y = conv(rx, h);
delay = L;
samples = y(delay:L:end);

% 5.1 
phi = 15*pi/180;  % Phase rotation
%  Before MF sampling 
rx_phase_beforeMF = rx * exp(1j*phi);  

samples_beforeMF = rx_phase_beforeMF(L:L:end);

% After MF sampling 
samples_phase_afterMF = samples * exp(1j*phi);  % MF output + phase rotation

figure;
subplot(1,2,1)
scatter(real(samples_beforeMF), imag(samples_beforeMF), 'b.');
grid on; axis square;
xlabel('In-Phase'); ylabel('Quadrature');
title('QPSK Constellation BEFORE MF (15° Phase)');

subplot(1,2,2)
scatter(real(samples_phase_afterMF), imag(samples_phase_afterMF), 'r.');
grid on; axis square;
xlabel('In-Phase'); ylabel('Quadrature');
title('QPSK Constellation AFTER MF (15° Phase)');


% 5.2
df = 0.01*Rb;                        % Frequency offset
t = (0:length(tx)-1).' * Ts;

rx_freq = tx .* exp(1j*2*pi*df*t) + sigma*(randn(size(tx)) + 1j*randn(size(tx)));

% Before MF
samples_freq_beforeMF = rx_freq(L:L:end);

% After MF
y_freq = conv(rx_freq, h);
samples_freq_afterMF = y_freq(delay:L:end);

figure;
subplot(1,2,1)
scatter(real(samples_freq_beforeMF), imag(samples_freq_beforeMF), 'b.');
grid on; axis square;
xlabel('In-Phase'); ylabel('Quadrature');
title('QPSK Constellation BEFORE MF (Freq Offset)');

subplot(1,2,2)
scatter(real(samples_freq_afterMF), imag(samples_freq_afterMF), 'r.');
grid on; axis square;
xlabel('In-Phase'); ylabel('Quadrature');
title('QPSK Constellation AFTER MF (Freq Offset)');

% 5.3 
EbN0_dB_eye = [5 20];

for k = 1:2
    EbN0 = 10^(EbN0_dB_eye(k)/10);
    sigma = sqrt(1/(2*EbN0));
    
    rx_eye = tx + sigma*(randn(size(tx)) + 1j*randn(size(tx)));
    y_eye = conv(rx_eye, h);
    
    % eyediagram creates its own figure automatically
    eyediagram(real(y_eye), 2*L);
    title(['Eye Diagram at E_b/N_0 = ' num2str(EbN0_dB_eye(k)) ' dB']);
end

%% 6
R = 1;      
C = 1;      
W = 5;      

f_pos = linspace(0, 20, 1000);
f = [-fliplr(f_pos), f_pos]; 

% High-Pass RC Filter 
H_RC_pos = (1j*2*pi*f_pos*R*C) ./ (1 + 1j*2*pi*f_pos*R*C);
H_RC = [fliplr(conj(H_RC_pos)), H_RC_pos];

% Ideal Low-Pass Filter
H_LPF = double(abs(f) <= W);

% Total Transfer Function
V_f = H_RC .* H_LPF;

% Magnitude
V_mag = abs(V_f);

% Plot |V(f)|
figure;
plot(f, V_mag, 'b', 'LineWidth', 2);
grid on;
xlabel('Frequency (Hz)');
ylabel('|V(f)|');
title('Two-Sided Total Transfer Function Magnitude |V(f)|');
xlim([min(f) max(f)]);
ylim([0 1.1]);

% Identify Effective Passband
passband_idx = find(V_mag >= 0.707);  % 3-dB criterion
f_passband = f([passband_idx(1), passband_idx(end)]);
hold on;
yline(0.707, 'r--', '3-dB Level');
text(f_passband(1), 0.75, sprintf('Passband ~ %.2f Hz to %.2f Hz', f_passband(1), f_passband(2)));

%% 7
Ac = 2;      
N0 = 0.1;     
fc = 1;       
R = 1; C = 1; 

% LPF Bandwidth Vector
Wvec = linspace(fc*1.001, 50, 2000); % W > fc

a = 2*pi*R*C;

% SNR = (A^2 / N0) * (pi*R*C*fc^2)/(1 + (2*pi*R*C*fc)^2) * (2*pi*R*C*W - atan(2*pi*R*C*W))
num = (pi*R*C*fc^2) / (1 + (a*fc)^2);
den = a*Wvec - atan(a*Wvec);
SNR_th = (Ac^2 / N0) * (num ./ den);

SNR_dB = 10*log10(SNR_th);

% Plot SNR vs W
figure('Name','Output SNR vs LPF Bandwidth W');
plot(Wvec, SNR_dB, 'b', 'LineWidth', 1.5); grid on;
xlabel('LPF Bandwidth W (Hz)');
ylabel('Output SNR (dB)');
title('Theoretical Output SNR vs LPF Bandwidth W');

% Find Optimal W (Maximum SNR)
[SNR_max, idxMax] = max(SNR_th);
W_opt = Wvec(idxMax);
SNR_opt_dB = 10*log10(SNR_max);

hold on;
plot(W_opt, SNR_opt_dB, 'ro', 'MarkerFaceColor','r','MarkerSize',6);
text(W_opt + 0.5, SNR_opt_dB, sprintf('W_{opt}=%.3f Hz', W_opt), 'VerticalAlignment','bottom');

fprintf('Optimal LPF bandwidth: W_opt = %.6f Hz\n', W_opt);
fprintf('Maximum SNR: %.4f dB\n', SNR_opt_dB);

%% 8 
Ac = 2;
N0 = 0.1;
R = 1;
C = 1;
fc = 1;                 
a = 2*pi*R*C;

% SNR for W = fc and W → ∞ 

W_fc = fc;
W_inf = 1e6;             % Large value approximating infinity

num = (pi*R*C*fc^2) / (1 + (a*fc)^2);

% SNR at W = fc
den_fc = (a*W_fc) - atan(a*W_fc);
SNR_fc = (Ac^2/N0) * (num / den_fc);
SNR_fc_dB = 10*log10(SNR_fc);

% SNR for very large W
den_inf = (a*W_inf) - atan(a*W_inf);
SNR_inf = (Ac^2/N0) * (num / den_inf);
SNR_inf_dB = 10*log10(SNR_inf);

fprintf('SNR at W = f_c        : %.3f dB\n', SNR_fc_dB);
fprintf('SNR for W → infinity : %.3f dB\n', SNR_inf_dB);

H_RC_fc = (2*pi*fc*R*C) / sqrt(1 + (2*pi*fc*R*C)^2);

% Signal power 
P_signal = (Ac^2/2) * H_RC_fc^2;

fprintf('Signal power at output: %.4f W\n', P_signal);

% 8.2 Output Noise PSD for W = 10 Hz

W = 10;
f = linspace(-20, 20, 10001);

% RC high-pass filter magnitude squared
H_RC_sq = ( (2*pi*f*R*C).^2 ) ./ (1 + (2*pi*f*R*C).^2 );

% Ideal LPF
H_LPF = double(abs(f) <= W);

% Total magnitude squared
V_sq = H_RC_sq .* H_LPF;

% Output noise PSD
S_out = (N0/2) * V_sq;

% Plot
figure;
plot(f, 10*log10(S_out + eps), 'LineWidth', 1.5);
grid on;
xlabel('Frequency (Hz)');
ylabel('Output Noise PSD (dB/Hz)');
title('Output Noise PSD after RC High-Pass and LPF (W = 10 Hz)');
xline(-W, '--'); xline(W, '--');
legend('S_{out}(f)', 'LPF Cutoff', 'Location', 'best');

%% 9
Ac  = 2;          
fc  = 1;          
N0  = 0.1;       
R   = 1; 
C   = 1;

Fs  = 200;       
Ts  = 1/Fs;
T   = 50;         
t   = 0:Ts:T-Ts;

Wvec = linspace(fc*1.05, 20, 30);  
Nmc  = 100;       % Monte Carlo runs
 
s = Ac*cos(2*pi*fc*t);

SNR_MC = zeros(size(Wvec));

for k = 1:length(Wvec)
    W = Wvec(k);
    SNR_runs = zeros(Nmc,1);

    for m = 1:Nmc

        % Generate AWGN 
        noise_var = N0*Fs/2;
        n = sqrt(noise_var)*randn(size(t));

        x = s + n;

        Nfft = length(t);
        f = (-Nfft/2:Nfft/2-1)*(Fs/Nfft);

        X = fftshift(fft(x));
        S = fftshift(fft(s));
        Nn = fftshift(fft(n));

        jwRC = 1j*2*pi*f*R*C;
        H_RC = jwRC ./ (1 + jwRC);

        H_LPF = double(abs(f) <= W);

        % Total transfer function
        V = H_RC .* H_LPF;

        % Filter signal and noise 
        Ys = ifft(ifftshift(S .* V));
        Yn = ifft(ifftshift(Nn .* V));

        Ps = mean(real(Ys).^2);
        Pn = mean(real(Yn).^2);

        SNR_runs(m) = Ps / Pn;
    end

    % Average over realizations
    SNR_MC(k) = mean(SNR_runs);
end

%  THEORETICAL SNR
a = 2*pi*R*C;
num = (pi*R*C*fc^2) / (1 + (a*fc)^2);
den = (a.*Wvec) - atan(a.*Wvec);
SNR_th = (Ac^2/N0) * (num ./ den);

figure;
plot(Wvec, 10*log10(SNR_MC), 'o', 'LineWidth', 1.2); hold on;
plot(Wvec, 10*log10(SNR_th), '-', 'LineWidth', 1.5);
grid on;
xlabel('LPF Bandwidth W (Hz)');
ylabel('SNR (dB)');
title('Monte Carlo vs Theoretical SNR');
legend('Monte Carlo SNR','Theoretical SNR','Location','best');

%% 10

R = 1; C = 1;          % RC parameters
W_LPF = 10;            % LPF bandwidth 
Fs = 200;              
f = linspace(0, 5*W_LPF, 1e5);  

H_RC_sq = (2*pi*f*R*C).^2 ./ (1 + (2*pi*f*R*C).^2);

H_LPF_sq = double(f <= W_LPF);

H_tot_sq = H_RC_sq .* H_LPF_sq;

Beq = trapz(f, H_tot_sq);  % Integrate |H(f)|^2 over frequency

fprintf('Equivalent Noise Bandwidth Beq = %.4f Hz\n', Beq);
