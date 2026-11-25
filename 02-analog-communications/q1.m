clc;
clear;
close all;
%%
fs = 1000;                 
Ts = 1/fs;
t = -1:Ts:10;             
fc = 100;                  
ka = 0.5;                
u4 = heaviside(t - 4);
u8 = heaviside(t - 8);
m = exp(-6*t).*(u4 - u8) + exp(6*t).*(u4 - u8);

%  TIME DOMAIN PLOT
figure;
plot(t, m, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Message Signal m(t)');
grid on;

%  FREQUENCY DOMAIN OF MESSAGE
N = length(m);
M_f = fftshift(fft(m));
f = (-N/2:N/2-1)*(fs/N);
M_mag = abs(M_f);
M_mag = M_mag / max(M_mag);
figure;
plot(f, M_mag, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
title('Spectrum of Message Signal');
grid on;
xlim([-20 20]);

%  BANDWIDTH ESTIMATION (−20 dB)
threshold = 0.1;                   % −20 dB threshold
idx = find(M_mag >= threshold);
B_20dB = max(abs(f(idx)));
fprintf('Estimated Bandwidth (−20 dB method): %.2f Hz\n', B_20dB);

%  ENERGY-BASED BANDWIDTH (95%)
E_total = sum(M_mag.^2);
E_cum = cumsum(M_mag.^2);
E_cum = E_cum / E_total;
idx95 = find(E_cum >= 0.95, 1);
B_95 = abs(f(idx95));
fprintf('Estimated Bandwidth (95%% energy): %.2f Hz\n', B_95);

%  AM MODULATION
carrier = cos(2*pi*fc*t);
s_am = (1 + ka*m).*carrier;
%  TIME DOMAIN AM SIGNAL
figure;
plot(t, s_am, 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('AM Modulated Signal');
grid on;
xlim([3.5 8.5]);
%  FREQUENCY DOMAIN OF AM SIGNAL
S_f = fftshift(fft(s_am));
S_mag = abs(S_f);
S_mag = S_mag / max(S_mag);
figure;
plot(f, S_mag, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
title('Spectrum of AM Signal');
grid on;
xlim([-200 200]);
%  DISPLAY AM BANDWIDTH
AM_Bandwidth = 2 * B_20dB;
fprintf('AM Signal Bandwidth ≈ %.2f Hz\n', AM_Bandwidth);

% 3
% ENVELOPE DETECTOR (RC)
RC = 0.05;                 % Time constant (seconds)
alpha = exp(-Ts/RC);       
env = zeros(size(s_am));
env(1) = 0;                
for n = 2:length(s_am)
    rectified = abs(s_am(n));  
    if rectified > env(n-1)
        env(n) = rectified;    
    else
        env(n) = alpha*env(n-1);
    end
end

% Time-domain envelope output
figure;
plot(t, env, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Envelope Detector Output');
grid on;
xlim([3.5 8.5]);
% Frequency-domain envelope output
ENV_f = fftshift(fft(env));
ENV_mag = abs(ENV_f)/max(abs(ENV_f));
figure;
plot(f, ENV_mag, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
title('Spectrum of Envelope Detector Output');
grid on;
xlim([-20 20]);

% 4
% ABS(.) + LPF DEMODULATION
rect_sig = abs(s_am);
figure;
plot(t, rect_sig, 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Rectified AM Signal');
grid on;
xlim([3.5 8.5]);

% Frequency domain after rectification
RECT_f = fftshift(fft(rect_sig));
RECT_mag = abs(RECT_f)/max(abs(RECT_f));
figure;
plot(f, RECT_mag, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
title('Spectrum After Rectification');
grid on;
xlim([-200 200]);

% Low-pass filter 
lpf_cutoff = 20;  
[b,a] = butter(4, lpf_cutoff/(fs/2));
m_rec_abs = filtfilt(b,a,rect_sig);

figure;
plot(t, m_rec_abs, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Recovered Signal After LPF');
grid on;
xlim([3.5 8.5]);

% Frequency domain after LPF
MABS_f = fftshift(fft(m_rec_abs));
MABS_mag = abs(MABS_f)/max(abs(MABS_f));
figure;
plot(f, MABS_mag, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
title('Spectrum After LPF');
grid on;
xlim([-20 20]);

%5
coherent_mix = s_am .* cos(2*pi*fc*t);

% Time domain: Mixer output
figure;
plot(t, coherent_mix, 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Coherent Detection: Mixer Output');
grid on;
xlim([3.5 8.5]);

% Frequency domain: Mixer output
COH_MIX_f = fftshift(fft(coherent_mix));
COH_MIX_mag = abs(COH_MIX_f) / max(abs(COH_MIX_f));
figure;
plot(f, COH_MIX_mag, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
title('Spectrum of Mixer Output');
grid on;
xlim([-200 200]);

% Low-pass filtering
m_coherent = filtfilt(b, a, coherent_mix);
figure;
plot(t, m_coherent, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Recovered Signal (Coherent Detection)');
grid on;
xlim([3.5 8.5]);
MCOH_f = fftshift(fft(m_coherent));
MCOH_mag = abs(MCOH_f) / max(abs(MCOH_f));
figure;
plot(f, MCOH_mag, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
title('Spectrum of Coherently Demodulated Signal');
grid on;
xlim([-20 20]);

% 6
% amdemod
ini_phase = 0;  
m_rec_amdemod = amdemod(s_am, fc, fs, ini_phase, ka);
% TIME DOMAIN PLOT
figure;
plot(t, m_rec_amdemod, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Recovered Signal using amdemod');
grid on;
xlim([3.5 8.5]);
% FREQUENCY DOMAIN PLOT
M_AMDEM_f = fftshift(fft(m_rec_amdemod));
M_AMDEM_mag = abs(M_AMDEM_f)/max(abs(M_AMDEM_f));
figure;
plot(f, M_AMDEM_mag, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
title('Spectrum of amdemod Recovered Signal');
grid on;
xlim([-20 20]);
% COMPARISON 
figure;
plot(t, env, 'r', 'LineWidth', 1); hold on;
plot(t, m_rec_abs, 'g', 'LineWidth', 1);
plot(t, m_coherent, 'b', 'LineWidth', 1);
plot(t, m_rec_amdemod, 'k--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Comparison of Demodulation Methods');
legend('Envelope Detector', 'Abs + LPF', 'Coherent Detection', 'amdemod');
grid on;
xlim([3.5 8.5]);
env_f = abs(fftshift(fft(env)))/max(abs(fftshift(fft(env))));
m_abs_f = abs(fftshift(fft(m_rec_abs)))/max(abs(fftshift(fft(m_rec_abs))));
m_coh_f = abs(fftshift(fft(m_coherent)))/max(abs(fftshift(fft(m_coherent))));
m_amdemod_f = abs(fftshift(fft(m_rec_amdemod)))/max(abs(fftshift(fft(m_rec_amdemod))));
figure;
plot(f, env_f, 'r', 'LineWidth', 1); hold on;
plot(f, m_abs_f, 'g-o', 'LineWidth', 1, 'MarkerIndices', 1:50:length(f)); % markers
plot(f, m_coh_f, 'b', 'LineWidth', 1);
plot(f, m_amdemod_f, 'k--', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
title('Frequency Domain Comparison of Demodulation Methods');
legend('Envelope Detector', 'Abs + LPF', 'Coherent Detection', 'amdemod');
grid on;
xlim([-20 20]);

% 7 
% DSB-SC MODULATION & COHERENT DEMODULATION
s_dsb = m .* cos(2*pi*fc*t);
% Time-domain
figure;
plot(t, s_dsb, 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Amplitude');
title('DSB-SC Modulated Signal'); grid on; xlim([3.5 8.5]);
% Frequency-domain
S_DSB_f = fftshift(fft(s_dsb));
S_DSB_mag = abs(S_DSB_f)/max(abs(S_DSB_f));
figure;
plot(f, S_DSB_mag, 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Normalized Magnitude');
title('Spectrum of DSB-SC Signal'); grid on; xlim([-200 200]);

m_dsb_coh = filtfilt(b,a, s_dsb .* 2 .* cos(2*pi*fc*t));
% Time-domain recovered
figure;
plot(t, m_dsb_coh, 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Amplitude');
title('Recovered Signal after DSB-SC Coherent Demodulation'); grid on; xlim([3.5 8.5]);
% Frequency-domain recovered
M_DSBCOH_f = fftshift(fft(m_dsb_coh));
M_DSBCOH_mag = abs(M_DSBCOH_f)/max(abs(M_DSBCOH_f));
figure;
plot(f, M_DSBCOH_mag, 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Normalized Magnitude');
title('Spectrum after DSB-SC Coherent Demodulation'); grid on; xlim([-20 20]);

% 8
% USSB MODULATION & COHERENT DEMODULATION 
m_analytic = hilbert(m);  % m(t) + j*H{m(t)}
s_ussb = real(m_analytic .* exp(1j*2*pi*fc*t));
% Time-domain modulated signal
figure;
plot(t, s_ussb, 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Amplitude');
title('USSB Modulated Signal'); grid on;
xlim([3.5 8.5]);
% Frequency-domain modulated signal
S_USSB_f = fftshift(fft(s_ussb));
S_USSB_mag = abs(S_USSB_f)/max(abs(S_USSB_f));
figure;
plot(f, S_USSB_mag, 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Normalized Magnitude');
title('Spectrum of USSB Signal'); grid on;
xlim([-200 200]);

m_ussb_rec = filtfilt(b,a, 2*s_ussb.*cos(2*pi*fc*t));
% Time-domain recovered signal
figure;
plot(t, m_ussb_rec, 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Amplitude');
title('Recovered Signal after USSB Demodulation'); grid on;
xlim([3.5 8.5]);
% Frequency-domain recovered signal
M_USSB_f = fftshift(fft(m_ussb_rec));
M_USSB_mag = abs(M_USSB_f)/max(abs(M_USSB_f));
figure;
plot(f, M_USSB_mag, 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Normalized Magnitude');
title('Spectrum after USSB Demodulation'); grid on;
xlim([-20 20]);

% 9 LSSB MODULATION & COHERENT DEMODULATION 
m_analytic = hilbert(m);  % m(t) + j*hilbert(m(t))
% LSSB modulation (lower sideband)
s_lssb = real(m_analytic .* exp(1j*2*pi*fc*t));
% Time-domain
figure;
plot(t, s_lssb, 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('LSSB Modulated Signal');
grid on;
xlim([3.5 8.5]);
% Frequency-domain
S_LSSB_f = fftshift(fft(s_lssb));
S_LSSB_mag = abs(S_LSSB_f)/max(abs(S_LSSB_f));
figure;
plot(f, S_LSSB_mag, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
title('Spectrum of LSSB Signal');
grid on;
xlim([-200 200]);
% Coherent demodulation
m_lssb_rec = filtfilt(b,a, 2*s_lssb.*cos(2*pi*fc*t));
% Time-domain recovered
figure;
plot(t, m_lssb_rec, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Amplitude');
title('Recovered Signal after LSSB Demodulation ');
grid on;
xlim([3.5 8.5]);
% Frequency-domain recovered
M_LSSB_f = fftshift(fft(m_lssb_rec));
M_LSSB_mag = abs(M_LSSB_f)/max(abs(M_LSSB_f));
figure;
plot(f, M_LSSB_mag, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
title('Spectrum after LSSB Demodulation');
grid on;
xlim([-20 20]);
%% 10 
%COMPARISON OF MODULATION METHODS
fprintf('\nComparison of Modulation Techniques:\n');
fprintf('1. AM (Envelope/Abs+LPF): Simple, requires no carrier synchronization, slightly less accurate.\n');
fprintf('2. AM (Coherent/AMDEM0D): Most accurate recovery, requires precise carrier phase knowledge.\n');
fprintf('3. DSB-SC: Efficient bandwidth usage (double sideband only), but needs coherent demodulation.\n');
fprintf('4. USSB: Only upper sideband transmitted → halves bandwidth compared to DSB-SC, requires coherent detection.\n');
fprintf('5. LSSB: Only lower sideband transmitted → also halves bandwidth, requires coherent detection.\n\n');

fprintf('Summary:\n');
fprintf('- AM is simple but less bandwidth efficient.\n');
fprintf('- DSB-SC is bandwidth efficient but needs synchronization.\n');
fprintf('- SSB (USSB/LSSB) is most bandwidth efficient, reduces power usage, but more complex to implement.\n');
fprintf('- Coherent methods produce most accurate recovery but require carrier phase knowledge.\n');

% Recovered signals 
rec_am_env     = env;              % Envelope detector
rec_am_abs     = m_rec_abs;        % Abs + LPF
rec_am_coh     = m_coherent;       % Coherent AM
rec_am_builtin = m_rec_amdemod;    % amdemod
rec_dsb        = m_dsb_coh;         % DSB-SC
rec_ussb       = m_ussb_rec;        % USSB
rec_lssb       = m_lssb_rec;        % LSSB

figure; hold on;
plot(t, rec_am_env,     'r-',  'LineWidth', 1.2);
plot(t, rec_am_abs,     'g--', 'LineWidth', 1.2);
plot(t, rec_am_coh,     'b-.', 'LineWidth', 1.2);
plot(t, rec_am_builtin, 'k:',  'LineWidth', 1.8);
plot(t, rec_dsb,        'm-',  'LineWidth', 1.2);
plot(t, rec_ussb,       'c--', 'LineWidth', 1.2);
plot(t, rec_lssb,       'y-.', 'LineWidth', 1.2);

xlabel('Time (s)');
ylabel('Amplitude');
title('Time-Domain Comparison of Recovered Signals');
legend('AM Envelope', 'AM Abs + LPF', 'AM Coherent', 'AM amdemod', ...
       'DSB-SC', 'USSB', 'LSSB', 'Location', 'best');
grid on;
xlim([3.5 8.5]);

figure; hold on;

signals = {rec_am_env, rec_am_abs, rec_am_coh, rec_am_builtin, ...
           rec_dsb, rec_ussb, rec_lssb};
styles = {'r-', 'g--', 'b-.', 'k:', 'm-', 'c--', 'y-.'};

for k = 1:length(signals)
    Sig_f = abs(fftshift(fft(signals{k})));
    Sig_f = Sig_f / max(Sig_f);
    plot(f, Sig_f, styles{k}, 'LineWidth', 1.2);
end

xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
title('Frequency-Domain Comparison of Recovered Signals');
legend('AM Envelope', 'AM Abs + LPF', 'AM Coherent', 'AM amdemod', ...
       'DSB-SC', 'USSB', 'LSSB', 'Location', 'best');
grid on;
xlim([-20 20]);

%% 11

clc;

% --- PARAMETERS ---
fs_audio = 16000;        % Sampling frequency (Hz)
record_duration = 5;     % Recording time (seconds)
fc_audio = 3000;         % Carrier frequency (Hz) < fs/2
ka_audio = 0.8;          % Modulation index

fprintf('Recording will start now...\n');

% --- RECORD VOICE ---
recObj = audiorecorder(fs_audio, 16, 1);
recordblocking(recObj, record_duration);
fprintf('Recording complete.\n');

voice = getaudiodata(recObj);
voice = voice / max(abs(voice));   % Normalize

t_audio = (0:length(voice)-1)/fs_audio;

% --- PLAY ORIGINAL VOICE ---
sound(voice, fs_audio);
pause(record_duration + 1);

% --- SAVE AUDIO (OPTIONAL) ---
audiowrite('recorded_voice.wav', voice, fs_audio);

% --- TIME DOMAIN: ORIGINAL VOICE ---
figure;
plot(t_audio, voice, 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Original Voice Signal (Time Domain)');
grid on;

% --- FREQUENCY DOMAIN: ORIGINAL VOICE ---
N = length(voice);
f_audio = (-N/2:N/2-1)*(fs_audio/N);
VOICE_f = fftshift(fft(voice));
VOICE_mag = abs(VOICE_f)/max(abs(VOICE_f));

figure;
plot(f_audio, VOICE_mag, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
title('Original Voice Signal (Frequency Domain)');
grid on;
xlim([-4000 4000]);

% --- AM MODULATION ---
carrier_audio = cos(2*pi*fc_audio*t_audio');
s_am_audio = (1 + ka_audio*voice) .* carrier_audio;

% --- TIME DOMAIN: AM SIGNAL ---
figure;
plot(t_audio, s_am_audio, 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('AM Modulated Voice Signal (Time Domain)');
grid on;

% --- FREQUENCY DOMAIN: AM SIGNAL ---
S_AM_f = fftshift(fft(s_am_audio));
S_AM_mag = abs(S_AM_f)/max(abs(S_AM_f));

figure;
plot(f_audio, S_AM_mag, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
title('AM Modulated Voice Signal (Frequency Domain)');
grid on;
xlim([-6000 6000]);

% --- AM DEMODULATION (COHERENT) ---
demod_mix = s_am_audio .* cos(2*pi*fc_audio*t_audio');

% Low-pass filter
lpf_cutoff = 4000;  % Hz (speech bandwidth)
[b_audio, a_audio] = butter(4, lpf_cutoff/(fs_audio/2));

voice_rec = filtfilt(b_audio, a_audio, demod_mix);
voice_rec = voice_rec / max(abs(voice_rec));  % Normalize

% --- PLAY RECOVERED VOICE ---
sound(voice_rec, fs_audio);

% --- TIME DOMAIN: DEMODULATED VOICE ---
figure;
plot(t_audio, voice_rec, 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Recovered Voice Signal (Time Domain)');
grid on;

% --- FREQUENCY DOMAIN: DEMODULATED VOICE ---
VOICE_REC_f = fftshift(fft(voice_rec));
VOICE_REC_mag = abs(VOICE_REC_f)/max(abs(VOICE_REC_f));

figure;
plot(f_audio, VOICE_REC_mag, 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
title('Recovered Voice Signal (Frequency Domain)');
grid on;
xlim([-4000 4000]);

% --- QUALITY EVALUATION ---
fprintf('\nPart 11 – Audio Quality Evaluation:\n');
fprintf('The recovered voice signal closely resembles the original recording.\n');
fprintf('Some high-frequency attenuation and slight distortion may be observed.\n');
fprintf('Overall speech intelligibility is preserved using AM coherent demodulation.\n');
