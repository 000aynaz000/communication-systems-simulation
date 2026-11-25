%% Q4: Sampling Theorem, Aliasing, and DTMF Decoding
clear; close all; clc;
%% (a) Create continuous-time like signal x(t)
frame_dur = 0.08;           
T = 1;                 
fs_high = 50000;       
t = 0:1/fs_high:(T-1/fs_high);
f0 = 300;               
x = cos(2*pi*f0*t);
N = length(x);
Nfft = 2^nextpow2(N);
X = fftshift(fft(x, Nfft));
f_axis = linspace(-fs_high/2, fs_high/2, Nfft);
figure;
plot(f_axis, abs(X));
xlim([-5000 5000]);
xlabel('Frequency (Hz)');
ylabel('|X(f)|');
title('Magnitude Spectrum of x(t) = cos(2\pi 300 t) (continuous representation)');
grid on;
%% (b) Case 1: Good sampling fs = 1000 Hz (fs > 2*f0)
fs1 = 1000;
sample_factor1 = round(fs_high / fs1); 
x1 = x(1:sample_factor1:end);           
n1 = 0:length(x1)-1;
t1 = n1 / fs1;
N1 = length(x1);
Nfft1 = 2^nextpow2(4*N1);               
X1 = fftshift(fft(x1, Nfft1));
f1_axis = linspace(-fs1/2, fs1/2, Nfft1);
figure;
plot(f1_axis, abs(X1));
xlim([-fs1/2 fs1/2]);
xlabel('Frequency (Hz)');
ylabel('|X_1(f)|');
title(sprintf('Magnitude Spectrum of sampled x1[n] (fs = %d Hz) from %d to %d Hz', fs1, -fs1/2, fs1/2));
grid on;
[~, idx_max1] = max(abs(X1));
peak_freq1 = f1_axis(idx_max1);
fprintf('Good sampling: strongest peak at %0.2f Hz (expected near +/-300 Hz if in band).\n', peak_freq1);
%% (d) Case 2: Bad sampling fs = 400 Hz (aliasing expected)
fs2 = 400;
sample_factor2 = round(fs_high / fs2);  
x2 = x(1:sample_factor2:end);
n2 = 0:length(x2)-1;
t2 = n2 / fs2;
% (e) Spectrum of bad sampling
N2 = length(x2);
Nfft2 = 2^nextpow2(4*N2);
X2 = fftshift(fft(x2, Nfft2));
f2_axis = linspace(-fs2/2, fs2/2, Nfft2);
figure;
plot(f2_axis, abs(X2));
xlim([-fs2/2 fs2/2]);
xlabel('Frequency (Hz)');
ylabel('|X_2(f)|');
title(sprintf('Magnitude Spectrum of sampled x2[n] (fs = %d Hz) from %d to %d Hz', fs2, -fs2/2, fs2/2));
grid on;
% Find strongest peaks
[~, idx_max2] = max(abs(X2));
peak_freq2 = f2_axis(idx_max2);
fprintf('Bad sampling: strongest peak at %0.2f Hz (expected alias at +/-100 Hz).\n', peak_freq2);
%% (f) Analysis: compute alias analytically
f_alias = abs(f0 - fs2); 
fprintf('Analytical aliased frequency: falias = |%d - %d| = %d Hz\n', f0, fs2, f_alias);
fprintf('Does the plotted spectrum show peaks at +/-%d Hz? Inspect the plot and the printed peak above.\n', f_alias);
%% (g)-(h) DTMF Decoding
[dtmf_sig, fs_dtmf] = audioread('dtmf_sequence.wav'); % ensure file present
% If stereo, convert to mono
if size(dtmf_sig,2) > 1
    dtmf_sig = mean(dtmf_sig, 2);
end
figure;
plot((0:length(dtmf_sig)-1)/fs_dtmf, dtmf_sig);
xlabel('Time (s)');
ylabel('Amplitude');
title('DTMF signal waveform');
xlim([0 (length(dtmf_sig)-1)/fs_dtmf]);
grid on;
% (h) Spectrogram
win_len = 0.03 * fs_dtmf;    
win = hamming(round(win_len));
noverlap = round(0.02 * fs_dtmf); 
nfft = 4096;
figure;
[S,F,T,P] = spectrogram(dtmf_sig, win, noverlap, nfft, fs_dtmf, 'yaxis'); 
imagesc(T, F, 20*log10(abs(S))); axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram of DTMF signal (dB)');
ylim([500 1600]);             
colorbar;
% The spectrogram shows bright horizontal lines near the DTMF frequencies for each tone.
%% Automated decoding procedure:
row_freqs = [697, 770, 852, 941];         
col_freqs = [1209, 1336, 1477];           
keys = {'1','2','3'; '4','5','6'; '7','8','9'; '*','0','#'};
frame_len = round(frame_dur * fs_dtmf);
hop = round(frame_len/2);
siglen = length(dtmf_sig);
decoded_seq = '';
i = 1;
min_energy_threshold = 0.001; 
while (i + frame_len - 1) <= siglen
    frame = dtmf_sig(i:(i+frame_len-1)) .* hamming(frame_len);
    if sum(frame.^2) < min_energy_threshold
        i = i + hop;
        continue;
    end
    
    Nf = 2^nextpow2(4*frame_len);
    Fframe = fft(frame, Nf);
    mag = abs(Fframe(1:Nf/2+1));
    f_ax = linspace(0, fs_dtmf/2, length(mag));
    
    band_idx = find(f_ax >= 650 & f_ax <= 1550);
    [peaks_sorted, idx_sorted] = sort(mag(band_idx), 'descend');
    
    if isempty(peaks_sorted)
        i = i + hop;
        continue;
    end
    
    peak1_idx = band_idx(idx_sorted(1));
    found = false;
    for k = 2:min(8,length(idx_sorted))
        cand_idx = band_idx(idx_sorted(k));
        if abs(f_ax(cand_idx) - f_ax(peak1_idx)) > 70
            peak2_idx = cand_idx;
            found = true;
            break;
        end
    end
    if ~found
        i = i + hop;
        continue;
    end
    
    f_peak1 = f_ax(peak1_idx);
    f_peak2 = f_ax(peak2_idx);
    % sort so low freq first
    fpair = sort([f_peak1, f_peak2]);
    
    if abs(fpair(1) - row_freqs(1)) < 200 
        lowf = fpair(1);
        highf = fpair(2);
    else
        d1_row = min(abs(fpair(1)-row_freqs));
        d1_col = min(abs(fpair(1)-col_freqs));
        d2_row = min(abs(fpair(2)-row_freqs));
        d2_col = min(abs(fpair(2)-col_freqs));
        if d1_row < d1_col
            lowf = fpair(1); highf = fpair(2);
        else
            lowf = fpair(2); highf = fpair(1);
        end
    end
    
    [~, ridx] = min(abs(row_freqs - lowf));
    [~, cidx] = min(abs(col_freqs - highf));
    
    detected_key = keys{ridx, cidx};
    
    if isempty(decoded_seq) || decoded_seq(end) ~= detected_key
        decoded_seq(end+1) = detected_key; 
    end
    
    i = i + round(0.08 * fs_dtmf); 
end
fprintf('Decoded DTMF sequence (automated): %s\n', decoded_seq);
%% Display instruction for manual verification:
fprintf('\nManual verification:\n');
fprintf(' - Inspect the spectrogram (y-axis frequencies) and for each vertical block read the two bright lines (one ~700-950 Hz, one ~1200-1500 Hz)\n');
fprintf(' - Use the keypad table to map the pair to a digit (table used in code).\n');