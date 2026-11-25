%% question 2
clear; close all; clc;
%% (a) Load and plot
[ y, fs ] = audioread('corrupted_speech.wav');  
if size(y,2) > 1
    y = mean(y,2);   
end
N = length(y);
t = (0:N-1)/fs;
figure('Name','Time-domain signal');
plot(t, y);
xlabel('Time (s)');
ylabel('Amplitude');
title('(a) Noisy signal y[n] (time domain)');
grid on;
%% (b) Spectral Analysis: PSD estimate using Welch
% Parameters for PSD
win = hann(4096);         
noverlap = round(length(win)*0.5);
nfft = 2^nextpow2(4*length(win));
[ Pxx, f ] = pwelch(y, win, noverlap, nfft, fs, 'power');  
Pxx_dB = 10*log10(Pxx);
figure('Name','PSD (Welch)');
plot(f, Pxx_dB);
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('(b) PSD of noisy signal y[n]');
xlim([0 fs/2]);
grid on;
%% (c) Identify interference & speech band
minProm = 10;  
[pks, locs] = findpeaks(Pxx_dB, f, 'MinPeakProminence', minProm, 'MinPeakDistance', 10);
if isempty(pks)
    [pks, locs] = findpeaks(Pxx_dB, f, 'MinPeakProminence', 6, 'MinPeakDistance', 5);
end
if ~isempty(pks)
    [~, idxmax] = max(pks);
    fwhine = locs(idxmax);
else
    fwhine = NaN;
end
cdf = cumsum(Pxx);
cdf = cdf / cdf(end);
f95idx = find(cdf >= 0.95, 1, 'first');
if isempty(f95idx)
    f95 = f(end);
else
    f95 = f(f95idx);
end
fprintf('Detected whine frequency f_whine = %.1f Hz\n', fwhine);
fprintf('Frequency containing 95%% of signal power: %.1f Hz\n', f95);
figure(gcf); hold on;
if ~isnan(fwhine)
    plot(fwhine, Pxx_dB(f==fwhine),'ro','MarkerFaceColor','r');
    text(fwhine, Pxx_dB(find(f==fwhine,1)), sprintf('  f_{whine}=%.1f Hz', fwhine),'Color','r');
end
yline(10*log10(mean(Pxx)), '--k', 'Mean PSD');
xline(f95, '--g', sprintf('f_{95}=%.0f Hz', f95));
hold off;
%% (d) Filter 1: Notch filter (IIR notch using iirnotch) - remove the whine
if ~isnan(fwhine)
    w0 = fwhine / (fs/2);   
    Q = 60;                
    [b_notch, a_notch] = iirnotch(w0, w0/Q);
    y_notched = filtfilt(b_notch, a_notch, y);
else
    warning('No clear whine detected automatically. Skipping notch step; set fwhine manually if needed.');
    y_notched = y;
end
% Plot short segment to visualize removal
figure('Name','After Notch (time domain)');
plot(t, y_notched);
xlabel('Time (s)');
ylabel('Amplitude');
title('(d) y_{notched}[n] (time domain) after notch filter');
grid on;
%% (e) Filter 2: Low-pass FIR filter to remove high-frequency hiss
fcut = min(f95 * 1.10, fs/2 * 0.95);  
fcut = max(fcut, 3000); 
fprintf('Chosen LPF cutoff fcut = %.1f Hz\n', fcut);
% Design a linear-phase FIR lowpass filter using a Hamming window
firOrder = 512;              
Wn = fcut / (fs/2);          
b_lp = fir1(firOrder, Wn, 'low', hamming(firOrder+1));
z = filtfilt(b_lp, 1, y_notched);
figure('Name','Cleaned (time domain)');
plot(t, z);
xlabel('Time (s)');
ylabel('Amplitude');
title('(e) Final cleaned signal z[n] (time domain) after LPF');
grid on;
%% (f) Verification (PSD): plot PSD of final signal and compare
[ Pxx_z, fz ] = pwelch(z, win, noverlap, nfft, fs, 'power');
Pxx_z_dB = 10*log10(Pxx_z);
figure('Name','PSD Comparison');
plot(f, Pxx_dB, 'Color', [0.7 0.7 0.7]); hold on;    
plot(fz, Pxx_z_dB, 'b', 'LineWidth', 1.2);           
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('(f) PSD: original (grey) vs cleaned z[n] (blue)');
legend('Original y[n]', 'Cleaned z[n]');
xlim([0 fs/2]);
grid on;
if ~isnan(fwhine)
    figure('Name','Zoom on whine frequency');
    idx_zoom = f > (fwhine-200) & f < (fwhine+200);
    plot(f(idx_zoom), Pxx_dB(idx_zoom), '--k'); hold on;
    plot(fz(idx_zoom), Pxx_z_dB(idx_zoom), 'b', 'LineWidth', 1.2);
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    title(sprintf('(f2) Zoom near f_{whine}=%.1f Hz (original vs cleaned)', fwhine));
    legend('Original', 'Cleaned');
    grid on;
end
%% (g) Verification (Spectrogram): original vs cleaned side-by-side
spec_win = 1024;
spec_nover = round(spec_win*0.75);
spec_nfft = 2048;
figure('Name','Spectrogram: Original and Cleaned','Units','normalized','Position',[0.1 0.1 0.8 0.7]);
subplot(2,1,1);
spectrogram(y, spec_win, spec_nover, spec_nfft, fs, 'yaxis');
title('(g) Spectrogram of original noisy signal y[n]');
ylim([0 8]);  
colorbar off;
subplot(2,1,2);
spectrogram(z, spec_win, spec_nover, spec_nfft, fs, 'yaxis');
title('(g) Spectrogram of cleaned signal z[n]');
ylim([0 8]);
colorbar off;
%% (h) Listen and save
audiowrite('cleaned_speech.wav', z, fs);
fprintf('Saved cleaned audio to cleaned_speech.wav\n');
%% End: Short printed summary for report
fprintf('\nSummary for report:\n');
fprintf('- Detected whine frequency (Hz): %.1f\n', fwhine);
fprintf('- LPF cutoff (Hz) used: %.1f\n', fcut);
fprintf('- FIR order used: %d; Notch Q: %.1f\n', firOrder, Q);
fprintf('- Cleaned file saved as cleaned_speech.wav\n');