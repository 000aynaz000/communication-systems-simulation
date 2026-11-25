% q3
clear; close all; clc;

%% (a) Load Original Signal
[x, fs] = audioread('my_voice.wav');   
if size(x,2) > 1
    x = mean(x,2);                     
end
x = x(:);                              
t = (0:length(x)-1)/fs;

% normalize 
x = x / max(abs(x));


%% (b) The Muffling Channel
h = [0.6, 0.3, 0.1];   % channel impulse response (FIR)
nfft = 4096;
[H, f_H] = freqz(h,1,nfft,fs);        % frequency response sampled over [0,fs)
Hmag = abs(H);

figure('Name','Channel Frequency Response');
plot(f_H, Hmag, 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('|H(f)|');
title('Magnitude Response of Channel h[n] = 0.6δ[n]+0.3δ[n-1]+0.1δ[n-2]');
grid on;
xlim([0 fs/2]);

figure('Name','Channel Frequency Response (dB)');
plot(f_H, 20*log10(Hmag + eps), 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Channel Magnitude Response (dB)');
grid on;
xlim([0 fs/2]);

fprintf('Observation: this channel attenuates high frequencies — it behaves like a low-pass filter.\n');

%% (c) Create Distorted (Muffled) Signal
y = conv(x, h, 'same');
y = y / max(abs(y));   

%% (d) Compare Spectrograms of x[n] and y[n]
win = 1024;                 
overlap = round(0.5*win);  
nfft_spec = 2048;

figure('Name','Spectrogram - Original vs Muffled','Position',[100 100 1200 400]);
subplot(1,2,1);
spectrogram(x, win, overlap, nfft_spec, fs, 'yaxis');
title('Spectrogram of Original x[n]');
ylim([0 8]); 
subplot(1,2,2);
spectrogram(y, win, overlap, nfft_spec, fs, 'yaxis');
title('Spectrogram of Muffled y[n]');

figs = findall(gcf,'Type','Axes');

%% (e,f) Design Equalizer 

Nf = 512;
fvec = linspace(0, fs/2, Nf)';   

% Compute H(f) at these frequencies
Hf = freqz(h,1,fvec,fs);         
Hf = Hf(:);

Gdesired = 1.0 ./ (Hf + 0.1);   
F = linspace(0,1,length(Gdesired));
A = abs(Gdesired)';   
order = 50;
g = firls(order, F, A);   % g will be of length order+1 = 51
[Gw, fGw] = freqz(g,1,nfft,fs);

figure('Name','Channel and Equalizer Frequency Responses','Position',[100 100 700 450]);
plot(f_H, 20*log10(abs(H)+eps), 'LineWidth', 1.2); hold on;
plot(fGw, 20*log10(abs(Gw)+eps), 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('|H(f)| (channel)','|G(f)| (designed equalizer)');
title('Channel and Equalizer Frequency Responses (dB)');
xlim([0 fs/2]);
grid on;

%% (g) Apply Equalizer to muffled signal
xrec = conv(y, g, 'same');
xrec = xrec / max(abs(xrec) + eps);

%% Diagnostic & Improved Final Comparison (replace previous (h) section)

N = min([length(x), length(y), length(xrec)]);
x1 = x(1:N);
y1 = y(1:N);
xr1 = xrec(1:N);

% Spectrogram parameters
win = 512; overlap = round(0.5*win); nfft_spec = 1024;

% Compute STFTs 
[Sx, F, T]   = spectrogram(x1, win, overlap, nfft_spec, fs);
[Sy, ~, ~]   = spectrogram(y1, win, overlap, nfft_spec, fs);
[Sr, ~, ~]   = spectrogram(xr1, win, overlap, nfft_spec, fs);

% Convert to dB (
Sx_dB = 20*log10(abs(Sx)+eps);
Sy_dB = 20*log10(abs(Sy)+eps);
Sr_dB = 20*log10(abs(Sr)+eps);

% Compute common dynamic range but keep sensitivity 
maxAll = max([Sx_dB(:); Sy_dB(:); Sr_dB(:)]);
cmin = maxAll - 60;
cmax = maxAll;

%Plot 
figure('Name','Spectrograms (diagnostic)','Position',[100 100 1400 420]);
subplot(1,3,1);
imagesc(T, F, Sx_dB);
axis xy; colormap hot; caxis([cmin cmax]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('Original x[n]');
subplot(1,3,2);
imagesc(T, F, Sy_dB);
axis xy; colormap hot; caxis([cmin cmax]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('Muffled y[n]');
subplot(1,3,3);
imagesc(T, F, Sr_dB);
axis xy; colormap hot; caxis([cmin cmax]);
xlabel('Time (s)'); ylabel('Frequency (Hz)'); title('Recovered x_{rec}[n]');
colorbar('Position',[0.93 0.11 0.02 0.8]);
drawnow;

%Difference spectrogram 
Diff_dB = Sr_dB - Sy_dB;
figure('Name','Difference Spectrogram (Recovered - Muffled)','Position',[200 200 800 400]);
imagesc(T, F, Diff_dB); axis xy; colormap jet;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title('Difference: 20log10|S_{rec}| - 20log10|S_{y}|');
caxis([-20 20]); colorbar;
drawnow;

%  Frequency-domain average magnitude 
nfft = 2^nextpow2(N);
X = fft(x1, nfft);
Y = fft(y1, nfft);
XR = fft(xr1, nfft);
faxis = (0:nfft-1)/nfft*fs;

% Compute single-sided magnitude 
Xmag = 20*log10(abs(X)+eps);
Ymag = 20*log10(abs(Y)+eps);
Rmag = 20*log10(abs(XR)+eps);

% Plot only 0..fs/2
half = 1:floor(nfft/2);
figure('Name','Average Magnitude Spectra','Position',[200 200 700 450]);
plot(faxis(half), Xmag(half),'k','LineWidth',1.2); hold on;
plot(faxis(half), Ymag(half),'r','LineWidth',1.0);
plot(faxis(half), Rmag(half),'g','LineWidth',1.0);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
legend('Original','Muffled','Recovered'); grid on;
xlim([0 fs/2]);
drawnow;

% Compute numeric metric: energy above cutoff frequency
fc = 3000; 
idx_hf = find(faxis(half) >= fc);

% Compute band energy 
E_orig_hf = sum(abs(X(half(idx_hf))).^2);
E_muff_hf = sum(abs(Y(half(idx_hf))).^2);
E_rec_hf  = sum(abs(XR(half(idx_hf))).^2);

fprintf('\nHigh-frequency energy (f >= %d Hz):\n', fc);
fprintf(' Original: %.4e\n', E_orig_hf);
fprintf(' Muffled : %.4e (ratio to original: %.3f)\n', E_muff_hf, E_muff_hf/E_orig_hf);
fprintf(' Recovered: %.4e (ratio to original: %.3f, relative to muffled: %.3f)\n', E_rec_hf, E_rec_hf/E_orig_hf, E_rec_hf/E_muff_hf);

%Compute MSE and correlation between original and recovered 
mse_xr = mean((x1 - xr1).^2);
corr_xr = corrcoef(x1, xr1);
if numel(corr_xr) >= 2, corr_xr = corr_xr(1,2); else corr_xr = NaN; end
fprintf('\nTime-domain metrics:\n MSE(original, recovered) = %.4e\n Corr(original, recovered) = %.4f\n', mse_xr, corr_xr);

%% Save figures 
exportgraphics(findall(0,'Name','Channel Frequency Response'),'fig_channel_mag.png','Resolution',300);
exportgraphics(findall(0,'Name','Channel Frequency Response (dB)'),'fig_channel_db.png','Resolution',300);
exportgraphics(findall(0,'Name','Channel and Equalizer Frequency Responses'),'fig_resp.png','Resolution',300);
exportgraphics(findall(0,'Name','Spectrogram - Original vs Muffled'),'fig_spec_xy.png','Resolution',300);
exportgraphics(findall(0,'Name','Spectrograms: Original - Muffled - Recovered'),'fig_spec_all.png','Resolution',300);

fprintf('Design complete. Figures exported (PNG). Use Publish or Live Script to create PDF report.\n');
fprintf('\nShort results summary:\n');
fprintf('- Channel h[n] is low-pass (attenuates high-frequency components).\n');
fprintf('- The muffled signal y[n] shows reduced energy at high frequencies in the spectrogram.\n');
fprintf('- The FIR equalizer (51 taps) approximates 1/H(f) (regularized by +0.1).\n');
fprintf('- Recovered xrec[n] typically restores some high-frequency content, but may amplify noise and create ringing/artifacts.\n');
fprintf('  Listen and inspect spectrograms to evaluate quality (possible trade-offs: noise amplification, pre/post ringing).\n');

