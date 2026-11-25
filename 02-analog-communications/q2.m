%% 2
fs  = 10e3;         
Ts  = 1/fs;
tend = 10;           
t = 0:Ts:tend;
fc = 500;            
fDelta = 60;         

% x(t) and its integral
[x, x_int] = generate_signal(t);
figure;
subplot(2,1,1)
plot(t,x,'LineWidth',1.5);
grid on; xlabel('t (s)'); ylabel('x(t)');
title('Message Signal x(t)');
subplot(2,1,2)
plot(t,x_int,'LineWidth',1.5);
grid on; xlabel('t (s)'); ylabel('\int x(\tau)d\tau');
title('Integral of x(t)');

% FM Modulation (Manual Implementation)
xc_manual = cos(2*pi*fc*t + 2*pi*fDelta*x_int);
twin = 0.1;
idx = t <= twin;
figure;
plot(t(idx), xc_manual(idx), 'LineWidth', 1.2);
grid on;
xlabel('t (s)'); ylabel('x_c(t)');
title(['FM Modulated Signal (Manual) — First ' num2str(twin) ' s']);

% FM Modulation using fmmod 
xc_fmmod = fmmod(x, fc, fs, fDelta);
figure;
plot(t(idx), xc_fmmod(idx), 'LineWidth', 1.2);
grid on;
xlabel('t (s)'); ylabel('x_c(t)');
title(['FM Modulated Signal (fmmod) — First ' num2str(twin) ' s']);

%Fourier of x(t) and x_c(t)
[freq, Xmag]   = compute_fft(x, fs);
[~,   Xc_mag]  = compute_fft(xc_manual, fs);
figure;
subplot(2,1,1)
plot(freq, Xmag,'LineWidth',1.2);
grid on; xlim([0 2000]);
xlabel('Frequency (Hz)'); ylabel('|X(f)|');
title('Spectrum of Message Signal x(t)');
subplot(2,1,2)
plot(freq, Xc_mag,'LineWidth',1.2);
grid on; xlim([0 2000]);
xlabel('Frequency (Hz)'); ylabel('|X_c(f)|');
title('Spectrum of FM Signal x_c(t)');


% Carson Bandwidth, Modulation Index β 
[maxDev, beta, fm_eff, BW_carson] = carson_bandwidth(x, fDelta, fs);
fprintf('\n===== FM PARAMETER RESULTS =====\n');
fprintf('Max value of x(t):           %.2f\n', max(x));
fprintf('Max frequency deviation:     %.2f Hz\n', maxDev);
fprintf('Effective highest baseband f_m: %.2f Hz\n', fm_eff);
fprintf('Modulation index (beta):     %.2f\n', beta);
fprintf('Carson bandwidth:            %.2f Hz\n', BW_carson);
fprintf('=================================\n\n');

% 6
twin = 0.1;
idx = t <= twin;
dxcdt = zeros(size(xc_manual));
dxcdt(2:end-1) = (xc_manual(3:end) - xc_manual(1:end-2)) ./ (2*Ts);
abs_sig = abs(dxcdt);
fc_lpf = max(1, 2*fm_eff);
yd = lowpass(abs_sig, fc_lpf, fs);

figure;
subplot(3,1,1)
plot(t(idx), dxcdt(idx),'LineWidth',1.2);
grid on; title('Derivative of x_c(t)');
xlabel('t (s)'); ylabel('d/dt\{x_c(t)\}');
subplot(3,1,2)
plot(t(idx), abs_sig(idx),'LineWidth',1.2);
grid on; title('| d/dt(x_c(t)) |');
xlabel('t (s)');
subplot(3,1,3)
plot(t(idx), yd(idx),'LineWidth',1.2);
grid on; title('LPF Output (Demodulated Approximation)');
xlabel('t (s)'); ylabel('y_d(t)');

% 7
x_demod_fm = fmdemod(xc_manual, fc, fs, fDelta);
twin = 0.1;
idx = t <= twin;
x_seg      = x(idx);
yd_seg     = yd(idx);
fm_seg     = x_demod_fm(idx);
x_norm      = x_seg      / max(abs(x_seg));
yd_norm     = yd_seg     / max(abs(yd_seg));
fm_norm     = fm_seg     / max(abs(fm_seg));

figure;
plot(t(idx), x_norm, 'k','LineWidth',1.4); hold on;
plot(t(idx), yd_norm,'r','LineWidth',1.2);
plot(t(idx), fm_norm,'b--','LineWidth',1.2);
grid on;
xlabel('t (s)');
title('Comparison of Demodulated Signals (First 0.1 s)');
legend('Original x(t)','Freq-Discriminator Output','fmdemod Output');


function [x, x_int] = generate_signal(t)

u = @(s) double(s>=0);

x =  4*u(t) ...
    -7*u(t-4) ...
    +3*u(t-5) ...
    -3*u(t-8);
x_int = cumtrapz(t, x);

end

function [f, mag] = compute_fft(x, fs)

N = length(x);
X = fftshift(fft(x));
mag = abs(X)/N;

f = linspace(-fs/2, fs/2, N);
idx = f>=0;
f = f(idx);
mag = mag(idx);

end

function [maxDev, beta, fm_eff, BW] = carson_bandwidth(x, fDelta, fs)

maxDev = fDelta * max(abs(x));   % Δf_max
% -3 dB cutoff 
[f, Xmag] = compute_fft(x, fs);
Xnorm = Xmag / max(Xmag);
idx = find(Xnorm >= 0.707, 1, 'last');
fm_eff = f(idx);                  

% Modulation index β = Δf / fm
beta = maxDev / fm_eff;

% Carson Bandwidth rule
BW = 2*(maxDev + fm_eff);

end
