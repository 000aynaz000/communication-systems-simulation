
% question1
clear; close all; clc;
% (a) Create the probe signal 
fs = 8000;            
duration = 1.0;       
t = (0:1/fs:duration-1/fs).';   
fstart = 100;         
fend = 4000;    

% linear chirp
x = chirp(t, 0, duration, fend, 'linear'); %  freq sweep linear from 0->fend
k = (fend - fstart) / duration;                 % chirp rate
phi = 2*pi*(fstart*t + 0.5*k*t.^2);             % instantaneous phase
x = cos(phi);
N = length(x);
figure;
plot(t, x , "b");
xlabel("Time (s)")
ylabel("Amplitude");
title("Chirp Signal");
xlim([0 .7]);

% (b) Autocorrelation Rxx[k] 
[Rxx, lags_Rxx] = xcorr(x, 'biased');  
% Plot Rxx vs lag k (samples)
figure;
plot(lags_Rxx, Rxx);
xlabel('Lag k (samples)');
ylabel('R_{xx}[k]');
title('Autocorrelation R_{xx}[k] of the linear chirp');
grid on;
xlim([-1000 1000]);  % zoom near zero 

idx0 = find(lags_Rxx==0);
fprintf('R_{xx}[0] = %.6g (large spike at zero lag)\n', Rxx(idx0));

% (c) Load room data if provided; 
y_provided = false;
if isfile('room_data.mat')
    s = load('room_data.mat');
    if isfield(s,'y') && isfield(s,'x')
        y = s.y(:);
        x = s.x(:);   
        N = length(x);
        y_provided = true;
        fprintf('Loaded room_data.mat (using provided x and y)\n');
    elseif isfield(s,'y')
        y = s.y(:);
        fprintf('Loaded room_data.mat (using provided y; using generated x)\n');
        y_provided = true;
    end
elseif isfile('room_data.npz')
    try
        s = load('room_data.npz','-mat');  
        if isfield(s,'y'), y = s.y(:); y_provided = true; end
    catch
        warning('room_data.npz exists but could not be loaded in MATLAB automatically. Convert to .mat or reupload.');
    end
end
if ~exist('y','var')
    hlen = 2001;
    h = zeros(hlen,1);
    delays = [0, 800, 1600];     % sample delays (example)
    amps   = [1.0, 0.5, 0.3];    % amplitudes
    for ii = 1:length(delays)
        d = delays(ii)+1;        
        h(d) = amps(ii);
    end
    y = conv(x, h);
    y = y + 0.05*randn(size(y));  % add small noise
    fprintf('No room_data file found: running synthetic demo with echoes at %s samples\n', mat2str(delays));
end

% (d) Time-domain plot of y[n] 
Ny = length(y);
n_y = 0:(Ny-1);
figure;
plot(n_y, y);
xlabel('n (samples)');
ylabel('Amplitude');
title('Received signal y[n] (time domain)');
xlim([0 min(Ny,4000)]); % zoom first portion for clarity
grid on;

% (e) Cross-correlation Ryx[k] = sum_m y[m] x[m-k]
[Ryx, lags_yx] = xcorr(y, x, 'none');  
% lags_yx index: lag k = lags_yx(i)
Ryx = Ryx / max(abs(Ryx));

%(f) Plot Ryx vs lag k 
figure;
plot(lags_yx, Ryx);
xlabel('Lag k (samples)');
ylabel('R_{yx}[k]');
title('Cross-correlation R_{yx}[k] (y with x)');
grid on;

% zoom to expected positive lags
xlim([-500 min(length(Ryx)-1, 5000)]);
mask = lags_yx >= 0;
lags_pos = lags_yx(mask);
Ryx_pos = Ryx(mask);
pthresh = 0.05 * max(abs(Ryx_pos));   % 5% of max
[pk_vals, pk_locs_rel] = findpeaks(abs(Ryx_pos), 'MinPeakHeight', pthresh, 'SortStr', 'descend');
if isempty(pk_vals)
    warning('No peaks detected. Try lowering MinPeakHeight or smoothing Ryx.');
else
    pk_indices = pk_locs_rel(1:min(10,length(pk_locs_rel)));  
    pk_lags = lags_pos(pk_indices);
    pk_heights = Ryx_pos(pk_indices); 
    [~, order] = sort(abs(pk_heights),'descend');
    order = order(1:min(3,length(order)));
    top_lags = pk_lags(order);
    top_vals = pk_heights(order);
    % normalize relative amplitudes
    rel_ampl = abs(top_vals) / max(abs(top_vals));
   
    % Display top-3
    fprintf('\nTop echoes found (synthetic/demo or from file):\n');
    num_echoes = min([length(top_lags), length(top_vals)]);
    % Trim to the same size
    top_lags = top_lags(1:num_echoes);
    top_vals = top_vals(1:num_echoes);
    rel_ampl = abs(top_vals) ./ max(abs(top_vals));
    T = table((1:num_echoes)', top_lags(:), top_lags(:)/fs, top_vals(:), rel_ampl(:), ...
        'VariableNames', {'Rank', 'Lag_samples', 'Lag_seconds', 'Amplitude', 'Relative_Amplitude'});
    disp(T);
end

% (g) Conclusion text printed to console 
fprintf('\nConclusion:\n');
fprintf('Because R_{xx}[k] is strongly peaked at k = 0 (i.e. acts like a delta),\n');
fprintf('R_{yx}[k] = R_{hx}[k] * R_{xx}[k] approximately equals h[k] (scaled).\n');
fprintf('Therefore the peaks in R_{yx}[k] correspond to the impulse response (echoes) of the room.\n');
% Save plots for report 
savefig(gcf, 'Ryx_plot.fig');     % last figure (Ryx)
h1 = figure(1);
exportgraphics(h1, 'Rxx_autocorr.pdf','ContentType','vector');
h2 = figure(2);
exportgraphics(h2, 'x_signal.pdf','ContentType','vector');
h3 = figure(3);
exportgraphics(h3, 'y_signal.pdf','ContentType','vector');
h4 = figure(4);
exportgraphics(h4, 'Ryx_signal.pdf','ContentType','vector');
