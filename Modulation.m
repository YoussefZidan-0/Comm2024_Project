% Constants
fc_mod = zeros(1, 5);  % Carrier frequencies
modulated_signal = cell(1, 5);  % Modulated signals
sentsignal = [];  % Combined FDM signal
interp_factor = 14;  % Interpolation factor
fs = 44100;  % Sampling frequency
fs_interpolated = interp_factor * fs;  % Interpolated sampling frequency

% Modulation process
for i = 1:length(audioSignals)
    % Retrieve the padded mono signal
    x = paddedSignals{i};

    % Take the first portion of the signal
    x = x(1:floor(length(x) / 1));  

    % Carrier frequency (in Hz)
    fc_mod(i) = 100e3 + (i - 1) * 50e3;

    % Interpolate the signal
    x_inter = interp(x, interp_factor);
    
    % Time vector for the interpolated signal
    t = (0:1/fs_interpolated:(length(x_inter) - 1) / fs_interpolated)';
    
    % Generate the carrier signal and perform modulation
    Mod_carrier = cos(2 * pi * fc_mod(i) * t);
    modulated_signal{i} = x_inter .* Mod_carrier;
end

% Combine all modulated signals into the FDM signal
for i = 1:length(modulated_signal)
    if isempty(sentsignal)
        sentsignal = modulated_signal{i};
    else
        sentsignal = sentsignal + modulated_signal{i};
    end
end

% Plot the combined FDM signal in the time domain
figure;
subplot(2, 1, 1);
t_total = (0:length(sentsignal) - 1) / fs_interpolated;  % Time vector for FDM signal
plot(t_total, sentsignal);
title('FDM Signal (Time Domain)');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot the combined FDM signal in the frequency domain
subplot(2, 1, 2);
N = length(sentsignal);
f = (-N/2:N/2-1) * fs_interpolated / N;  % Frequency axis
spectrum_fdm = fftshift(fft(sentsignal, N));  % Compute and shift the FFT
plot(f, abs(spectrum_fdm) / N);
title('FDM Signal (Frequency Domain)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

% Play the FDM signal
%sound(sentsignal, fs_interpolated);
clearvars -except BWs   sentsignal interp_factor fs