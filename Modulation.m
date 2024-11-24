fc_mod = zeros(1,6);  % Carrier frequencies
t = cell(1, 5);  % Time vectors for each signal
modulated_signal = cell(1, 5);  % Modulated signals
sentsignal=[];
for i = 1:length(audioSignals)
    % Retrieve the current mono signal
    x_in = monoSignals{i};
    fs = samplingRates(i);  % Sampling rate for this signal

    % Take quarter of the length of the signal
    quarter_length = floor(length(x_in) / 4);
    x = x_in(1:quarter_length);  % Use the first quarter of the signal

    % Time vector for the current signal
    t{i} = (0:length(x)-1) / fs;  % Time vector, matching the length of the signal

    % Carrier frequency (in Hz)
    fc_mod(i) = 100e3 + (i - 1) * 50e3;  % Carrier frequency for the i-th signal

    % Check if the carrier frequency exceeds Nyquist frequency
    if fc_mod(i) > fs / 2
        % Interpolate both the signal and the carrier
        interp_factor = 5;  % Increase sampling rate by a factor of 5
        x_inter = interp1(1:length(x), x, linspace(1, length(x), length(x) * interp_factor), 'linear');
        t_interpolated = (0:length(x_inter)-1) / (interp_factor * fs);  % New time vector
        Mod_carrier = cos(2 * pi * fc_mod(i) * t_interpolated);  % Interpolated carrier signal
        modulated_signal{i} = x_inter .* Mod_carrier;  % Modulation with the interpolated carrier
    else
        % No interpolation, just modulate
        Mod_carrier = cos(2 * pi * fc_mod(i) * t{i});  % Carrier signal
        modulated_signal{i} = x .* Mod_carrier;  % Modulation
    end
end
for i = 1:length(modulated_signal)
    % Ensure the signal is a column vector before padding
    modulated_signal{i} = modulated_signal{i}(:);  % Convert to column vector if not already
    
    % Pad the shorter signals with zeros
    modulated_signal{i} = [modulated_signal{i}; zeros(max_length - length(modulated_signal{i}), 1)];
    
    % Combine into the FDM signal
    if isempty(sentsignal)
        sentsignal = modulated_signal{i};
    else
        sentsignal = sentsignal + modulated_signal{i};
    end
end

% Plot the combined FDM signal in the time domain
figure;
subplot(2, 1, 1);
t_total = (0:length(sentsignal)-1) / fs;  % Time vector for the combined signal
plot(t_total, sentsignal);
title('FDM Signal (Time Domain)');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot the combined FDM signal in the frequency domain
subplot(2, 1, 2);
N = length(sentsignal);  % Number of samples in the combined signal
f = (-N/2:N/2-1) * fs / N;  % Frequency axis
spectrum_fdm = fftshift(fft(sentsignal, N));  % Compute and shift the FFT
plot(f, abs(spectrum_fdm) / N);
title('FDM Signal (Frequency Domain)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
%sound(sentsignals,fs)