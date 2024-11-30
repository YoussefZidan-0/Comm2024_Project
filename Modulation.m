<<<<<<< HEAD
fc_mod = zeros(1,5);  % Carrier frequencies
t = cell(1, 5);  % Time vectors for each signal
modulated_signal = cell(1, 5);  % Modulated signals
sentsignal=[];
max_length=0;
interp_factor = 14;
fs=44100;
fs_interpolated=interp_factor*fs;
for i = 1:length(audioSignals)
    % Retrieve the current mono signal
    x = paddedSignals{i};

    % Take quarter of the length of the signal
    half_length = floor(length(x) / 1);
    x = x(1:half_length);  % Use the first quarter of the signal

    % Carrier frequency (in Hz)
    fc_mod(i) = 100e3 + (i - 1) * 50e3;  % Carrier frequency for the i-th signal



        x_inter=interp(x,interp_factor);
                 t = (0:1/fs_interpolated:((size(x_inter, 1)-1)/fs_interpolated))';
                 t = t(:, ones(1, size(x_inter, 2)));
        Mod_carrier = cos(2 * pi * fc_mod(i) * t);  
        modulated_signal{i} = x_inter .* Mod_carrier;  
        
        %modulated_signal{i}=ammod(x_inter,fc_mod(i),fs*interp_factor);
           %     t_interpolated = (0:length(x_inter)-1) / (interp_factor * fs);  % New time vector

        max_length = max(max_length, length(modulated_signal{i}));

end

for i = 1:length(modulated_signal)
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
t_total = (0:length(sentsignal)-1) / (10*fs);  % Time vector for the combined signal
plot(t_total, sentsignal);
title('FDM Signal (Time Domain)');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot the combined FDM signal in the frequency domain
subplot(2, 1, 2);
N = length(sentsignal);  % Number of samples in the combined signal
f = (-N/2:N/2-1) * 10*fs / N;  % Frequency axis
spectrum_fdm = fftshift(fft(sentsignal, N));  % Compute and shift the FFT
plot(f, abs(spectrum_fdm) / N);
title('FDM Signal (Frequency Domain)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
%sound(sentsignal, 10 * fs);
=======
fc_mod = zeros(1,5);  % Carrier frequencies
t = cell(1, 5);  % Time vectors for each signal
modulated_signal = cell(1, 5);  % Modulated signals
sentsignal=[];
max_length=0;
interp_factor = 14;
fs=44100;
fs_interpolated=interp_factor*fs;
for i = 1:length(audioSignals)
    % Retrieve the current mono signal
    x = paddedSignals{i};

    % Take quarter of the length of the signal
    half_length = floor(length(x) / 1);
    x = x(1:half_length);  % Use the first quarter of the signal

    % Carrier frequency (in Hz)
    fc_mod(i) = 100e3 + (i - 1) * 50e3;  % Carrier frequency for the i-th signal



        x_inter=interp(x,interp_factor);
                 t = (0:1/fs_interpolated:((size(x_inter, 1)-1)/fs_interpolated))';
                 t = t(:, ones(1, size(x_inter, 2)));
        Mod_carrier = cos(2 * pi * fc_mod(i) * t);  
        modulated_signal{i} = x_inter .* Mod_carrier;  
        
        %modulated_signal{i}=ammod(x_inter,fc_mod(i),fs*interp_factor);
           %     t_interpolated = (0:length(x_inter)-1) / (interp_factor * fs);  % New time vector

        max_length = max(max_length, length(modulated_signal{i}));

end

for i = 1:length(modulated_signal)
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
t_total = (0:length(sentsignal)-1) / (10*fs);  % Time vector for the combined signal
plot(t_total, sentsignal);
title('FDM Signal (Time Domain)');
xlabel('Time (s)');
ylabel('Amplitude');

% Plot the combined FDM signal in the frequency domain
subplot(2, 1, 2);
N = length(sentsignal);  % Number of samples in the combined signal
f = (-N/2:N/2-1) * 10*fs / N;  % Frequency axis
spectrum_fdm = fftshift(fft(sentsignal, N));  % Compute and shift the FFT
plot(f, abs(spectrum_fdm) / N);
title('FDM Signal (Frequency Domain)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
%sound(sentsignal, 10 * fs);
>>>>>>> bdfa72de2de8dc9a03a901e43a42f503bc13dc5a
