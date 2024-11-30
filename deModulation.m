<<<<<<< HEAD
% Constants
offset = 0;
IF = 25e3 ;  % Intermediate frequency (25 kHz)

fs_interpolated = interp_factor * fs;  % Sampling frequency after interpolation
%receivedSignal = awgn(sentsignal,10);  
receivedSignal = sentsignal; 

% Initialize arrays for demodulated signals
demodulated_signals = cell(1, length(audioSignals));

% User input for channel selection
fprintf("Choose one of these channels:\n");
fprintf("1. Short_BBCArabic2\n");
fprintf("2. Short_FM9090\n");
fprintf("3. Short_QuranPalestine\n");
fprintf("4. Short_RussianVoice\n");
fprintf("5. Short_SkyNewsArabia\n");
i = input("Choose: ");

% Carrier frequency for the given signal
fc_mod = 100e3 + (i - 1) * 50e3;  % Carrier frequency for the selected channel
bw = BWs(i);  % Bandwidth for the selected signal
if bw / 2 > 25000
    bw = 24999;  % Cap bandwidth to avoid exceeding Nyquist limit
end

% Filter specifications for RF stage
F_stop1 = fc_mod - bw - bw/4;
F_pass1 = fc_mod - bw;
F_pass2 = fc_mod + bw;
F_stop2 = fc_mod + bw + bw/4;
A_pass = 1; 
A_stop1 = 40; 
A_stop2 = 40;

% Step 1: Bandpass filter to isolate the modulated signal
rf_bandpass_filter = fdesign.bandpass(F_stop1, F_pass1, F_pass2, F_stop2, A_stop1, A_pass, A_stop2, fs_interpolated);
rf_bandpass_filter = design(rf_bandpass_filter, 'equiripple');
filtered_signal = filter(rf_bandpass_filter, receivedSignal);
%filtered_signal=receivedSignal;


% Plot RF stage after BPF
figure('Name', 'RF Stage: Filtered Signal');
rf_spectrum = fftshift(fft(filtered_signal));
N = length(rf_spectrum);
f = (-N/2:N/2-1) * fs_interpolated / N;
subplot(2, 1, 1);
plot(f, abs(rf_spectrum));
title("RF Stage After BPF for Audio " + num2str(i));
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

% Step 2: Mix with a local oscillator to shift the signal to IF
t = (1:length(filtered_signal))' / fs_interpolated;  % Time vector
carrier_LO = cos(2 * pi * (fc_mod + IF + offset) * t);  % Local oscillator signal
mixed_signal_IF = filtered_signal .* carrier_LO;  % Mixer output

% Plot IF stage
IF_spectrum = fftshift(fft(mixed_signal_IF));
subplot(2, 1, 2);
plot(f, abs(IF_spectrum));
title("IF Stage for Audio " + num2str(i));
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

% Step 3: Bandpass filter the mixed signal to isolate IF
F_stop1_IF = IF - bw - bw/4;
F_pass1_IF = IF - bw;
F_pass2_IF = IF + bw;
F_stop2_IF = IF + bw + bw/4;
if_bandpass_filter = fdesign.bandpass(F_stop1_IF, F_pass1_IF, F_pass2_IF, F_stop2_IF, A_stop1, A_pass, A_stop2, fs_interpolated);
if_bandpass_filter = design(if_bandpass_filter, 'equiripple');
filtered_signal_IF = filter(if_bandpass_filter, mixed_signal_IF);

% Step 4: Demodulate the IF signal by multiplying with an IF carrier
carrier_IF = cos(2 * pi * IF * t);  % IF carrier signal
baseband_signal = filtered_signal_IF .* carrier_IF;  % Demodulation

% Step 5: Low-pass filter to recover the baseband signal
F_pass = bw;
F_stop = bw + 5000;
lowpass_filter = fdesign.lowpass(F_pass, F_stop, A_pass, A_stop2, fs_interpolated);
lowpass_filter = design(lowpass_filter, 'butter');
Base_Band_received_signal_LPF = filter(lowpass_filter, baseband_signal);

% Resample and playback the demodulated signal
Base_Band_received_signal_LPF = 4 * resample(Base_Band_received_signal_LPF, 1, interp_factor);
sound(Base_Band_received_signal_LPF, fs);

% Plot baseband spectrum after LPF
figure('Name', 'Baseband Signal After LPF');
baseband_spectrum = fftshift(fft(Base_Band_received_signal_LPF));
N = length(baseband_spectrum);
f = (-N/2:N/2-1) * fs_interpolated / N;
subplot(2, 1, 2);
plot(f, abs(baseband_spectrum));
title("Baseband Spectrum After LPF for Audio " + num2str(i));
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;
=======
% Constants
offset = 0;
IF = 25e3 + offset;  % Intermediate frequency (25 kHz)

fs_interpolated = interp_factor * fs;  % Sampling frequency after interpolation
receivedSignal = awgn(sentsignal,10);  

% Initialize arrays for demodulated signals
demodulated_signals = cell(1, length(audioSignals));

% User input for channel selection
fprintf("Choose one of these channels:\n");
fprintf("1. Short_BBCArabic2\n");
fprintf("2. Short_FM9090\n");
fprintf("3. Short_QuranPalestine\n");
fprintf("4. Short_RussianVoice\n");
fprintf("5. Short_SkyNewsArabia\n");
i = input("Choose: ");

% Carrier frequency for the given signal
fc_mod = 100e3 + (i - 1) * 50e3;  % Carrier frequency for the selected channel
bw = BWs(i);  % Bandwidth for the selected signal
if bw / 2 > 25000
    bw = 24999;  % Cap bandwidth to avoid exceeding Nyquist limit
end

% Filter specifications for RF stage
F_stop1 = fc_mod - bw - bw/4;
F_pass1 = fc_mod - bw;
F_pass2 = fc_mod + bw;
F_stop2 = fc_mod + bw + bw/4;
A_pass = 1; 
A_stop1 = 40; 
A_stop2 = 40;

% Step 1: Bandpass filter to isolate the modulated signal
rf_bandpass_filter = fdesign.bandpass(F_stop1, F_pass1, F_pass2, F_stop2, A_stop1, A_pass, A_stop2, fs_interpolated);
rf_bandpass_filter = design(rf_bandpass_filter, 'equiripple');
filtered_signal = filter(rf_bandpass_filter, receivedSignal);

% Plot RF stage after BPF
figure('Name', 'RF Stage: Filtered Signal');
rf_spectrum = fftshift(fft(filtered_signal));
N = length(rf_spectrum);
f = (-N/2:N/2-1) * fs_interpolated / N;
subplot(2, 1, 1);
plot(f, abs(rf_spectrum));
title("RF Stage After BPF for Audio " + num2str(i));
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

% Step 2: Mix with a local oscillator to shift the signal to IF
t = (1:length(filtered_signal))' / fs_interpolated;  % Time vector
carrier_LO = cos(2 * pi * (fc_mod + IF) * t);  % Local oscillator signal
mixed_signal_IF = filtered_signal .* carrier_LO;  % Mixer output

% Plot IF stage
IF_spectrum = fftshift(fft(mixed_signal_IF));
subplot(2, 1, 2);
plot(f, abs(IF_spectrum));
title("IF Stage for Audio " + num2str(i));
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;

% Step 3: Bandpass filter the mixed signal to isolate IF
F_stop1_IF = IF - bw - bw/4;
F_pass1_IF = IF - bw;
F_pass2_IF = IF + bw;
F_stop2_IF = IF + bw + bw/4;
if_bandpass_filter = fdesign.bandpass(F_stop1_IF, F_pass1_IF, F_pass2_IF, F_stop2_IF, A_stop1, A_pass, A_stop2, fs_interpolated);
if_bandpass_filter = design(if_bandpass_filter, 'equiripple');
filtered_signal_IF = filter(if_bandpass_filter, mixed_signal_IF);

% Step 4: Demodulate the IF signal by multiplying with an IF carrier
carrier_IF = cos(2 * pi * IF * t);  % IF carrier signal
baseband_signal = filtered_signal_IF .* carrier_IF;  % Demodulation

% Step 5: Low-pass filter to recover the baseband signal
F_pass = bw;
F_stop = bw + 5000;
lowpass_filter = fdesign.lowpass(F_pass, F_stop, A_pass, A_stop2, fs_interpolated);
lowpass_filter = design(lowpass_filter, 'butter');
Base_Band_received_signal_LPF = filter(lowpass_filter, baseband_signal);

% Resample and playback the demodulated signal
Base_Band_received_signal_LPF = 4 * resample(Base_Band_received_signal_LPF, 1, interp_factor);
sound(Base_Band_received_signal_LPF, fs);

% Plot baseband spectrum after LPF
figure('Name', 'Baseband Signal After LPF');
baseband_spectrum = fftshift(fft(Base_Band_received_signal_LPF));
N = length(baseband_spectrum);
f = (-N/2:N/2-1) * fs_interpolated / N;
subplot(2, 1, 2);
plot(f, abs(baseband_spectrum));
title("Baseband Spectrum After LPF for Audio " + num2str(i));
xlabel("Frequency (Hz)");
ylabel("Magnitude");
grid on;
>>>>>>> bdfa72de2de8dc9a03a901e43a42f503bc13dc5a
