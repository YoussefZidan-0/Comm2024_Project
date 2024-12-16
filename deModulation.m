% Constants
offset = 0;
IF = 25e3;  % Intermediate frequency (25 kHz)
fs_interpolated = interp_factor * fs;  % Interpolated sampling rate
receivedSignal = sentsignal;  % Received FDM signal
Rf_enable = 0;  % Enable/Disable RF filter

% Relaxed Filter Constraints
A_pass = 2;
A_stop1 = 20;
A_stop2 = 20;

fprintf("\nSignal Processing Initialized.\n");
while true
    % User input for channel selection
    fprintf("\nChoose one of these channels (Enter 0 to Exit):\n");
    fprintf("1. Short_BBCArabic2\n");
    fprintf("2. Short_FM9090\n");
    fprintf("3. Short_QuranPalestine\n");
    fprintf("4. Short_RussianVoice\n");
    fprintf("5. Short_SkyNewsArabia\n");
    i = input("Choose: ");
    
    % Exit condition
    if i == 0
        fprintf("Exiting program...\n");
        break;
    elseif i < 1 || i > 5
        fprintf("Invalid input! Please choose a number between 1 and 5, or 0 to exit.\n");
        continue;
    end
    
    % Carrier frequency and bandwidth for the selected channel
    fc_mod = 100e3 + (i - 1) * 50e3;  % Carrier frequency for the selected channel
    bw = BWs(i);  % Bandwidth (fixed for all channels)
    
    % Define proper increasing frequencies and cap at Nyquist
    F_stop1 =  fc_mod - bw ;  % Stopband start
    F_pass1 = fc_mod - bw/2;                % Passband start
    F_pass2 = fc_mod + bw/2;                % Passband end
    F_stop2 =  fc_mod + bw ;  % Stopband end

    % RF Bandpass Filter
    if Rf_enable == 1
        if F_stop2 > fs_interpolated / 2
            error("Filter frequencies exceed Nyquist limit. Increase fs_interpolated.");
        end
        rf_bandpass_filter = designfilt('bandpassiir', ...
            'StopbandFrequency1', F_stop1, 'PassbandFrequency1', F_pass1, ...
            'PassbandFrequency2', F_pass2, 'StopbandFrequency2', F_stop2, ...
            'StopbandAttenuation1', A_stop1, 'PassbandRipple', A_pass, ...
            'StopbandAttenuation2', A_stop2, 'SampleRate', fs_interpolated);
        fprintf("RF Filter Order: %d\n", filtord(rf_bandpass_filter));
        filtered_signal = filter(rf_bandpass_filter, receivedSignal);
    else
        filtered_signal = receivedSignal;
    end
    
    % Plot RF Stage Spectrum
    figure(1); clf;
    rf_spectrum = fftshift(fft(filtered_signal, 2^nextpow2(length(filtered_signal))));
    N = length(rf_spectrum);
    f = (-N/2:N/2-1) * fs_interpolated / N;
    subplot(2, 1, 1);
    plot(f, abs(rf_spectrum));
    title("RF Stage After BPF for Audio " + num2str(i));
    xlabel("Frequency (Hz)"); ylabel("Magnitude"); grid on;

    % Step 2: Mix with Local Oscillator to Shift to IF
    t = (1:length(filtered_signal))' / fs_interpolated;  % Time vector
    carrier_LO = cos(2 * pi * (fc_mod + IF) * t);  % Local oscillator signal
    mixed_signal_IF = filtered_signal .* carrier_LO;  % Mixer output

    % Plot IF Stage Spectrum
    IF_spectrum = fftshift(fft(mixed_signal_IF, 2^nextpow2(length(mixed_signal_IF))));
    subplot(2, 1, 2);
    plot(f, abs(IF_spectrum));
    title("IF Stage for Audio " + num2str(i));
    xlabel("Frequency (Hz)"); ylabel("Magnitude"); grid on;

    % Step 3: Bandpass Filter to Isolate IF
    F_stop1_IF = max(0, IF - bw - 5e3);  % Stopband start
    F_pass1_IF = IF - bw;                % Passband start
    F_pass2_IF = IF + bw;                % Passband end
    F_stop2_IF = min(fs_interpolated / 2 - 1, IF + bw + 5e3);  % Stopband end
    
    if_bandpass_filter = designfilt('bandpassiir', ...
        'StopbandFrequency1', F_stop1_IF, 'PassbandFrequency1', F_pass1_IF, ...
        'PassbandFrequency2', F_pass2_IF, 'StopbandFrequency2', F_stop2_IF, ...
        'StopbandAttenuation1', A_stop1, 'PassbandRipple', A_pass, ...
        'StopbandAttenuation2', A_stop2, 'SampleRate', fs_interpolated);
    fprintf("IF Filter Order: %d\n", filtord(if_bandpass_filter));
    filtered_signal_IF = filter(if_bandpass_filter, mixed_signal_IF);

    % Step 4: Demodulate IF Signal
    carrier_IF = cos(2 * pi * IF * t);  % IF carrier signal
    baseband_signal = filtered_signal_IF .* carrier_IF;  % Demodulation

    % Step 5: Apply Low-Pass Filter to Recover Baseband Signal
    F_pass = bw; F_stop = bw + 5e3;
    lowpass_filter = designfilt('lowpassiir', ...
        'PassbandFrequency', F_pass, 'StopbandFrequency', F_stop, ...
        'PassbandRipple', A_pass, 'StopbandAttenuation', A_stop1, ...
        'SampleRate', fs_interpolated);
    Base_Band_received_signal_LPF = filter(lowpass_filter, baseband_signal);

    % Resample and Playback the Demodulated Signal
    Base_Band_received_signal_LPF = 4 * resample(Base_Band_received_signal_LPF, 1, interp_factor);
    sound(Base_Band_received_signal_LPF, fs);

    % Plot Baseband Spectrum After LPF
    figure(2); clf;
    baseband_spectrum = fftshift(fft(Base_Band_received_signal_LPF, 2^nextpow2(length(Base_Band_received_signal_LPF))));
    N = length(baseband_spectrum);
    f = (-N/2:N/2-1) * fs_interpolated / N;
    plot(f, abs(baseband_spectrum));
    title("Baseband Spectrum After LPF for Audio " + num2str(i));
    xlabel("Frequency (Hz)"); ylabel("Magnitude"); grid on;

    fprintf("Audio %d processed successfully.\n", i);
end
fprintf("Signal Processing Completed.\n");
