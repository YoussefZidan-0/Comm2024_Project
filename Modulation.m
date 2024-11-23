fc_mod = zeros(1,6);
Omega = zeros(1,6);
t = cell(1, 5);
modulated_signal = cell(1,5);

for i = 1:length(audioSignals)
    % Retrieve the current right channel signal (mono signal assumed)
    x_in = monoSignals{i};
    fs = samplingRates(i);  % Sampling rate for this signal

    % Take quarter of the length of the signal
    half_length = floor(length(x_in) / 4);
    x = x_in(1:half_length);  % Use the first quarter of the signal

    % Time vector for the current signal
    t{i} = (0:length(x)-1) / fs;  % Time vector, matching the length of the signal

    % Carrier frequency (in Hz)

        fc_mod(i) = (100e3 + (i - 1) * 50e3);  

    % Check if the carrier frequency exceeds Nyquist frequency
    if fc_mod(i) * 2 > fs
        % Interpolate both the signal and the carrier
        interp_factor = 5;  % Increase sampling rate by a factor of 5
        x_interpolated = interp1(1:length(x), x, linspace(1, length(x), length(x) * interp_factor), 'linear');
        t_interpolated = linspace(0, length(x) / fs, length(x) * interp_factor);
        Mod_carrier = cos(2 * pi * fc_mod(i) * t_interpolated);  % Interpolated carrier signal
        modulated_signal{i} = x_interpolated .* Mod_carrier;  % Modulation with the interpolated carrier
    else
        % No interpolation, just modulate
        Mod_carrier = cos(2 * pi * fc_mod(i) * t{i});  % Carrier signal
        modulated_signal{i} = x .* Mod_carrier;  % Modulation
    end
end
