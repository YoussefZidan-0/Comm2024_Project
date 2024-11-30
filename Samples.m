%% Section 1: Load Data
addpath 'Audio Signals'\
[c1, fs1] = audioread("Short_BBCArabic2.wav");
[c2, fs2] = audioread("Short_FM9090.wav");
[c3, fs3] = audioread("Short_QuranPalestine.wav");
[c4, fs4] = audioread("Short_RussianVoice.wav");
[c5, fs5] = audioread("Short_SkyNewsArabia.wav");

audioSignals = {c1, c2, c3, c4, c5};
samplingRates = [fs1, fs2, fs3, fs4, fs5];
Left_channel = cell(1, length(audioSignals));
Right_channel = cell(1, length(audioSignals));
monoSignals = cell(1, 5);
paddedSignals = cell(1, 5);
maxLength = 0;

% Extract left, right, and mono channels; calculate max signal length
for i = 1:length(audioSignals)
    Left_channel{i} = audioSignals{i}(:, 1);
    Right_channel{i} = audioSignals{i}(:, 2);
    monoSignals{i} = sum(audioSignals{i}, 2);  % Combine left and right channels
    maxLength = max(maxLength, length(monoSignals{i}));
end

% Pad signals to match the maximum length
for i = 1:length(monoSignals)
    paddedSignals{i} = [monoSignals{i}; zeros(maxLength - length(monoSignals{i}), 1)];
end

%% Section 2: Plot and Extract Data

BWs = zeros(1, 5);  % Initialize array for bandwidths

% Plot frequency spectra and calculate bandwidths
for i = 1:length(audioSignals)
    x = audioSignals{i};
    fs = samplingRates(i);

    % FFT for frequency spectrum
    N = length(x);
    X = fft(x, N);
    f = (-N/2:N/2-1) * fs / N;

    % Calculate bandwidth
    BWs(i) = obw(monoSignals{i}, fs);

    % Plot frequency spectrum
    figure;
    plot(f, abs(fftshift(X)) / N);
    title(['Frequency Spectrum of Signal ', num2str(i)]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on;
end

% Convert bandwidths to Hz and kHz
BWs_Hz = BWs;
BWs_kHz = BWs_Hz / 1000;
clearvars -except BWs paddedSignals audioSignals
