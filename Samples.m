%%section1: load data

[c1,fs1]=audioread("Short_BBCArabic2.wav");
[c2,fs2]=audioread("Short_FM9090.wav");
[c3,fs3]=audioread("Short_QuranPalestine.wav");
[c4,fs4]=audioread("Short_RussianVoice.wav");
[c5,fs5]=audioread("Short_SkyNewsArabia.wav");
audioSignals = {c1, c2, c3, c4, c5};
samplingRates = [fs1, fs2, fs3, fs4, fs5];
Left_channel = cell(1, length(audioSignals));
Right_channel = cell(1, length(audioSignals));
monoSignals = cell(1, 5);
paddedSignals = cell(1, 5);
maxLength = 0;

for i=1:length(audioSignals)
    Left_channel{i}=audioSignals{i}(:,1);
    Right_channel{i} = audioSignals{i}(:, 2);
    monoSignals{i}=sum(audioSignals{i}, 2);
    maxLength = max(maxLength, length( monoSignals{i}));

end
% sound(monoSignals{2},fs2)

% Pad signals to the maximum length
for i = 1:length(monoSignals)
    x = monoSignals{i};
    % Pad the signal with zeros to match the maximum length
    paddedSignals{i} = [x; zeros(maxLength - length(x), 1)];
end
%% section 2 : plot and extract data

% plotting for the all original signal
BWs = zeros(1,5);
for i = 1:length(audioSignals)
    x = audioSignals{i};
    fs = samplingRates(i);

    % FFT
    N = length(x);
    X = fft(x, N);
    f = (-N/2:N/2-1) * fs / N;

    %[pxx, f_pxx] = periodogram(x, [], length(x), fs);
    %  avg_pxx = mean(pxx, 2);

    %     pxx=abs(x).^2;
    %     tot_sum=cumsum(pxx);
    %     normlizedsum=tot_sum/tot_sum(end);

    BWs(i)=obw(monoSignals{i},fs);

    % Plot
    grid on
    figure;
    plot(f, abs(fftshift(X)) / N);
    title(['Frequency Spectrum of Signal ', num2str(i)]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
end
grid off


BWs_Hz = BWs ;
BWs_kHz = BWs_Hz / 1000;
