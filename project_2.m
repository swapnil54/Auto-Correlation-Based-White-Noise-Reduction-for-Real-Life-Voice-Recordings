clc;
clear all;
close all;

% Load the clean voice recording
[voice, Fs] = audioread("C:\Users\moham\Desktop\matlab 2nd assignment\voice.wav");

% Ensure the signal is mono
if size(voice, 2) > 1
    voice = mean(voice, 2); % Convert to mono
end

% Fourier Transform to observe frequency components
L = length(voice);          % Length of the signal
Y = fft(voice);             % Compute Fourier Transform
P = abs(Y / L);             % Normalize the magnitude
f = Fs * (0:(L/2)) / L;     % Frequency vector (positive frequencies only)
P = P(1:L/2+1);             % Single-sided spectrum

% Plot the frequency components
figure;
plot(f, P);
title('Frequency Components of Clean Voice Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Optional: Play the clean voice
sound(voice, Fs);

%% 3.2 Preprocessing the Noisy Signal 
% Load the noisy voice recording
[noisy_voice, Fs] = audioread("C:\Users\moham\Desktop\matlab 2nd assignment\noise with voice.wav");

% Ensure the noisy signal is mono
if size(noisy_voice, 2) > 1
    noisy_voice = mean(noisy_voice, 2); % Convert to mono
end

% Normalize the noisy signal to the range [-1, 1]
noisy_voice = noisy_voice / max(abs(noisy_voice));

% Plot the normalized noisy signal
figure;
plot(noisy_voice);
title('Normalized Noisy Voice Signal');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;

%% 3.3 Auto-Correlation Calculation 
% Limit the length of the signal to avoid excessive memory usage
duration = 10; % seconds to process
segment_length = min(length(noisy_voice), Fs * duration);
noisy_voice_segment = noisy_voice(1:segment_length);

% Limit the lag range
max_lag = Fs * 0.1; % Analyze lags up to 0.1 seconds (100ms)
[R, lags] = xcorr(noisy_voice_segment, max_lag, 'coeff'); % Normalized auto-correlation

% Ensure sizes of lags and R match
if length(lags) ~= length(R)
    error('Mismatch between lags and R sizes. Check the xcorr output.');
end

% Plot the Auto-Correlation Function
figure;
plot(lags, R, 'b', 'LineWidth', 1.5);
title('Auto-Correlation of Noisy Signal');
xlabel('Lag (samples)');
ylabel('Correlation Coefficient');
grid on;

% Highlight peaks in ACF
[maxPeakValue, maxPeakIndex] = max(R); % Main peak at lag = 0
threshold = 0.05; % Example threshold
secondaryPeaks = find(R > threshold & lags ~= 0); % Significant peaks

% Annotate the peaks
hold on;
plot(lags(maxPeakIndex), maxPeakValue, 'ro', 'MarkerSize', 10); % Main peak

% Only plot secondary peaks if they exist
if ~isempty(secondaryPeaks)
    validSecondaryPeaks = secondaryPeaks(secondaryPeaks <= length(lags)); % Avoid indexing errors
    plot(lags(validSecondaryPeaks), R(validSecondaryPeaks), 'g*'); % Secondary peaks
else
    disp('No secondary peaks found in the Auto-Correlation Function.');
end
legend('Auto-Correlation', 'Primary Peak', 'Secondary Peaks');

%% 3.4 Noise Estimation and Thresholding 
% Thresholding to isolate periodic components (speech vs noise)
significant_peaks = R > threshold; % Identify significant peaks
noise_estimation = R .* (~significant_peaks); % Suppress noise components
speech_components = R .* significant_peaks; % Retain periodic components

% Plot the thresholded ACF
figure;
plot(lags, R, 'b', 'LineWidth', 1.5);
hold on;
if any(significant_peaks)
    plot(lags(significant_peaks), speech_components(significant_peaks), 'r.', 'MarkerSize', 10); % Speech components
    plot(lags(~significant_peaks), noise_estimation(~significant_peaks), 'g.', 'MarkerSize', 10); % Noise components
end
title('Thresholding in Auto-Correlation Function');
xlabel('Lag');
ylabel('Correlation Coefficient');
legend('ACF', 'Significant Peaks (Speech)', 'Suppressed Noise');

%% 3.5 Signal Reconstruction 
% Reconstruct the clean speech signal
reconstructed_signal = ifft(fft(noisy_voice) .* fft(speech_components, length(noisy_voice)));

% Normalize the reconstructed signal to the range [-1, 1]
reconstructed_signal = reconstructed_signal / max(abs(reconstructed_signal));

% Optional: Manually clip the signal to ensure no values exceed [-1, 1]
reconstructed_signal(reconstructed_signal > 1) = 1;
reconstructed_signal(reconstructed_signal < -1) = -1;

% Plot the original noisy signal and reconstructed signal
figure;
subplot(2, 1, 1);
plot(noisy_voice);
title('Noisy Voice Signal');
xlabel('Sample Index');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(real(reconstructed_signal));
title('Reconstructed Clean Voice Signal');
xlabel('Sample Index');
ylabel('Amplitude');

% Save the reconstructed signal to a file
audiowrite("C:\Users\moham\Desktop\matlab 2nd assignment\reconstructed_voice.wav", real(reconstructed_signal), Fs);

%% 3.6 Performance Evaluation 
% Compute Signal-to-Noise Ratio (SNR)
original_clean_power = sum(voice.^2) / length(voice); % Power of the clean signal
noise_power = sum((noisy_voice - real(reconstructed_signal)).^2) / length(noisy_voice); % Power of the residual noise
SNR_before = 10 * log10(original_clean_power / (sum(noisy_voice.^2) / length(noisy_voice))); % Before noise removal
SNR_after = 10 * log10(original_clean_power / noise_power); % After noise removal

% Display SNR results
fprintf('SNR Before Noise Removal: %.2f dB\n', abs(SNR_before));
fprintf('SNR After Noise Removal: %.2f dB\n', abs(SNR_after));


% Visual Comparison of Signals
figure;
subplot(3, 1, 1);
plot(noisy_voice);
title('Noisy Voice Signal');
xlabel('Sample Index');
ylabel('Amplitude');

subplot(3, 1, 2);
plot(real(reconstructed_signal));
title('Reconstructed Clean Signal');
xlabel('Sample Index');
ylabel('Amplitude');

subplot(3, 1, 3);
spectrogram(real(reconstructed_signal), 256, [], [], Fs, 'yaxis');
title('Spectrogram of Reconstructed Signal');
