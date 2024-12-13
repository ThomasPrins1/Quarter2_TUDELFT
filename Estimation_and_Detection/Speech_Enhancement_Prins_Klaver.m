% Load Data.mat
load('Data.mat')
% Variables:

freq=16000; % Sample frequency
Data_length = 601655; % Total length of data given
FL = 0.020; % Frame length time based
L = freq*FL; % Frame length sample based
overlap = L/2; % Overlap samples with 50% overlap
i = 0; % Counter
f_k = freq/L; % in hz, frequency bin
nrframes = floor(Data_length/160); % nr of frames over all samples
K = 1024; % FFT length

result = zeros(Data_length,1); % Initialize results per frame
fft_result = zeros(Data_length,1);
Framed_matrix = zeros(Data_length,nrmics); % Intialize the framed and windowed matrix composed of results matrix
A = ones(nrmics,1); % Given matrix comprised of ones (weights)
R_w = zeros(nrmics,nrmics); % Covariance matrix
Variance = eye(nrmics,1); % Variance per microphone
S_lk = zeros(K,nrframes,nrmics);
S_hat_lk = zeros(K,nrframes,nrmics);
S_hat = zeros(K,nrframes,nrmics);
temp2 = zeros(nrmics,1);
var_emp = zeros(nrmics,1);
fft_data_slice = zeros(nrmics,1);

% Create hanning window
Window = hann(L);
% Find the covariance of the noise
R_w = cov(Data(1:freq,:));

%% exercise 1: Estimator
% Loop though microphones
    
    % This uses stft to frame and fft the data using a 320 hanning window:
    fft_framed_data = stft(Data,freq,Window=hann(L),OverlapLength=overlap,FFTLength=K);
    for l=1:L
        for k=1:K
            fft_data_slice = squeeze(fft_framed_data(k, l, :));  % Get a row vector for all microphones
            S_hat(k,l,:) = inv(A'*inv(R_w)*A)*A'*inv(R_w)*fft_data_slice;
        end
    end
    s_hat = ifft(S_hat);
%% OLD
    % Multiply frame with each sample
    % At indexes 1:161:321:481...
    %for frame=1:overlap+1:Data_length-L
    %    Framed_data(:,1) = Data(frame:frame+L-1,index);
    %    temp = Framed_data.*Window;
    %    result(frame:frame+L-1,1) = result(frame:frame+L-1,1) + temp;
    %end
    % Framed_matrix(:,index) = result;
    % Fft of DATA per microphone
    % BLUE estimator between microphones
    % S_hat = inv(A'*inv(R_w)*A)*A'*inv(R_w)*Framed_matrix';
    % question is why BLUE? This should be equal to MVU

    %% Exercise 2: emperical variance calculation:
    % fft is taken over all data with a hanning window:

    fft_data = stft(Data,freq,FFTLength=K);
    for index = 1:nrmics
        for l=1:L
            for k=1:K
                temp2(index,1) = temp2(index,1) + (abs(S_hat(k,l)-fft_data(k,l,index))^2);
            end
        end
        var_emp(index,1) = sum(temp2(1:index,1)/(K*L))/index;
    end
% Plot results
figure;
plot(1:nrmics, var_emp);
xlabel('Number of Microphones');
ylabel('Variance of Estimator');
title('Variance of Target Estimator as a Function of Microphones');
