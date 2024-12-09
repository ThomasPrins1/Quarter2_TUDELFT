% Load Data.mat
load('Data.mat')
% Variables:

freq=16000; % Sample frequency
Data_length = 601655; % Total length of data given
FL = 0.020; % Frame length time based
L = freq*FL; % Frame length sample based
overlap = L/2; % Overlap samples with 50% overlap
i = 0; % Counter
result = zeros(Data_length,1); % Initialize results per frame
Framed_matrix = zeros(Data_length,nrmics); % Intialize the framed and windowed matrix composed of results matrix
fft_Framed_matrix = zeros(Data_length,nrmics); % Fft of framed matrix
A = ones(nrmics,1); % Given matrix comprised of ones (weights)
R_w = zeros(nrmics,nrmics); % Covariance matrix
Variance = eye(nrmics,1); % Variance per microphone
% Create hanning window
Window = hann(L+1);

% Loop though microphones
for index=1:nrmics
    % Find variance for each microphone (not used currently):
    Variance(index,1) = var(Data(1:freq+1,index));
    % Multiply frame with each sample
    % At indexes 1:161:321:481...
    for frame=1:overlap:Data_length-L
        Framed_data(:,1) = Data(frame:frame+L,index);
        temp = Framed_data.*Window;
        result(frame:frame+L,1) = result(frame:frame+L,1) + temp;
        %i = i + 1; %frame number
    end
    % Calculate covariance matrix
    R_w = cov(Data(1:freq+1,:));
    Framed_matrix(:,index) = result;
    % Fft of DATA per microphone
    fft_Framed_matrix(:,index) = fft(result,Data_length);
    % BLUE estimator between microphones
    S_hat = inv(A'*inv(R_w)*A)*A'*inv(R_w)*Framed_matrix';
end
