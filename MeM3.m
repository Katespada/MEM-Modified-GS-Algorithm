clear; clc; close all;

% Parameters for MEM and phase retrieval
max_iterations = 1e5 ; % Number of iterations
p = 0.1; % Pause time
beta = 0.1;            % Feedback parameter for MEM regularization
lambda = 0.01; 

[sample_amplitude, true_phase, dp_amplitude] = makeDiffractionPattern;
sample_amplitude = rot90(rot90(sample_amplitude));
N = size(dp_amplitude,2);

% Calculating the baseline entropy of the true phase profile
% Normalize the pixel intensities
normalized_phase = true_phase / sum(true_phase(:));
% Avoid numerical issues with log(0) by adding a small epsilon
epsilon = 1e-12;
normalized_phase = normalized_phase + epsilon;
% Calculate Shannon entropy
baseline_entropy = -sum(normalized_phase(:) .* log(normalized_phase(:)));
baseline_entropy_array = baseline_entropy * ones(1,max_iterations);

%Initialize metrics
    entropy_object = zeros(max_iterations, 1); % Entropy of the object
    phase_rmse = zeros(max_iterations, 1);    % RMSE between true and reconstructed phase


% Guess for the complex-valued wave distribution
% distribution
phase = (2*rand(N,N) - 1)*pi;
FT2 = dp_amplitude .* exp(1i*phase);

% Start the iterative phase retrieval with MEM
for ii = 1:max_iterations 
    fprintf('Iteration: %d\n', ii)

    if strcmp(get(gcf,'currentCharacter'),char(8)) % hit the backspace key to stop
        break;
    end

    % updating complex-valued sample distribution
    IFT2_prime = IFT2Dc(FT2);

    % replace amplitude keep phase
    IFT2 = sample_amplitude .* exp(1i*angle(IFT2_prime));

    % transform to detector
    FT2_prime = FT2Dc(IFT2);

    % replace amplitude keep phase
    FT2 = dp_amplitude .* exp(1i*angle(FT2_prime));

    if mod(ii,10)==0
        tiledlayout(2,2)
        nexttile(1)
        imagesc(abs(IFT2_prime));   %amplitude before update
        axis image xy;

        nexttile(2)
        imagesc(angle(IFT2_prime)); %reconstructed phase
        axis image xy;

        nexttile(3)
        imagesc(sample_amplitude);   %true amplitude
        axis image xy;
        nexttile(4)
        imagesc(true_phase);         %true phase
        axis image xy;
        drawnow;

    end

end
