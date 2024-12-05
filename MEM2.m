
clear; clc;

% Parameters for MEM and phase retrieval
max_iterations = 200 ; % Number of iterations
%beta = 0.9;            % Feedback parameter for MEM regularization
%lambda = 0.01;         % Lagrange multiplier for entropy contribution
N = 256; % Number of pixels
p = 0.1; % Pause time
beta = 0.1;            % Feedback parameter for MEM regularization
lambda = 0.01; 

% reading diffraction pattern
  fid = fopen('a_dp.bin', 'r');
  dp = fread(fid, [N, N], 'real*4');
  fclose(fid);   
  imagesc(rot90(log(dp + 10))), colorbar;
  colormap(gray);
  dp_amplitude = sqrt(dp);

% amplitude distribution in the sample plane
    fid = fopen('a_sample_amplitude.bin', 'r');
    sample_amplitude = fread(fid, [N, N], 'real*4');
    fclose(fid); 


% reading in the true (simulated) phase of the sample
    fid_phase = fopen('a_sample_phase.bin', 'r');
    true_phase = fread(fid_phase, [N, N], 'real*4');   
    fclose(fid_phase);

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
field_detector = dp_amplitude .* exp(i*phase);

% Start the iterative phase retrieval with MEM
for ii = 1:max_iterations 
    fprintf('Iteration: %d\n', ii)

 % updating complex-valued sample distribution
      sample_phase_updated = angle(IFT2Dc(field_detector));
      sample_amplitude_updated = abs(IFT2Dc(field_detector));
      sample_phase_updated = sample_phase_updated.*sample_amplitude_updated;
  
      % showing updated phase distribution of the sample
      subplot(1,2,1);
      %imshow(flipud(rot90(sample_phase_updated)), []);
      %imshow(flipud(rot90(sample_phase_updated)), [0, max(sample_phase_updated(:))]);
      sample_phase_shifted = sample_phase_updated-min(sample_phase_updated);
      imshow(flipud(rot90(sample_phase_shifted)),[]);
      title('MEM Updated phase distribution of the sample')
      xlabel({'x / px'})
      ylabel({'y / px'})
      axis on
      set(gca,'YDir','normal')
      colormap('gray')
      colorbar;

       % showing updated amplitude distribution of the sample
      subplot(1,2,2);
      imshow(flipud(rot90(sample_amplitude_updated)), []);
      title('Sample Amplitude')
      xlabel({'x'})
      ylabel({'y'})
      axis on
      set(gca,'YDir','normal')
      colormap('gray')
      colorbar;
      
      pause(p);
      sample_updated = sample_amplitude.*exp(i*sample_phase_updated);

    % Enforce object domain constraints with MEM
    % Calculate entropy term (maximize uniformity)
    magnitude = abs(sample_updated);
    magnitude_normalized = magnitude / sum(magnitude(:)); % Normalize
    entropy_term = -magnitude_normalized .* log(magnitude_normalized + eps); % Avoid log(0)
    entropy_term(isnan(entropy_term)) = 0; % Handle NaNs
    
    % Update the object with MEM regularization
    sample_updated_MEM=sample_updated .* (1 - beta) + ...
             beta * exp(lambda * sum(entropy_term(:))) .* abs(sample_updated);
      
     % updating complex-valued wavefront distribution in the detector plane
      field_detector_updated = FT2Dc(sample_updated_MEM);

    % replacing updated amplitude for the measured amplitude
      field_detector = dp_amplitude.*exp(i*angle(field_detector_updated));  

    % Calculate metrics
    % 1. Entropy of the object
    normalized_amplitude = sample_amplitude_updated / sum(sample_amplitude_updated(:));
    entropy_object(ii) = -sum(normalized_amplitude(:) .* log(normalized_amplitude(:) + eps));
    
    % 2. RMSE between true and reconstructed phase
    phase_difference = true_phase - sample_phase_updated;
    phase_residuals = angle(exp(1i * phase_difference)); % Wrap residuals to [-pi, pi]
    phase_rmse(ii) = sqrt(mean(phase_residuals(:).^2)); % RMSE of phase residuals
   
end

% showing reconstructed phase of the sample
      figure();
      subplot(1,3,1);
      imshow(flipud(rot90(sample_phase_updated)), [0, max(sample_phase_updated(:))]);
      title('(MEM) Sample Phase Distribution at Iteration: ',num2str(max_iterations));
      xlabel({'x / px'})
      ylabel({'y / px'})
      axis on
      set(gca,'YDir','normal')
      colormap('gray')
      colorbar; 
% showing reconstructed amplitude of the sample
      subplot(1,3,2);
      imshow(flipud(rot90(sample_amplitude_updated)), []);
      title('Sample Amplitude Distribution')
      xlabel({'x'})
      ylabel({'y'})
      axis on
      set(gca,'YDir','normal')
      colormap('gray')
      colorbar; 
      
  % showing the real phase distribution of the sample
      subplot(1,3,3);
      imshow(flipud(rot90(true_phase)), []);
      title('True Phase Distribution')
      xlabel({'x / px'})
      ylabel({'y / px'})
      axis on
      set(gca,'YDir','normal')
      colormap('gray')
      colorbar;

  % Entropy Plot
    figure();
    plot(1:max_iterations, entropy_object, 'b-', ...
         1:max_iterations, baseline_entropy_array, 'r--', 'LineWidth', 2);
    legend('Entropy (Object)', 'Baseline Entropy');
    title('Entropy of the Object vs Iteration');
    xlabel('Iteration');
    ylabel('Entropy');
    grid on;

 % Log-Scale Error Plot
    figure();
    semilogy(1:max_iterations, phase_rmse, 'LineWidth', 2);
    grid on;
    title('MEM Log-Scale Error: RMSE of Phase vs Iteration');
    xlabel('Iteration');
    ylabel('RMSE (Log Scale)');

      