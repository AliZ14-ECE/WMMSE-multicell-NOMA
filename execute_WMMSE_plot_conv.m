% convergence plot

alpha_rng_length = 10; % number of alpha values to compute.
NC = 4; % #Cells which equals #BSs
NU = 10; % #USers in each cell.
P = 16;
nvar = 1.9905e-08; % Noise Variance 
epsilon = 1e-5; % For convergence test.
inner_radius = 500; 
minR_ratio = 0.01;
numIter = 2000;
num_reals = 1000; %# of channel realizations
alpha_rng = [1,alpha_rng_length]; % the range of the alpha computation indices
seed = 1;

%% system model
% load the presaved channel parameters.
file_name = sprintf('WMMSE_for_conv/WMMSE_%dx%dpower%dabs.mat', NC, NU,P);
if (exist(file_name, 'file'))
    load(file_name, "conv", "WR_vs_iter");
else
    channelsFileName = sprintf('channels_for_NU/Channels%dx%dpower%d.mat', NC, NU, P);
    load(channelsFileName,'H', 'in', 'D'); 
    executedFrom = '1';
    execute_WMMSE
end

plot_conv(alpha_rng, conv, WR_vs_iter);


function [] = plot_conv(alpha_rng, conv, WR_vs_iter)
  "plotting convergence"
  
  NIter = zeros(1, 10);
  S = zeros(1, 10);
  for alpha_idx = alpha_rng
    % find the sample corresponds to the average number of iterations for each alpha:
    
    NIter(alpha_idx) = mean(conv(alpha_idx, :)); % avg number of iterations.
    
    % compute the index that corresponds to the value which is the closest
    % to the avg number of iterations.
    [~, S(alpha_idx)]= min(abs(conv(alpha_idx, :) - NIter(alpha_idx).'))
    NIter(alpha_idx) = conv(alpha_idx, S(alpha_idx))
  end
  
  maxNIter = max(NIter);
  figure; hold on; grid on;
  	  plot(0:maxNIter/2, WR_vs_iter(S(1 ), 1:maxNIter/2+1, 1), 'm<-', 'linewidth', 2)
      plot(0:maxNIter/2, WR_vs_iter(S(10), 1:maxNIter/2+1, 10), 'ro-', 'linewidth', 2)

  xlabel("number of iterations", 'FontSize', 15);
  ylabel("Weighted sum rate", 'FontSize', 15);
  legend('distance-based alpha', 'uniform alpha', 'FontSize', 15);
end
