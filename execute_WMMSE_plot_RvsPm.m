alpha_rng_length = 10; % number of alpha values to compute.
NU = 10; % #USers in each cell.
NC = 4;
nvar = 1.9905e-08; % Noise Variance
epsilon = 1e-5; % For convergence test.
inner_radius = 500; 
minR_ratio = 0.01;
numIter = 2000;
num_reals = 1000;
alpha_rng = [1, alpha_rng_length];
P_max_idx = 20;
seed = 1;

% generate and save the data if not exist
for P=0:P_max_idx
    fileName1  = sprintf('WMMSE_for_powers/WMMSE_%dx%dpower%dalpha%dabs.mat', NC, NU,P, 1 );
    fileName10 = sprintf('WMMSE_for_powers/WMMSE_%dx%dpower%dalpha%dabs.mat', NC, NU,P, 10);
    if (~exist(fileName1, 'file') || ~exist(fileName10, 'file'))
    clear H in D;
    fileName = sprintf('channels_for_powers/Channels%dx%dpower%d.mat', NC, NU, P);
    load(fileName,'H', 'in', 'D'); 
    executedFrom = '2';
    execute_WMMSE
    end
end

clear RR RR_max tdma WRR WRR_max convv;
RR = zeros(10, P_max_idx+1);
RR_max = zeros(10, P_max_idx+1);
WRR = zeros(10, P_max_idx+1);
WRR_max = zeros(10, P_max_idx+1);
tdma = zeros(10, P_max_idx+1);
fdma = zeros(10, P_max_idx+1);
convv = zeros(10, P_max_idx+1); 
for alpha_idx = alpha_rng
    for P = 1:P_max_idx+1
        file_name = sprintf('WMMSE_for_powers/WMMSE_%dx%dpower%dalpha%dabs.mat', NC, NU,P-1, alpha_idx);
        load(file_name, 'Powers', 'conv', 'R_sums', 'Rmax_sums', 'WR_sums', 'WRmax_sums', 'tdma_rates');

        RR(alpha_idx,P) = mean(R_sums);
        RR_max(alpha_idx, P) = mean(Rmax_sums);
        WRR(alpha_idx, P) = mean(WR_sums);
        WRR_max(alpha_idx, P) = mean(WRmax_sums);
        tdma(alpha_idx, P) = mean(tdma_rates);
        if (mod(alpha_idx,2) == 1)
            fdma(alpha_idx, P) = mean(tdma_rates); % TDMA, not FDMA. (although the same result would be optained.)
        end
    end
end

figure; hold on; grid on;
plot(0:20, RR(10, :),'m<-', 'linewidth',2);
plot(0:20, RR(1, :), 'bo-', 'linewidth',2);
% plot(0:20, RR(2, :), 'bo--', 'linewidth',2);
% plot(0:20, RR(3, :), 'r*-', 'linewidth',2);
% plot(0:20, RR(4, :), 'r*--', 'linewidth',2);
% plot(0:20, RR(5, :), 'g+-', 'linewidth',2);
% plot(0:20, RR(6, :), 'g+--', 'linewidth',2);
% plot(0:20, RR(7, :), 'cv-', 'linewidth',2);
% plot(0:20, RR(8, :), 'cv--', 'linewidth',2);
% plot(0:20, RR(9, :), 'm<-', 'linewidth',2);


plot(0:20, RR_max(1, :), 'ks-.', 'linewidth',2);
plot(0:20, fdma(1, :), 'gd--', 'linewidth',2);

xlabel('Maximum power for each BS (dBW)', 'FontSize', 15);
ylabel('Sum rate (bits/s/Hz)', 'FontSize', 15);

legend('uniformly-distributed-alpha WMMSE', 'distance-based-alpha WMMSE', 'uniform power allocation', 'OMA', 'FontSize', 15);