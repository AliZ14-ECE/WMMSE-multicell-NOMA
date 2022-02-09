% Apply the WMMSE algorithim AND SAVE THE RESULTS TO FILE!!, BE CAREFULL!
% The parameters should be assigned earlier.

Pmax = 10^(P/10); % Maximum power for each cell.
Pm = Pmax*ones(1,NC);

h = zeros(NU,NC,NC,num_reals);
d = zeros(NU,NC,NC,num_reals);

for c=1:NC
    num_cell_reals = sum(in(:,c),1);
    cell_indecies = find(in(:,c)==1);
    if(num_cell_reals < num_reals)
        i = 0;
        f = H(:,c,:,cell_indecies);
        dis = D(:,c,:,cell_indecies);
        while (i*num_cell_reals < num_reals)
            start_point = i*num_cell_reals+1;
            end_point = min((i+1)*num_cell_reals, num_reals);
            h(:,c,:,start_point:end_point) = f(:,:,:,1:end_point-start_point+1);
            d(:,c,:,start_point:end_point) = dis(:,:,:,1:end_point-start_point+1);
            i = i + 1;
        end
    else
        h(:,c,:,:) = H(:, c, :, cell_indecies(1:num_reals)); d(:,c,:,:) = D(:, c, :, cell_indecies(1:num_reals));
    end
end
R_vs_iter = zeros(num_reals, numIter, alpha_rng_length);
WR_vs_iter= zeros(num_reals, numIter, alpha_rng_length);

% The (numIter -1) scalling factor is for the case where the convergence
% condition not satisfied, so this channel takes the maximum number of
% iterations (minus 1, to avoid errors after the while loop).
conv= (numIter -1) * ones(alpha_rng_length,num_reals);
 
for alpha_idx = alpha_rng
Powers = zeros(NU,NC,num_reals);
for s = 1:num_reals
        
    % print each 500 realizations
    if mod(s, 500) == 0
        to_disp = sprintf('\t%d. %d', P, s);
        disp(to_disp)
    end

    HH = abs(h(:,:,:,s)); % the abs value of h for each channel realization
    dd = d(:,:,:,s); % the distance for each channel realization
    
    % Next, define weights. Compute alpha depending on alpha_idx.
    alpha = alpha_computation(HH, dd, NU, NC, alpha_idx);

    g = ones(NU, NC);      % Receivers' gains.
    w = ones(NU, NC);      % Weights of the WMMSE problem.
    A = zeros(NU, NC);     % An intermediate variable.
    lambda = zeros(1, NC); % The Lagrange multipliers

    v_init = sqrt(Pmax*rand(NU,NC)); % Initialize v's randomly
    vs = v_init;
    
    % Calculating the corresponding (initial ?) values of g & w for all users in all cells.
    for c=1:NC
        for u=1:NU
            % Compute the intra-cell interference.
            intra = HH(u, c, c)^2*sum(abs(vs(1:u-1, c)).^2);
            
            %Compute the inter-cell interference.
            inter = 0;
            for k=1:NC
                if k~=c
                    inter = inter + HH(u, c, k)^2*sum(abs(vs(:, k)).^2);
                end
            end
            
            g(u, c) = conj(h(u, c, c, s)*vs(u, c))/(abs(h(u, c, c, s) * vs(u, c))^2 + inter + intra + nvar);
            w(u, c) = 1/(abs(g(u, c)* h(u, c, c, s)*vs(u, c)-1)^2 + abs(g(u, c))^2*(inter+intra+nvar));
            A(u, c) = alpha(u, c) * w(u, c) * g(u, c);
        end
    end

    vnew = 0; vold = 0; iter=0;
    while(iter <= numIter)
        iter = iter+1;

        % update what follows for the complex v !!!
        
        % This loop is to compute the v's for each cell
        for c = 1:NC
            
            % Compute the intercell interference for the v-equation
            inter = 0;
            for k=1:NC
                if k~=c
                    inter = inter + sum(A(:, k) .* conj(g(:, k)) .*HH(:, k, c).^2);
                end
            end

            % compute the total power and check if it is less or equal to the maximum power.
            % If the power budget constraint is satisfied, then this value of v is correct and lambda = 0.
            % Else, compute the value of lambda and then compute v according to it.
            
            for u=1:NU
                vs(u,c) = conj( A(u, c) * h(u, c, c, s) ) / (sum(A(u:NU, c) .* conj(g(u:NU, c)).*HH(u:NU,c,c).^2) + inter);
            end
            compPower = sum(abs(vs(:,c)).^2) - Pm(c); % must be <= 0, for lambda = 0.
            
            if compPower > 0
                l = bisection(Pmax, A, h(:, :, :, s), g, inter, c, NU); % Finding lambda using bisection search.
                if numel(l(l>=0)) > 1
                    error("numel(lambda(c))> 1");
                end
                if numel(l(l>=0)) == 0
                    "there is no lambda that satisfies the KKT conditions for this iteration and this cell." %%!
                    error("numel(lambda(c))== 0");
                end
                
                lambda(c) = l(l>=0);
                for u=1:NU
                    vs(u,c) = conj( A(u, c) * h(u, c, c, s) ) / (sum(A(u:NU, c) .* conj(g(u:NU, c)).*HH(u:NU,c,c).^2) + inter + lambda(c));
                end

                compPower = sum(abs(vs(:,c)).^2) - Pm(c) ; % <= 0
                if abs(compPower) > 1e-7
                    error("the solution of the quatric equation is wrong!");
                end
            end
        end
        
        % calculate the rate after computing the powers.
        R_vs_iter(s, iter, alpha_idx) = sum(rate(NC, NU, HH, vs, nvar), 'all'); % the sum rate
        WR_vs_iter(s, iter, alpha_idx)= sum(Wrate(NC, NU, HH, vs, alpha, nvar), 'all'); % the sum weighted rate

        % check if the algorithm comverges.
        vold = vnew;
        vnew = sum(log2(w),'all');
        if abs(vnew-vold) < epsilon
           conv(alpha_idx, s) = iter;
           break;
        end
        
        % Compute g and w.
        for c=1:NC
            for u=1:NU
              intra = HH(u, c, c)^2*sum(abs(vs(1:u-1, c)).^2);
            
            %Compute the inter-cell interference.
            inter = 0;
            for k=1:NC
                if k~=c
                    inter = inter + HH(u, c, k)^2*sum(abs(vs(:, k)).^2);
                end
            end
            
            g(u, c) = conj(h(u, c, c, s)*vs(u, c))/(abs(h(u, c, c, s) * vs(u, c))^2 + inter + intra + nvar);
            w(u, c) = 1/(abs(g(u, c)* h(u, c, c, s)*vs(u, c)-1)^2 + abs(g(u, c))^2*(inter+intra+nvar));
            A(u, c) = alpha(u, c) * w(u, c) * g(u, c);
            end
        end
    end
    
     R_vs_iter(s, conv(alpha_idx, s)+1:end, alpha_idx) =  R_vs_iter(s, conv(alpha_idx, s), alpha_idx);
    WR_vs_iter(s, conv(alpha_idx, s)+1:end, alpha_idx) = WR_vs_iter(s, conv(alpha_idx, s), alpha_idx);

    R = rate(NC, NU, HH, vs, nvar);
    Rmax = rate(NC, NU, HH, v_init, nvar);
    
    WR = Wrate(NC, NU, HH, vs, alpha, nvar);
    WRmax = Wrate(NC, NU, HH, v_init, alpha, nvar);
    
    tdma_rate = tdmaRates(NC,NU,HH,alpha,Pmax,nvar);


    R_sum    = sum(R   , 'all');
    Rmax_sum = sum(Rmax, 'all');
    
    WR_sum    = sum(WR   , 'all');
    WRmax_sum = sum(WRmax, 'all');
    
    tdma_sum = sum(tdma_rate, 'all');
    
    R_sums(s) = R_sum;
    Rmax_sums(s) = Rmax_sum;
    
    WR_sums(s) = WR_sum;
    WRmax_sums(s) = WRmax_sum;
    
    tdma_rates(s) = tdma_sum;
    
    Powers(:,:,s) = vs; 
end

if (executedFrom == '3')
    file_name = sprintf('WMMSE_for_NU/WMMSE_%dx%dpower%dalpha%dabs.mat', NC, NU, P, alpha_idx);
    save(file_name, 'Powers', 'conv', 'R_sums', 'Rmax_sums', 'WR_sums', 'WRmax_sums', 'tdma_rates');
elseif (executedFrom == '2')
    file_name = sprintf('WMMSE_for_powers/WMMSE_%dx%dpower%dalpha%dabs.mat', NC, NU, P, alpha_idx);
    save(file_name, 'Powers', 'conv', 'R_sums', 'Rmax_sums', 'WR_sums', 'WRmax_sums', 'tdma_rates');
end

end

if (executedFrom == '1')
    file_name = sprintf('WMMSE_for_conv/WMMSE_%dx%dpower%dabs.mat', NC, NU,P);
    save(file_name, "conv", "WR_vs_iter");
end

