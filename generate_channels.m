% This file is to generate and save channels that satisfy the SIC sufficient condition 

% % BE CAREFUL OF THE 'SAVE' STATEMENT!!!

NS = 10^4;
inner_radius = 500; minR_ratio = 0.01; seed = 1;
Pmax = 16; % dBw
nvar_dBm = -147;
BW = 10 * 10^6;
nvar = 10^(nvar_dBm/10) * 10^-3 * BW; % linear scale (Watts)
rng(seed);
Pm = 10^(Pmax/10);
tic;

for Pmax = 0:20
for NC = 4:6
for NU = 2:20
  %execute only if:
  if( (Pmax == 16) || (NC == 4 && NU == 10) )
    disp(Pmax);
    A = zeros(NC,NU,2,NS);
    R = inner_radius - minR_ratio*inner_radius;      % effective cell radius
    for c = 1:NC
        for u = 1:NU
            d = sum(rand(2,NS),1) .* R;              % user distributed in the cell uniformly
            d(d>R) = 2*R - d(d>R);
            A(c,u,1,:) = d + minR_ratio*inner_radius;               % real MS location
            A(c,u,2,:) = 2*pi*rand(1,NS);
        end               
    end
    
    A = permute(A, [4,2,3,1]);
    H = zeros(NU,NC,NC,NS);
    D = zeros(NU,NC,NC,NS);
    in = ones(NS, NC);

    for s=1:NS
        [h, distances, ms, Cell] = generate_IBC_channel(NU, inner_radius, NC, minR_ratio, seed, squeeze(A(s,:,1,:)), squeeze(A(s,:,2,:)), 0);
        H(:,:,:,s) = h;
        D(:,:,:,s) = distances;
        p = test_SIC_Suff_condition(NC, NU, inner_radius, minR_ratio,h, Pm, nvar, seed);
        for c = 1:NC 
            if (p(c)~=0)
                in(s,c) = 0;
                H(:, c, :, s) = NaN;
            end
        end
    end

    if (Pmax == 16)
        file_name = sprintf('channels_for_NU/Channels%dx%dpower%d.mat', NC, NU, Pmax);
        save(file_name, 'H', 'in','D', 'nvar');
    end
    if (NC == 4 && NU == 10)
        file_name = sprintf('channels_for_powers/Channels%dx%dpower%d.mat', NC, NU, Pmax);
        save(file_name, 'H', 'in','D', 'nvar');
    end
        
    toc;
  end
end
end
end
