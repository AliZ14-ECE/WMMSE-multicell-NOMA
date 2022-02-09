function alpha = alpha_computation(HH, dd, NU, NC, alpha_idx)
%%% computes the wieghts (alpha) in different ways acording to the [alpha_idx]

% distance-based alpha
if alpha_idx ==1
    for c=1:NC
        alpha_sum = 0;
        for u=1:NU
            alpha(u,c) = (dd(u,c,c)/((sum(dd(u,c,:))-dd(u,c,c))))^2;
            alpha_sum = alpha_sum+alpha(u,c);
        end
        alpha(:,c) = alpha(:,c)/alpha_sum;
    end
    
%----------
elseif alpha_idx ==2
    for c=1:NC
        alpha_sum = 0;
        for u=1:NU
            alpha(u,c) = (dd(u,c,c)/((sum(dd(u,c,:))-dd(u,c,c))))^2;
        end
        alpha(:,c) = alpha(:,c);
    end
    
elseif alpha_idx ==3
    for c=1:NC
        alpha_sum = 0;
        for u=1:NU
            alpha(u,c) = (HH(u,c,c)/(sum(HH(u,c,:))-HH(u,c,c)))^-2;
            alpha_sum = alpha_sum+alpha(u,c);
        end
        alpha(:,c) = alpha(:,c)/alpha_sum;
    end
    
elseif alpha_idx ==4
    for c=1:NC
        alpha_sum = 0;
        for u=1:NU
            alpha(u,c) = (HH(u,c,c)/(sum(HH(u,c,:))-HH(u,c,c)))^-2;
        end
        alpha(:,c) = alpha(:,c);
    end
    
elseif alpha_idx ==5
    for c=1:NC
        alpha_sum = 0;
        for u=1:NU
            alpha(u,c) = (HH(u,c,c)/(sum(HH(u,c,:))-HH(u,c,c)))^-1;
            alpha_sum = alpha_sum+alpha(u,c);
        end
        alpha(:,c) = alpha(:,c)/alpha_sum;
    end
        
elseif alpha_idx ==6
    for c=1:NC
        alpha_sum = 0;
        for u=1:NU
            alpha(u,c) = (HH(u,c,c)/(sum(HH(u,c,:))-HH(u,c,c)))^-1;
        end
        alpha(:,c) = alpha(:,c);
    end
    
    
elseif alpha_idx ==7
    for c=1:NC
        alpha_sum = 0;
        for u=1:NU
            alpha(u,c) = (dd(u,c,c)/((sum(dd(u,c,:))-dd(u,c,c))))^1;
            alpha_sum = alpha_sum+alpha(u,c);
        end
        alpha(:,c) = alpha(:,c)/alpha_sum;
    end
    
elseif alpha_idx == 8
    for c=1:NC
        alpha_sum = 0;
        for u=1:NU
            alpha(u,c) = (dd(u,c,c)/((sum(dd(u,c,:))-dd(u,c,c))))^1;
        end
        alpha(:,c) = alpha(:,c);
    end
    
elseif alpha_idx == 9
    alpha = ones(NU, NC);
    alpha = alpha./sum(alpha, 1); %?????
    
% uniform alpha
elseif alpha_idx == 10
    alpha = 0.05*ones(NU, NC);
%----
elseif alpha_idx ==11
    for c=1:NC
        alpha_sum = 0;
        for u=1:NU
            alpha(u,c) = (HH(u,c,c)/(sum(HH(u,c,:))-HH(u,c,c)))^-1.5;
            alpha_sum = alpha_sum+alpha(u,c);
        end
        alpha(:,c) = alpha(:,c)/alpha_sum;
    end
    
elseif alpha_idx ==12
    for c=1:NC
        alpha_sum = 0;
        for u=1:NU
            alpha(u,c) = (HH(u,c,c)/(sum(HH(u,c,:))-HH(u,c,c)))^-1.5;
        end
        alpha(:,c) = alpha(:,c);
    end
    
end

    
