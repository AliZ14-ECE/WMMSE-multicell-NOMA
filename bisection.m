function lambda = bisection(Pm, A, h, g, inter, c, NU)
  %%% why are you passing NU?? it is the first dimension of h, g, and A.
  %%% why are you passing c?? why not just pass Pm(c), A(:, c), ... etc.?
    lambda = 0;
    start_point = 0;
    end_point = 100;
    while true
        x = sum_power(end_point,A, h, g, inter, c, NU);
        if x > Pm
          % start_point = end_point;       %%% isn't this better?? so the search converges faster???
            end_point = end_point*2;
        else
            
            break;
        end
    end

    while true
        lambda = (end_point+start_point)/2;
        s = sum_power(lambda,A, h, g, inter, c, NU);
        if abs(s-Pm) < 1e-7
            break;
        end
        if s>Pm
            start_point = lambda;
        else
            end_point = lambda;
        end
    end
end

function s = sum_power(lambda,A, h, g, inter, c, NU)
    vs = zeros(size(h,1),1);
    for u=1:NU
        vs(u) = conj( A(u, c) * h(u, c, c) ) / (sum(A(u:NU, c) .* conj(g(u:NU, c)).*abs(h(u:NU,c,c)).^2) + inter + lambda);
    end
    s = sum(abs(vs).^2);
end
