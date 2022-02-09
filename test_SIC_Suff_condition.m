function p = test(NC, NU, inner_radius, minR_ratio,h, Pm, nv, seed)
h = abs(h).^2;
h = permute(h, [3,2,1]);
h_cell = zeros(NC,NU);
for c=1:NC
    h_cell(c,:) = h(c,c,:);
end
Pmax = Pm*ones(1,NC);
U = NU*ones(1,NC);
nvar = nv*ones(NC,NU);
p = 10*ones(1,NC);
n = 10*ones(1,NC);
for c=1:NC
[p,n] = func_suff_condition_optSIC(NC,U,NU,h,h_cell,Pmax,nvar);
end