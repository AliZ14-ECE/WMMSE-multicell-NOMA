function [perc_pair_unsatisfy,num_pairs_unsatisfy]=func_suff_condition_optSIC(B,U,M,h,h_cell,P_max,N0)
SIC_satisfy=zeros(B,M,M); %%SIC_Theo(b,i,k)  k \to i in cell b
Q=zeros(B,B,M,M);
for b=1:B
    for i=1:U(b)
        for k=1:U(b)
            for j=1:B
                if j~=b && (h_cell(b,k)/h_cell(b,i)) < (h(j,b,k)/h(j,b,i))
                    Q(j,b,i,k)=1;
                end
            end
            if h_cell(b,k)/N0(b,k) > h_cell(b,i)/N0(b,i) && h_cell(b,k)/N0(b,k) - h_cell(b,i)/N0(b,i) > ...
                    sum(Q(:,b,i,k).*P_max(1,:)' .* (h(:,b,k) .*h_cell(b,i) - h(:,b,i) .*h_cell(b,k))./(N0(b,k)*N0(b,i)))
                SIC_satisfy(b,i,k)=1;
                SIC_satisfy(b,k,i)=-1;
            end
        end
    end
end
%Total number of user pairs in each cell: M!/(M-2)!2!
for b=1:B
    Tot_pairs(b)=U(b)*(U(b)-1)/2;
end
%number of nonzero elements
for b=1:B
    num_pairs_unsatisfy(b)=Tot_pairs(b)-nnz(SIC_satisfy(b,:,:))./2;
end
perc_pair_unsatisfy=100*(num_pairs_unsatisfy./Tot_pairs); %percentage of total user pairs