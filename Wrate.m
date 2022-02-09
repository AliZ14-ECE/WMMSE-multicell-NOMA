% Please update this for the complex v !!!

function R = Wrate(NC, NU, H, v, alpha, nvar)
R = zeros(NU, NC);
intra = 0;
for c = 1:NC
  for i = 1:NU
      inter = 0;
      for k=1:NC
          if k~=c
              temp = H(i,c,k)^2*sum(abs(v(:, k)).^2);
              inter = inter + temp;
          end
      end
      if i == 1
          intra = 0;
      else
          intra = sum((H(1:i-1 , c, c).*abs(v(1:i-1, c))).^2);
      end
      R(i, c) = alpha(i, c) * log2(1+(H(i, c, c)*abs(v(i, c)))^2/(inter+intra+nvar));
  end
end

end