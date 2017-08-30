function y=cell2_fn(x,k1,k2,temp)
   y=zeros(k2-k1+1,1);
   for j=k1:k2
      tmp=x(temp==j);
      y(j-k1+1,1)=tmp(1);

   end

return
end