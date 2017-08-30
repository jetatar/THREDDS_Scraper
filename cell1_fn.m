function y=cell1_fn(x,k1,k2,temp)
   tmp=x(temp==k1);
   y=ones(k2-k1+1,1)*tmp(1);

return
end