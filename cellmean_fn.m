function y=cellmean_fn(x,k1,k2,temp)

for i=k1:k2
   y(i-k1+1)=mean(x(temp==i));
end

return
end