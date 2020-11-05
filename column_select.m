function C = column_select(A,c,k)
[U,S,V] = svd(A,'econ');
ps = sum(V(:,1:k).^2,2)/k;
cps = min(1,c*ps);
probs = rand(size(cps));
inds = find(probs<cps);
C = A(:,inds);
end
