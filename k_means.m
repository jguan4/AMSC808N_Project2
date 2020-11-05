function means = k_means(A,k,tol)
[N,d] = size(A);
k_inds = randsample(N,k);
means = A(k_inds,:);
err = norm(means);
while err>tol
    old_means = means;
    dists = pdist2(means,A,'euclidean');
    [min_dists,labels] = min(dists,[],1);
    for i = 1:k
        cluster_i = find(labels==i);
        means(i,:) = mean(A(cluster_i,:),1);
    end
    err = norm(old_means-means);
end
end
