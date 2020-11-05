function [X,Y,errs] = low_rank_factorization(A,lambda,Omega,k,max_iter)
[N,d] = size(A);
Omega_t = Omega';
A_omega = Omega.*A;
A_t_omega = Omega_t.*A';
X = A(:,1:k);
Y = eye(d,k);
R = Omega.*(A-X*Y');
err0 = norm(R,'fro')^2;
errs = [err0];
iter = 0;
while iter<max_iter
    for i = 1:N
        Omega_i = Omega(i,:);
        inds = find(Omega_i>0);
        Y_omegai = Y(inds,:);
        a_omegai = A_omega(i,inds)';
        xi = (Y_omegai'*Y_omegai+lambda*eye(k))\(Y_omegai'*a_omegai);
        X(i,:) = xi';
    end
    
    for i = 1:d
        Omega_i = Omega(i,:);
        inds = find(Omega_i>0);
        X_omegai = X(inds,:);
        a_omegai = A_omega(i,inds)';
        yi = (X_omegai'*X_omegai+lambda*eye(k))\(X_omegai'*a_omegai);
        Y(i,:) = yi';
    end
    R = Omega.*(A-X*Y');
    errs = [errs, norm(R,'fro')^2];
    iter = iter + 1;
end
end