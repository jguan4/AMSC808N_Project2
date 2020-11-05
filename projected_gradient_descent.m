function [W,H,errs] = projected_gradient_descent(A,k,tol,max_iter)
alpha = 0.001;
[N,d] = size(A);
rng(1);
% W = rand(N,k);
% H = rand(k,d);
W = max(A(:,1:k),0);
H = eye(k,d);
H = H + 0.1*rand(k,d);
R = A-W*H;
err0 = norm(R,'fro')^2;
errs = [err0];
err = err0;
iter = 0;
while err>tol*err0 && iter<max_iter
    W_new = W+alpha*R*H';
    W_new = max(W_new,0);
    H_new = H+alpha*W'*R;
    H_new = max(H_new,0);
    W = W_new;
    H = H_new;
    R = A - W*H;
    err = norm(R,'fro')^2;
    errs = [errs,err];
    iter = iter+1;
end
end