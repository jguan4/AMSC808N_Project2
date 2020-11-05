function [W,H,errs] = projected_lee_seung(A,k,tol,max_iter)
alpha = 0.001;
[N,d] = size(A);
perturb = 1;
rng(1);
% W = rand(N,k);
% H = rand(k,d);
W = max(A(:,1:k),0);
H = eye(k,d);
H = H + 0.1*rand(k,d);
R = A-W*H;
err0 = norm(R,'fro')^2;
errs = [err0];

W_new = W+alpha*R*H';
W_new = max(W_new,0);
H_new = H+alpha*W'*R;
H_new = max(H_new,0);
W = W_new;
H = H_new;
R = A - W*H;
err = norm(R,'fro')^2;
errs = [errs,err];
iter =1;
while err>tol*err0 && iter<max_iter
    W_new = W.*(A*H')./(W*H*H'+perturb*rand(N,k));
    H_new = H.*(W'*A)./(W'*W*H+perturb*rand(k,d));
    W = W_new;
    H = H_new;
    R = A - W*H;
    err = norm(R,'fro')^2;
    errs = [errs,err];
    iter = iter + 1;
end
end
