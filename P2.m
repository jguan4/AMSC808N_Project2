close all;
clear all;

A = readtable('MovieRankings36.csv');
A = table2array(A);
[N,d] = size(A);
inds = find(isnan(A));
C = A;
C(inds) = 0;
inds_n = find(~isnan(A));
Omega = zeros(N,d);
Omega(inds_n) = 1;

k = 5;
lambda = 0.0001;
tol = 1e-4;
max_iter = 50;

[X,Y,errs] = low_rank_factorization(C,lambda,Omega,k,max_iter);
figure;
plot(errs);
set(gca,'Fontsize',14);
xlabel('k','Fontsize',14);
ylabel('norm(M-CUR)_F','Fontsize',14);

[X,Y,errs] = nuclear_norm(C,lambda,Omega,k,max_iter);
figure;
plot(errs);
set(gca,'Fontsize',14);
xlabel('k','Fontsize',14);
ylabel('norm(M-CUR)_F','Fontsize',14);