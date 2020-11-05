close all;
clear all;

A = readtable('MovieRankings36.csv');
A = table2array(A);
[N,d] = size(A);
inds = find(isnan(A));
inds_c = setdiff(1:N*d,inds);
[row,col] = ind2sub([N,d],inds);
row_c = setdiff(1:N, row);
B = A(row_c,:);
C = A;
C(inds) = 0;

k = 5;
tol = 1e-4;
max_iter = 1000;

[W,H,errs] = projected_gradient_descent(C,k,tol,max_iter);
figure;
plot(errs);

[W,H,errs] = lee_seung(C,k,tol,max_iter);
figure;
plot(errs);

[W,H,errs] = projected_lee_seung(C,k,tol,max_iter);
figure;
plot(errs);