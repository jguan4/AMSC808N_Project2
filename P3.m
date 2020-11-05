close all;
clear all;

[M,y] = readdata();
realizationN = 100;
[Um,Sm,Vm] = svd(M,'econ');

ks = 2:10;
as = 1:8;
res_errs = [];
relative_errs = [];
errs = [];
for k = ks
    Mk = Um(:,1:k)*Sm(1:k,1:k)*Vm(:,1:k)';
    res_err = norm(M-Mk,'fro');
    res_errs = [res_errs,res_err];
    for a = as
        c = a*k;
        r = c;
        relative_err_vec = zeros(1,realizationN);
        errs_vec = zeros(1,realizationN);
        tic
        for i = 1:realizationN
            [C,U,R] = CUR(M,c,r,k);
            diff = M-C*U*R;
            errs_vec(i) = norm(M-C*U*R,'fro');
            relative_err = errs_vec(i)/res_err;
            relative_err_vec(i) = relative_err;
        end
        toc
        relative_err_ave = mean(relative_err_vec);
        errs_ave = mean(errs_vec);
        relative_errs = [relative_errs;[k,a,relative_err_ave]];
        errs = [errs;[k,a,errs_ave]];
    end
end

figure
[n,~] = size(relative_errs);
for i = 1:length(as)
    a = as(i);
    inds = i:as(end):n;
    relative_err_i = relative_errs(inds,3);
    plot(ks,relative_err_i,'o-','LineWidth',2,'DisplayName',['a=',num2str(a)]);
    hold on;
end
legend;
set(gca,'Fontsize',14);
xlabel('k','Fontsize',14);
ylabel('norm(M-CUR)_F/norm(M-M_k)_F','Fontsize',14);
figure
[n,~] = size(errs);
for i = 1:length(as)
    a = as(i);
    inds = i:as(end):n;
    err_i = errs(inds,3);
    plot(ks,err_i,'o-','LineWidth',2,'DisplayName',['a=',num2str(a)]);
    hold on;
end
legend;
set(gca,'Fontsize',14);
xlabel('k','Fontsize',14);
ylabel('norm(M-CUR)_F','Fontsize',14);
figure
plot(ks,res_errs,'o-','LineWidth',2);
set(gca,'Fontsize',14);
xlabel('k','Fontsize',14);
ylabel('norm(M-M_k)_F','Fontsize',14);
