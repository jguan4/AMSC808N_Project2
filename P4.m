close all;
clear all;

[M,y] = readdata();
M_logic = M;
inds = find(M_logic>0);
M_logic(inds) = 1;

iminus = find(y == -1);
iplus = setdiff((1:size(M,1))',iminus);
M_f = M_logic(iminus,:);
M_l = M_logic(iplus,:);

fid=fopen('features_idx.txt','rt');
words = [];
while ~feof(fid)
    l=fgetl(fid);
    if isempty(l)
        continue;
    elseif strfind(l,'#')
        continue;
    else
        word=split(l);
        word = string(word{2});
        words = [words;word];
    end
end
fid=fclose(fid);

%% before preprocessing, print out high leverage words
k = 5;
a = 8;
c = a*k;
r = c;
toggle = 'k_norm';
inds = leverage_score(M,k,k,toggle);

lev_words = words(inds);
fprintf("The %d maximal leverage words are:\n",k)
for i = 1:length(lev_words)
    fprintf("%s\n",lev_words(i));
end
test_proj(M(:,inds), y)

%% preprocessing and plot projection
% take out words appearing often in both classes
pp_f = sum(M(iminus,:),1)/71;
pp_l = sum(M(iplus,:),1)/68;
pps = pp_f.*pp_l;
t = find(pps>0.4);
M1 = M;
M1(:,t) = [];
words1 = words;
words1(t) = [];
M_f(:,t) = [];
M_l(:,t) = [];
M_logic(:,t) = [];

pp_f = sum(M_f,1)/71;
pp_l = sum(M_l,1)/68;
pps = pp_f.*pp_l;
t = find(pps>0.03);
M1(:,t) = [];
words1(t) = [];
M_f(:,t) = [];
M_l(:,t) = [];
M_logic(:,t) = [];

col_sum_m1 = sum(M1,1)./sum(M_logic,1).^2;
col_sum_m1_std = (col_sum_m1-mean(col_sum_m1))/std(col_sum_m1);
temp_inds = find(col_sum_m1_std>1);
words2 = words1;
words2(temp_inds) = [];
M1(:,temp_inds) = [];
M_f(:,temp_inds) = [];
M_l(:,temp_inds) = [];
M_logic(:,temp_inds) = [];

col_max = max(M1)./sum(M_logic,1);
col_max_std = (col_max-mean(col_max))/std(col_max);
temp_inds = find(col_max_std>7);
words2(temp_inds) = [];
M1(:,temp_inds) = [];
M_f(:,temp_inds) = [];
M_l(:,temp_inds) = [];
M_logic(:,temp_inds) = [];

[U,S,V] = svd(M_f,'econ');
s_f = sum(V(:,1:k).^2,2)/k;
[U,S,V] = svd(M_l,'econ');
s_l = sum(V(:,1:k).^2,2)/k;

prob_f = zeros(1,size(M1,2));
prob_i = zeros(1,size(M1,2));
for i = 1:size(M1,2)
    Mi = M1(:,i);
    inds = find(Mi>0);
    finds = find(ismember(inds,iminus));
    prob_f(i) = length(finds)/length(inds);
    prob_i(i) = length(inds)/length(Mi);
end
prob_l = 1-prob_f;
p_f = length(iminus)/length(M1);
p_l = length(iplus)/length(M1);
temp_f = (prob_f-0.5).*s_f';
temp_l = (prob_l-0.5).*s_l';
[~,inds_f] = maxk(temp_f,5000);
[~,inds_l] = maxk(temp_l,5000);

inds = [inds_f,inds_l];
inds = sort(inds);
words3 = words2(inds);
M_test = M1(:,inds);

test_proj(M_test, y)


%% functions
function test_proj(M, y)
svd_projection(M,y);
leverage_projection(M,y);
end

function svd_projection(M,y)
[U,S,V] = svd(M,'econ');
p = V(:,1:2);
proj = M*p;
plot_and_color(proj,y,2)
end

function leverage_projection(M,y)
inds = leverage_score(M,5,5,'k_norm');
A = M(:,inds);
[U,S,V] = svd(A,'econ');
p = V(:,1:2);
proj = A*p;
plot_and_color(proj,y,2)
end

function plot_and_color(P,y,toggle)
[N,d] = size(P);
iminus = find(y == -1);
iplus = setdiff((1:N)',iminus);
switch toggle
    case 2
        figure;
        plot(P(iminus,1),P(iminus,2),'Linewidth',2,'Linestyle','none','Marker','s','color','r','DisplayName','Indiana');
        hold on;
        plot(P(iplus,1),P(iplus,2),'Linewidth',2,'Linestyle','none','Marker','<','color','b','DisplayName','Florida');
        legend;
        set(gca,'Fontsize',14);
        xlabel('p_1','Fontsize',14);
        ylabel('p_2','Fontsize',14);
    case 3
        figure;
        plot3(P(iminus,1),P(iminus,2),P(iminus,3),'Linewidth',2,'Linestyle','none','Marker','s','color','r','DisplayName','Indiana');
        hold on;
        plot3(P(iplus,1),P(iplus,2),P(iplus,3),'Linewidth',2,'Linestyle','none','Marker','<','color','b','DisplayName','Florida');
        legend;
        set(gca,'Fontsize',14);
        xlabel('p_1','Fontsize',14);
        ylabel('p_2','Fontsize',14);
end
end

function inds = leverage_score(M,k,m,toggle)
switch toggle
    case 'k_norm'
        [U,S,V] = svd(M,'econ');
        s = sum(V(:,1:k).^2,2)/k;
    case 'lev_score'
        MTM = M'*M;
        temp = MTM\M';
        s = diag(M*temp);
end
[B,inds] = maxk(s,m);
end

