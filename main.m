clear all
close all
load('data.mat');
load('J.mat');
load('leadfield.mat');
load('sourcemodel.mat');
load('VertConn.mat');
%% Whiten measurements and lead field matrix
StimTime=100;%time_noise
Bpre = B(:,1:StimTime);
Cov_n = (Bpre - repmat(mean(Bpre,2),1,StimTime)) * (Bpre - repmat(mean(Bpre,2),1,StimTime))'/(StimTime - 1);
rnkC_noise = rank(single(Cov_n));
variance = diag(Cov_n);
isPca = 1;
if isequal(Cov_n, diag(variance))
    isPca = 0;
end
[VV,D] = eig(Cov_n);
D = diag(D);
[D,I] = sort(D,'descend');
VV = VV(:,I);
if ~isPca
    D = 1./D;
    W = diag(sqrt(D)) * VV';
else
    D = 1 ./ D;
    D(rnkC_noise+1:end) = 0;
    W = diag(sqrt(D)) * VV';
    W = W(1:rnkC_noise,:);
end
% clear VV D I
L = W*Gain;
B = W*B;
%
[s_TS_Champ] = TS_Champ(L,B,VertConn,'tau',0.01);
%% plot
alpha=0.1;
xx=abs(s_TS_Champ(:,155));
xx=xx.*(xx>alpha*max(xx));
xxx=xx/max(xx);
figure;
ft_plot_mesh(sourcemodel_fig,'vertexcolor',xxx,'colormap','jet');%'*RdBu'
%
xx=abs(J(:,155));
xxx=xx/max(xx);
figure;
ft_plot_mesh(sourcemodel_fig,'vertexcolor',xxx,'colormap','jet');%'*RdBu'

