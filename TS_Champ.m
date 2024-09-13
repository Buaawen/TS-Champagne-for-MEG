function [S] = TS_Champ(L,B,VertConn,varargin)
% This script is to solve MEG source imaging problem under TS_Champagne
% =================================================================
% [S] = TS_Champ(L,B,VertConn,varargin) 
% Input:
%   -L: lead-field matrix
%   -B: MEG measurements
%   -VertConn:  Connectivity between dipoles
% Output:
%   -S: reconstructed source
%==================================================================
% Author: Wen Li
% Data: 2024/9/13


%% Initialize 
[~,nSource] = size(L); 

k=5;
p=6;
delta=0.01;
W=eye(nSource);


if nargin>3
   for inar = 1:2:length(varargin)
       Param = lower(varargin{inar});
       Value = varargin{inar+1};
       switch Param
           case 'delta'
               delta = Value;
           case 'k'
               k = Value;    
           case 'p'
               p = Value;
       end
   end
end



[S0] = WEN_champ(B,L,'maxiter',50);%First time using the Champagne method


QS0=diag(S0*S0');
indexs=find(QS0>delta*max(QS0));

Q=VertConn+speye(nSource);
QZ={};
QQ=speye(nSource);
for i=1:15%11
    QZ{i}=double(QQ>0); 
    QQ=QQ*Q;
end
Q5=QZ{k+1};%Neighborhoods of the first k orders
indz=[];
for i=1:length(indexs)
    ind=find(Q5(indexs(i),:)>0);
    indz=[indz,ind];
end
indz=unique(indz);

W2=sparse(nSource,(p+1)*length(indz));%calculate W2
num=0;
for i=1:p+1
    for j=1:length(indz)
        num=num+1;
        W2(:,num)=QZ{i}(:,indz(j));
    end
end


%%


[S01] = WEN_champ(B,L*W2);%Second time using the Champagne method
S=W2*S01;

end


