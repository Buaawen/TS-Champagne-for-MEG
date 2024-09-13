function [s]=   WEN_champ(B,L,varargin)
%==========================================================================
% Discrition:

% Input: B :             MEG measurements
%        L :             Lead field matrix

%        varargin: 
%                  max_iter :   maximum iterations
%                  epsilon  :   stop criterion


% Output: s ---------- estimated source

% example:  [s]=WEN_champ(B,L,'epsilon',1e-4,'max_iter',50);

%==========================================================================

[nSensor,nSource]=size(L);                      % number of sensors and sources
nSnap=size(B,2);                                % number of time points  
Y=B;


% Default Control Parameters
epsilon                 = 1e-4;             % stop criterion
MAX_ITER                = 50;               % maximum iteration numbers


if nargin>2
   for inar = 1:2:length(varargin)
       Param = lower(varargin{inar});
       Value = varargin{inar+1};
       switch Param
           case 'maxiter'
               MAX_ITER = Value;
           case 'epsilon'
               epsilon = Value;    
       end
   end   
end









l0          = ones(nSensor,1);  % initial values of system noise's convariance

q0          = ones(nSource,1);  % initial values of system source's convariance
cost_old   = 0;
s0         = zeros(nSource,nSnap);                 % initial state estimates

evidence   = zeros(MAX_ITER,1);                    % evidence values

%%===================================================================
%                      iteration
%==========================================================================
    fprintf('\nRunning Champagne for MEG source localization...\n');
for iter = 1:MAX_ITER


    %-----------estimate s0 s1----------------------%

    sigma_y0=diag(l0)+L*(q0.*L');
    invsigma_y0=inv(sigma_y0);
    s0=q0.*L'*invsigma_y0*(B);

    %--------------------------------update q0 q1---------------------------------%
    
    
    for i = 1:nSource
        %----------Gradient descent----------%
        q0(i) = norm(s0(i,:),'fro')/sqrt(nSnap*L(:,i)'*invsigma_y0*L(:,i));
    end
%--------------------------update l0---------------------------------%
     d_y=Y-L*s0;
    for i = 1:nSensor
        %----------Gradient descent----------%
        l0(i) = norm(d_y(i,:),'fro')/sqrt(nSnap*invsigma_y0(i,i));
    end

%--------------------------check stop condition---------------------------%

%% identical {a_i}
cost = -0.5*(trace(B*B'/sigma_y0)+ nSnap*log(det(sigma_y0)) + nSensor*nSnap*log(2*pi));% log p(B,a,q)

MSE = ((cost-cost_old)/cost);
cost_old = cost;
if 1
    disp(['Champagne Iteration: ',num2str(iter),'       MSE: ',num2str(MSE),' ])%,'  ])
end
if abs(MSE)<epsilon
    break;
end
evidence(iter) = cost;







 end
 fprintf('\n Champagne Finished!\n');
 disp(' ');
s=s0;
end