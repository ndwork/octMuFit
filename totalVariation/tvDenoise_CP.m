function [xn,costs] = tvDenoise_CP(b,params)
% This is based on deblur_mixedTV_CP_orVersion.m.
% orVersion stands for "overrelaxed version" of Chambolle-Pock.

% We solve minimize ||x - b||_1 + gamma*|| Gx ||
% subject to 0 <= x_{ij}  for all i,j.
% A convolves with mask using replicate boundary conditions.
% G is a discrete gradient operator using symmetric boundary conditions.
% A = P + S , G = D + R.

maxIter = params.maxIter;
showTrigger = params.showTrigger;
gamma = params.gamma;
theta = params.theta;
overRelax = params.overRelax;

evalIsoNorm = @(y) sum(sum(sqrt(y(:,:,1).^2 + y(:,:,2).^2)));

[numRows,numCols] = size(b);
numPix = numRows*numCols;
%%%%%%%%%%%% Set up operators.

applyA = @(x) x;
applyATrans = @(x) x;

applyG = @(x) computeGrad2D_neumannBCs(x);
applyGTrans = @(x) -computeDiv_neumannBCs(x);

applyK = @(x) cat(3,applyA(x),applyG(x));
applyKTrans = @(z) applyATrans(z(:,:,1)) + applyGTrans(z(:,:,2:3));

%%%%%%%%%%%%%%%%%%%

nrmA = 1;
nrmK = sqrt(8 + nrmA^2);

% nrmK_est = estimateNormByPowerIteration(applyK,applyKTrans,b,500);

% tau = params.tau;
% sigma = 1/(tau*nrmK^2);

tau = 1/nrmK;
sigma = 1/nrmK;

xn = b;
yn = applyK(xn);
xnBar = xn;

costs = [];
figure('Name','xn inside Chambolle-Pock')

for n = 1:maxIter 
    
    ATransyn_1 = applyATrans(yn(:,:,1));
    GTransyn_2 = applyGTrans(yn(:,:,2:end));
    KTransyn = ATransyn_1 + GTransyn_2;          
    
    % compute proxTauG, which simply projects onto C.
%     xnp1 = min(xn - tau*KTransyn,1);
    xnp1 = xn - tau*KTransyn;
    xnp1 = max(xnp1,0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
            
    z = yn + sigma*applyK(2*xnp1 - xn);
    
    % compute proxSigmaFStar evaluated at z.
    z1 = z(:,:,1);
    z2 = z(:,:,2:3);
    
    term = z1 - sigma*b;    
    ynp1_1 = min(term,1);
    ynp1_1 = max(ynp1_1,-1);
        
    ynp1_2 = gamma*projectOntoDualIsoNormUnitBall(z2/gamma);
    ynp1 = cat(3,ynp1_1,ynp1_2); % This is dual feasible.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    xn = overRelax*xnp1 + (1-overRelax)*xn;
    yn = overRelax*ynp1 + (1-overRelax)*yn; 
        
    Axnp1 = applyA(xnp1);
    Gxnp1 = applyG(xnp1);    
            
    cost = sum(sum(abs(Axnp1-b))) + gamma*evalIsoNorm(Gxnp1);
    costs = [costs,cost];    
    
    if mod(n,showTrigger) == 0
        imshow(xn,[])
        disp(['Chambolle-Pock iteration is: ',num2str(n)])
        disp(['primal cost is: ',num2str(cost)])         
        keyboard        
    end     
    
end

end