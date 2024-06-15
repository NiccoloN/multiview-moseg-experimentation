%% function to do co-regularization
%
%   Input: K - cell arrays of kernels
%          nMotion - number of motions
%          lambda - co-regularization parameter
%
function [U_CoReg, itr ,loss, exitinfo,  RedLapTrace] = func_CoRegularize_eig(K,nMotion,lambda,epsilon,MaxItr)

%% Convergence Parameters
if ~exist('MaxItr','var')
    MaxItr = 30;
end

if ~exist('epsilon','var')
    epsilon = 1e-8;
end

nKernel = size(K,3);

%% Intialize Each Spectral Embedding
U = [];

U_tmp_temp = zeros(size(K,1), size(K,2));

for k_i = 1:nKernel
    
    D = diag(sum(K(:,:,k_i),1));
    L(:,:,k_i) = eye(size(D)) - D^-0.5 * K(:,:,k_i) *D^-0.5;
    
    [U_tmp,S,V] = svd(L(:,:,k_i));
    
    U(:,:,k_i) = U_tmp(:,end-nMotion+1:end);

    U_tmp_temp = U_tmp_temp + U_tmp;
    
end

%% Do Co-Regularization
exitinfo.reason = '';

U_CoReg = U;

UUt = zeros(size(K));
for k_i = 1:nKernel
    UUt(:,:,k_i) = U(:,:,k_i)*U(:,:,k_i)';
end

itr = 1;
RedLapTrace = 0;
L_CoReg = zeros(size(K));
L_CoReg_Red = zeros(nMotion,nMotion,nKernel);



while true
    
    for k_i = 1:nKernel
        
        %%% Compute CoReg Kernels
        
        UUt_k = sum(UUt,3) - U_CoReg(:,:,k_i)*U_CoReg(:,:,k_i)';
%         L_CoReg(:,:,k_i) = L(:,:,k_i) - lambda/nMotion * UUt_k;
        L_CoReg(:,:,k_i) = L(:,:,k_i) - lambda * UUt_k;

%         [U_tmp,S,V] = svd(L_CoReg(:,:,k_i));
        
        [Vec,D] = eig((L_CoReg(:,:,k_i)'+L_CoReg(:,:,k_i))/2);
        
        [~,idx] = sort(diag(D),'ascend');
        
        U_CoReg(:,:,k_i) = Vec(:,idx(1:nMotion));

        L_CoReg_Red(:,:,k_i) = L_CoReg(idx(1:nMotion), idx(1:nMotion), k_i);
        
%         U_CoReg(:,:,k_i) = U_tmp(:,end-nMotion+1:end);
        
        UUt(:,:,k_i) = U_CoReg(:,:,k_i)*U_CoReg(:,:,k_i)';
        
    end    
    
    %% Evaluate Loss
    loss(itr) = Loss(U_CoReg,L_CoReg);
    
    %% Check Convergence
    %%% Check Loss Change
    if itr > 1 && abs(loss(end-1)-loss(end)) < epsilon
        exitinfo.reason = 'converge';
        break;
    end
    
    %%% Exceed Max Iteration
    if itr > MaxItr
        exitinfo.reason = 'timeout';
        break;
    end
    
    %% Display Results
    if itr == 1
        string = sprintf('loss = %.5f; change of loss = %.5f\n',loss(itr),0);
        % fprintf('%s',string);
    else
        %         for i = 1:length(string)
        %             fprintf('\b');
        %         end
        string = sprintf('loss = %.5f; change of loss = %.5f\n',loss(itr),loss(itr)-loss(itr-1));
        % fprintf('%s',string);
    end
    itr = itr+1;

end

%% Normalize embeddings to output the right reduced Laplacian trace

U_All = [];
for k_i = 1:size(U_CoReg,3)
    U_All = [U_All func_L2Normalize(U_CoReg(:,:,k_i),2)];
end

% Normalize
U_All_Normalized = func_L2Normalize(U_All,2);

for k_i = 1:nKernel
    U_Coreg_Normalized = U_All_Normalized(:,(k_i-1)*nMotion+1:k_i*nMotion);
    RedLapTrace = RedLapTrace + trace(U_Coreg_Normalized'*L_CoReg(:,:,k_i)*U_Coreg_Normalized);
end

function L = Loss(U_CoReg,L_CoReg)

L = 0;

for k_i = 1:size(U_CoReg,3)
    
    L = L + trace(U_CoReg(:,:,k_i)'*L_CoReg(:,:,k_i)*U_CoReg(:,:,k_i));
    
end





