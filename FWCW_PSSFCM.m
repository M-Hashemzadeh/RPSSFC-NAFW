function [Cluster_elem,M,Z]=FWCW_PSSFCM(X,M,k,t_max,N,fuzzy_degree,d,f,b,balance_parm,a_coefficient,b_coefficient,T_pow, alpha_1, alpha_2, landa, Neig, NR, p, q)
%
%[Cluster_elem,M,W,Z]=FWCW_FCM(X,M,k,p_init,p_max,p_step,t_max,beta_memory,N,fuzzy_degree,d,q,landa)
%
%A.Najafi, N.Aghazadeh, M.Hashemzadeh and A.Golzari oskouei, "A Robust Possibilistic Semi-Supervised Fuzzy
%Clustering Algorithm with Neighborhood-Aware Feature Weighting"
%
%Function Inputs
%===============
%
%X is an Nxd data matrix, where each row corresponds to an instance.
%
%M is a kxd matrix of the initial cluster centers. Each row corresponds to a center.
%
%k is the number of clusters.
%
%t_max is the maximum number of iterations.
%
% Other parameters
%
%Function Outputs
%================
%
%Cluster_elem is a kxd matrix containing the final cluster assignments.
%
%M is a kxd matrix of the final cluster centers. Each row corresponds to a center.
%
%z is a kxd matrix of the final weights of each fatuter in each cluster.
%
%Courtesy of A. Golzari Oskouei

%--------------------------------------------------------------------------
%Weights are uniformly initialized.
Z=ones(k,d)/d;  %initial faeture weights
w=ones(N,NR)/NR;  %initial faeture weights
Cluster_elem = ones(k,N)/k;

%Other initializations.
Iter=1; %Number of iterations.
O_F_old=inf; %Previous iteration objective (used to check convergence).
%--------------------------------------------------------------------------

fprintf('\nStart of fuzzy C-means clustering method based on feature-weight and cluster-weight learning iterations\n');
fprintf('----------------------------------\n\n');

%The proposed iterative procedure.
while 1
    
    %Update the cluster assignments.
    for j=1:k
        distance(j,:,:) = (1-tanh((-1.*repmat(landa,N,1).*((X-repmat(M(j,:),N,1)).^2))));
        WBETA = transpose(Z(j,:).^q);
        WBETA(WBETA==inf)=0;
        dNK(:,j) = reshape(distance(j,:,:),[N,d]) * WBETA   ;
        
        cc = (1-Cluster_elem(j,:)).^fuzzy_degree;
        dNKprim = dNK(:,j)';
        dNK_neig(:,j)= ((a_coefficient+alpha_2).*dNK(:,j)) + (alpha_1/NR).*sum(   (cc(Neig).*dNKprim(Neig)) .* w.^p   ,2);
        
    end
    
    tmp1 = zeros(N,k);
    tmp3 = zeros(N,k);
    tmp5 = zeros(N,k);
    for j=1:k
        tmp2 = (dNK_neig./repmat(dNK_neig(:,j),1,k)).^(1/(fuzzy_degree-1));
        tmp2(tmp2==inf)=0;
        tmp2(isnan(tmp2))=0;
        tmp1=tmp1+tmp2;
        
        tmp4 = ((alpha_2.*b.*f.*dNK)./repmat(dNK_neig(:,j),1,k)).^(1/(fuzzy_degree-1));
        tmp4(tmp4==inf)=0;
        tmp4(isnan(tmp4))=0;
        tmp3=tmp3+tmp4;
        
        tmp6(:,j) = ((alpha_2.*b.*f(:,j).*dNK(:,j))./dNK_neig(:,j)).^(1/(fuzzy_degree-1));
        tmp6(tmp6==inf)=0;
        tmp6(isnan(tmp6))=0;
        tmp5(:,j)=tmp5(:,j)+tmp6(:,j);
    end
    Cluster_elem = ((1+tmp3-tmp5)./tmp1)';
    
    
    %Update gama.
    for j=1:k
        distance(j,:,:) = (1-tanh((-1.*repmat(landa,N,1).*((X-repmat(M(j,:),N,1)).^2))));
        WBETA = transpose(Z(j,:).^q);
        WBETA(WBETA==inf)=0;
        dNK(:,j) = reshape(distance(j,:,:),[N,d]) * WBETA ;
    end
    
    gama = balance_parm * 2 * sum(dNK .* transpose(Cluster_elem.^fuzzy_degree)) ./ sum(Cluster_elem.^fuzzy_degree,2)';
    
    %Update the T.
    T = (1+(repmat(2*b_coefficient./gama, N,1) .* dNK).^(1/(T_pow-1))).^-1;
    
    %Calculate the fuzzy C-means clustering method based on feature-weight and cluster-weight learning objective.
    O_F=object_fun_FWCW_PSSFCM(N,d,k,Cluster_elem,landa,M,fuzzy_degree,Z,X,gama,T,a_coefficient,b_coefficient,T_pow,alpha_1, alpha_2,b,f,Neig,NR,w,p, q);
    
    if ~isnan(O_F)
        fprintf('The clustering objective function is %f\n\n',O_F);
    end
    
    %Check for convergence. Never converge if in the current (or previous)
    %iteration empty or singleton clusters were detected.
    if Iter>=t_max || ~isnan(O_F) && ~isnan(O_F_old) && (abs(1-O_F/O_F_old) < 1e-6 )
        
        fprintf('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n');
        fprintf('The final objective function is =%f.\n',O_F);
        fprintf('\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n');
        
        break;
        
    end
    
    O_F_old=O_F;
    
    %Update the cluster centers.
    mf1 = Cluster_elem.^fuzzy_degree;       % MF matrix after tanhonential modification
    tf = T'.^T_pow;
    mf2 = (Cluster_elem-(b.*f)').^fuzzy_degree;
    mf = (a_coefficient.*mf1) + (b_coefficient.*tf) + (alpha_2 * mf2);
    
    for j=1:k
        cc = (1-Cluster_elem(j,:)).^fuzzy_degree;
        a1  = (X .* (sech((-1.*repmat(landa,N,1)).*((X-repmat(M(j,:),N,1)).^2))).^2);
        b1  = sum((a1(Neig) .* cc(Neig)) .* w.^p,2);
        
        a2  = ((sech((-1.*repmat(landa,N,1)).*((X-repmat(M(j,:),N,1)).^2))).^2);
        b2  = sum((a2(Neig) .* cc(Neig)) .* w.^p,2);
        
        M(j,:) = ((mf(j,:) * (X .* (sech((-1.*repmat(landa,N,1)).*((X-repmat(M(j,:),N,1)).^2))).^2)) + ((alpha_1/NR) * Cluster_elem(j,:) .^fuzzy_degree *b1) )./ ((((mf(j,:)*(sech((-1.*repmat(landa,N,1)).*((X-repmat(M(j,:),N,1)).^2))).^2)))+ ((alpha_1/NR) * Cluster_elem(j,:) .^fuzzy_degree *b2) ); %new center
    end
    
    %Update the feature weights.
    for j=1:k
        distance(j,:,:) = (1-tanh((-1.*repmat(landa,N,1).*((X-repmat(M(j,:),N,1)).^2))));
        a2 = (1-tanh((-1.*repmat(landa,N,1).*((X-repmat(M(j,:),N,1)).^2))));
        b2  = sum(a2(Neig) .* cc(Neig) .* w.^p ,2);
        dWkm(j,:) = ((((a_coefficient.*(Cluster_elem(j,:).^fuzzy_degree)) + (b_coefficient.*(T(:,j).^T_pow))') + (alpha_2 * ((Cluster_elem(j,:)-(b.*f(:,j))').^fuzzy_degree)))* reshape(distance(j,:,:),[N,d]))  +   (((alpha_1/NR) * Cluster_elem(j,:) .^fuzzy_degree *b2));
    end
    
    tmp1 = zeros(k,d);
    for j=1:d
        tmp2 = (dWkm./repmat(dWkm(:,j),1,d)) .^ (1/(q-1));
        tmp2(tmp2==inf)=0;
        tmp2(isnan(tmp2))=0;
        tmp1=tmp1+tmp2;
    end
    Z = 1./tmp1;
    Z(isnan(Z))=1;
    Z(Z==inf)=1;
    
    if nnz(dWkm==0)>0
        for j=1:k
            if nnz(dWkm(j,:)==0)>0
                Z(j,find(dWkm(j,:)==0)) = 1/nnz(dWkm(j,:)==0);
                Z(j,find(dWkm(j,:)~=0)) = 0;
            end
        end
    end
    
    
    
    %Update the sample weights.
    for j=1:k
        a2 = (1-tanh((-1.*repmat(landa,N,1).*((X-repmat(M(j,:),N,1)).^2))));
        WBETA = transpose(Z(j,:).^q);
        WBETA(WBETA==inf)=0;
        b2  = (a2 * WBETA);
        mf = (Cluster_elem(j,:)).^fuzzy_degree;
        cc = (1-Cluster_elem(j,:)).^fuzzy_degree;
        dWknr(j,:,:) = mf' .* b2(Neig) .* cc(Neig);
    end
    dWnr = reshape(sum(dWknr),[N,NR]);
    
    tmp1 = zeros(N,NR);
    for j=1:NR
        tmp2 = (dWnr./repmat(dWnr(:,j),1,NR)).^(1/(p-1));
        tmp2(tmp2==inf)=0;
        tmp2(isnan(tmp2))=0;
        tmp1=tmp1+tmp2;
    end
    w = 1./tmp1;
    w(isnan(w))=1;
    w(w==inf)=1;
    
    if nnz(dWnr==0)>0
        for j=1:N
            if nnz(dWnr(j,:)==0)>0
                w(j,find(dWnr(j,:)==0)) = 1/nnz(dWnr(j,:)==0);
                w(j,find(dWnr(j,:)~=0)) = 0;
            end
        end
    end
    
    Iter=Iter+1;
end
end



