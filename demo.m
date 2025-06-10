% This demo shows how to call the Neighborhood based Sample and Feature Weighted PFCM algorithm described in the paper:
% A.Najafi, N.Aghazadeh, M.Hashemzadeh and A.Golzari oskouei,
% "A Robust Possibilistic Semi-Supervised Fuzzy Clustering Algorithm with Neighborhood-Aware Feature Weighting"
% For the demonstration, the letter dataset of the above paper is used.

clc
clear all
close all

%Load the dataset. The last column of dataset is true labels.
X=load('iris.mat');
X=X.iris;

class=X(:,end);
X(:,end)=[];    %delete last column (true labels) in clustering process

[N,d]=size(X);
X=(X(:,:)-min(X(:)))./(max(X(:)-min(X(:)))); %Normalize data between 0 and 1 (optinal)

%Algorithm parameters.
%---------------------
k=size(unique(class),1);  %number of clusters.
t_max=100;                %maximum number of iterations.
Restarts=10;              %number of algorithm restarts.
fuzzy_degree=2;           %fuzzy membership degree
q = 2;
balance_parm = 1;         %balance parameter among to terms of loss function
T_pow = 2;                %power of T
a_coefficient = 1;        %coefficient of u
b_coefficient = 1;        %coefficient of T

alpha_1 = 4;
alpha_2 = 1;
p = 2;

I = 0.0001;     %The value of this parameter is in the range of (0 and 1]

landa=I./var(X);
landa(landa==inf)=1;

NR = 7;         % size of window for finding neighbors
Neig = Find_Neighbors(NR, X);

f = double(class == 1:max(class));% convert class to onehot encoding
f(:,sum(f)==0)=[];
labeled_rate = 20;                % rate of labeled data (0-100)

% Variable to store the execution time of the first iteration
elapsed_time_first_repeat = 0;
%---------------------
%Cluster the instances using the propsed procedure.
%---------------------------------------------------------
for repeat=1:Restarts
    fprintf('========================================================\n')
    fprintf('proposed clustering algorithm: Restart %d\n',repeat);
    
    % label indicator vector
    rand('state',repeat)
    b = zeros(N,1);
    tmp1=randperm(N);
    b(tmp1(1:N*labeled_rate/100))=1;
    
    %initialize with labeled data.
    if labeled_rate==0
        %Randomly initialize the cluster centers.
        tmp2=randperm(N);
        M=X(tmp2(1:k),:);
    else
        M = ((b.*f)'*X)./repmat(sum(b.*f)',1,d);
        if sum(isnan(M))>=1
            tmp2=randperm(N);
            tem3= X(tmp2(1:k),:);
            M(isnan(M))=tem3(isnan(M));
        end
    end
    
    %Execute proposed clustering algorithm.
    %Get the cluster assignments, the cluster centers and the cluster weight and feature weight.
    time1 = clock;
    [Cluster_elem,M,Z]=FWCW_PSSFCM(X,M,k,t_max,N,fuzzy_degree,d, f,b,balance_parm,a_coefficient,b_coefficient,T_pow, alpha_1, alpha_2, landa, Neig, NR, p, q);
    time2 = clock;
    simtime(repeat) = etime(time2, time1);
    
    [~,unsupervised_Cluster]=max(Cluster_elem,[],1); %Hard clusters. Select the largest value for each sample among the clusters, and assign that sample to that cluster.
    semisupervised_Cluster  = unsupervised_Cluster;
    
    semisupervised_Cluster(tmp1(1:N*labeled_rate/100))=class(tmp1(1:N*labeled_rate/100));
    
    % Evaluation metrics
    
    EVAL = Evaluate(class,semisupervised_Cluster');
    Accurcy_semisupervised(repeat)=EVAL(1);
    F1_scores_semisupervised(repeat) = EVAL(2);
    NMI_semisupervised(repeat)=EVAL(3);
    
    fprintf('End of Restart %d\n',repeat);
    fprintf('========================================================\n\n')
end

fprintf('Average semisupervised accurcy over %d restarts: %f.\n',Restarts,mean(Accurcy_semisupervised));
fprintf('Average F1_score over %d restarts: %f.\n', Restarts, mean(F1_scores_semisupervised));
fprintf('Average semisupervised NMI over %d restarts: %f.\n',Restarts,mean(NMI_semisupervised));
fprintf('Average Execution time for the first repeat: %f seconds.\n', mean(simtime));




