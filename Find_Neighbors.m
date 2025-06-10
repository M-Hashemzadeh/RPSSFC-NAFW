function result = Find_Neighbors(NR, X)

    function D2 = naneucdist(XI,XJ)
        %NANEUCDIST Euclidean distance ignoring coordinates with NaNs
        n = size(XI,2);
        sqdx = (XI-XJ).^2;
        nstar = sum(sqdx,2); % Number of pairs that do not contain NaNs
        D2 = 1- tanh(-1.*nstar);
    end

D = pdist(X,@naneucdist);
D = squareform(D);
[~,result] = sort(D, 2);
result = result(:,2:NR+1);

end