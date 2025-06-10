function j_fun = object_fun_FWCW_PSSFCM(N,d,k,Cluster_elem,landa,M,fuzzy_degree,Z,X,gama,T,a_coefficient,b_coefficient,T_pow,alpha_1, alpha_2,b,f,Neig,NR, w, p, q)

j_fun3 = 0;

for j=1:k
    distance(j,:,:) = (1-tanh((-1.*repmat(landa,N,1).*((X-repmat(M(j,:),N,1)).^2))));
    WBETA = transpose(Z(j,:).^q);
    WBETA(WBETA==inf)=0;
    dNK(:,j) = reshape(distance(j,:,:),[N,d]) * WBETA ;
    
    %term 2
    cc = (1-Cluster_elem(j,:)).^fuzzy_degree;
    dNKprim = dNK(:,j)';
    j_fun3=j_fun3+sum((Cluster_elem(j,:).^fuzzy_degree)'.* sum(   (cc(Neig).*dNKprim(Neig)) .* w.^p   ,2));
    
end

j_fun1 = 2 * sum(sum(dNK .* ((a_coefficient.*transpose(Cluster_elem.^fuzzy_degree)) + (b_coefficient.*(T.^T_pow)))));
j_fun2 = sum(sum(dNK .* transpose((Cluster_elem-(b.*f)').^fuzzy_degree)));
j_fun4 = sum(sum((1-T).^T_pow).* gama);
j_fun = j_fun1 +(alpha_2 * j_fun2)+(alpha_1/NR)*j_fun3 + j_fun4 ;

end

