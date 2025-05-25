function [ A1Complete, U, delta] = covTransformerDecomposer(ZIGMAtemp)
%This function does eigenvalue decomposition

[V,HH,W] = eig(ZIGMAtemp);
[h,ind] = sort(diag(HH),'descend');
delta = HH(ind,ind);
U = V(:,ind);

A1Complete = U*delta.^(1/2);

end

