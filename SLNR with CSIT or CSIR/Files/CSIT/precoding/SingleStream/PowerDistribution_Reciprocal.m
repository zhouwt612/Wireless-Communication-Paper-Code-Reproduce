% reference is a row vector
function p = PowerDistribution_Reciprocal(reference)

[row_num, col_num]=size(reference);

p=zeros(1,col_num);

for i1=1:col_num,
    p(1,i1)=col_num*(1/reference(1,i1))/sum(1./reference);
end