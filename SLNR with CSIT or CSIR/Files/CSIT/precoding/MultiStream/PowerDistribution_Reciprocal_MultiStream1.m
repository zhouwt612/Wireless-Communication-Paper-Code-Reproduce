% reference is a row vector
% For the multistream case, the least value in reference is used
function p = PowerDistribution_Reciprocal_MultiStream1(K,S,reference)

p=zeros(1,K*S);

temp=zeros(1,K);
for i1=1:K,
    temp(1,i1)=min(reference(1,(i1-1)*S+1:i1*S));
end

p_temp=zeros(1,K);
p_temp=PowerDistribution_Reciprocal(temp);

for i1=1:K,
    for i2=1:S,
        p(1,(i1-1)*S+i2)=p_temp(1,i1);
    end
end