function codebook = RVQ(MT,R,L)
% construct a codebook which contains L codes of R*MT dimension

codebook = zeros(L*R,MT);

for i1=1:L,
    temp1=(1/sqrt(2))*complex(randn(R,MT),randn(R,MT));
    temp1=temp1/(sqrt(trace(temp1'*temp1)));
    
    codebook((i1-1)*R+1:i1*R,:)=temp1;
end