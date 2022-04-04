function A = gram_schmidt(B)

[row_num, col_num]=size(B);
A=zeros(row_num,col_num);

A(:,1)=B(:,1);
for i1=2:col_num,
    temp=zeros(row_num,1);
    for i2=1:(i1-1),
        temp=temp+conj((B(:,i1)'*A(:,i2))/(A(:,i2)'*A(:,i2)))*A(:,i2);
    end
    A(:,i1)=B(:,i1)-temp;
end

for i1=1:col_num,
    A(:,i1)=A(:,i1)/sqrt(A(:,i1)'*A(:,i1));
end