function u = MMSE(N,Mi,S,K,H,sigma2)

u=zeros(N,K*S);
u=H'*inv(H*H'+sigma2*eye(K*Mi));

for i1=1:K,
    u(:,i1)=u(:,i1)/(sqrt((u(:,i1)')*u(:,i1)));
end