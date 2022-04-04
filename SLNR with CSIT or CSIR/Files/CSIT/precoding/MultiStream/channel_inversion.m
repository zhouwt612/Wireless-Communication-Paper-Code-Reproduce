function u = channel_inversion(N,S,K,H)

u=zeros(N,K*S);
u=H'*inv(H*H');

for i1=1:K*S,
    u(:,i1)=u(:,i1)/(sqrt((u(:,i1)')*u(:,i1)));
end