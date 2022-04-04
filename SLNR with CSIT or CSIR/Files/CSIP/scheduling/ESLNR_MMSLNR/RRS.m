function scheduled = RRS(Kall,K,scheduled_old)

scheduled=zeros(1,K);

temp1=scheduled_old(1,K);

scheduled(1,1)=mod(temp1,Kall)+1;

for i1=2:K,
    scheduled(1,i1)=mod(scheduled(1,i1-1),Kall)+1;
end