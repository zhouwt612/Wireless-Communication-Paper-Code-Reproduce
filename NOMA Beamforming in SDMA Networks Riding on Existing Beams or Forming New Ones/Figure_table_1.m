clear all
%close all

Rp=5;
P0=0.5;
sigma = 10^(-4);
K=2; %number of users
N=K;
A=complex(sqrt(0.5)*randn(N,K),sqrt(0.5)*randn(N,K));
[x1,x2,x3]=svd(A);
H =  complex(sqrt(0.5)*randn(N,K),sqrt(0.5)*randn(N,K)); %N by K
% H=[  -0.0723 - 0.6116i   0.2257 - 0.1166i
%   -0.1707 - 0.0212i   0.2212 + 0.4439i];
H = [1 0; 0 1]+0.0001* complex(sqrt(0.5)*randn(N,K),sqrt(0.5)*randn(N,K));
Px = H*inv(H'*H);
P = Px/trace(Px*Px');
g =  [sin(pi*45/1/180) cos(pi*45/1/180)]'+ 0.0001*complex(sqrt(0.5)*randn(N,1),sqrt(0.5)*randn(N,1));;%H(:,1)+0.0001*complex(sqrt(0.5)*randn(N,1),sqrt(0.5)*randn(N,1)); 
%[sin(pi*45/1/180) cos(pi*45/1/180)]'

G = g*g'; 

for k = 1 : K
    tau(k) = max(0,abs(H(:,k)'*P(:,k))^2/(2^Rp-1)-sigma);
end

gP = g'*P*P'*g;  
%Strategy one
%Case I
cvx_begin quiet sdp 
   cvx_precision best
   variable W(N,N) hermitian
   dual variables c{K+2} 
   maximize( trace(W*G) )
   subject to
      for i = 1 : K
          real(trace(W*H(:,i)*H(:,i)')) <=tau(i) : c{i};
      end
      c{K+1} : real(trace(W)) <=P0;
      c{K+2} : W>=0;
cvx_end

svd_case1 = svd(W);
[U,di,V]=svd(W);
w_case1 = U(:,1)*sqrt(di(1,1));
rate_case1 = log2(1+abs(g'*w_case1)^2/(gP+sigma));

%Case II - contains two cases
a0=real(gP+sigma);
for i = 1 : K
    ai(i) = real(H(:,i)'*P*P'*H(:,i)+sigma);
end 

cvx_begin quiet sdp 
   cvx_precision best
   variable z(1,1) nonnegative
   variable W(N,N) hermitian
   dual variables c{K+3}
   maximize( z )
   subject to
      c{1} : real(trace(W*G)) >=a0*z;
      c{2} : real(trace(H(:,1)*H(:,1)'*W)) >=z*ai(1);
      for i = 2 : K
          real(trace(W*H(:,i)*H(:,i)')) <=tau(i) : c{i+1};
      end
      c{K+2} : real(trace(W)) <=P0;
      c{K+3} : W>=0;
cvx_end
svd_case21 = svd(W);
[U,di,V]=svd(W);
w_case21 = U(:,1)*sqrt(di(1,1));
rate_case21 = log2(z+1);  

%Case II - second  
a0=real(gP+sigma);
for i = 1 : K
    ai(i) = real(H(:,i)'*P*P'*H(:,i)+sigma);
end 

cvx_begin quiet sdp 
   cvx_precision best
   variable z(1,1) nonnegative
   variable W(N,N) hermitian
   dual variables c{K+3}
   maximize( z )
   subject to
      c{1} : real(trace(W*G)) >=a0*z;
      c{2} : real(trace(H(:,2)*H(:,2)'*W)) >=z*ai(2);
      for i = 1 : 1
          real(trace(W*H(:,i)*H(:,i)')) <=tau(i) : c{i+2};
      end
      c{K+2} : real(trace(W)) <=P0;
      c{K+3} : W>=0;
cvx_end
svd_case22 = svd(W);
[U,di,V]=svd(W);
w_case22 = U(:,1)*sqrt(di(1,1));
rate_case22 = log2(z+1);  


%Case III
cvx_begin quiet sdp 
   cvx_precision best
   variable z(1,1) nonnegative
   variable W(N,N) hermitian
   dual variables c{K+3}
   maximize( z )
   subject to
      c{1} : real(trace(W*G)) >=a0*z;
      for i = 1 : K
          c{i+1} : real(trace(H(:,i)*H(:,i)'*W)) >= z*ai(i);
      end
      c{K+2} : real(trace(W)) <=P0;
      c{K+3} : W>=0;
cvx_end 

svd_case3 = svd(W);
[U,di,V]=svd(W);
w_case3 = U(:,1)*sqrt(di(1,1));
rate_case3 = log2(z+1);  

Rate1=real([rate_case1 rate_case21 rate_case22   rate_case3  ])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Strategy 2
%Case I

for i = 1 : K
    alpha1(i) = min([P0/(P(:,i)'*P(:,i)) tau(i)/abs(H(:,i)'*P(:,i))^2]);
end
rx_case1 = log2(1+alpha1.*abs(g'*P).^2/(gP+sigma));
 

for i = 1 : K
    alpha2 = P0/(P(:,i)'*P(:,i));
    rx_case2(i) = min(log2(1+alpha2*abs(g'*P(:,i))^2/(sigma+(g'*P*P'*g))) ,...
        log2(1+alpha2*abs(H(:,i)'*P(:,i))^2/(sigma+H(:,i)'*P*P'*H(:,i))));
end 

real([rx_case1; rx_case2])
[svd_case1 svd_case21 svd_case22 svd_case3]

 


