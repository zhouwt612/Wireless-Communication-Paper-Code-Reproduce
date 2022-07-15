clear all
%close all

Rp=1;
disc = 3;
al=3;
sigma2_dbm= -94;%+10*log10(BW)+Nf; %Thermal noise in dBm -90+10+10
sigma=10^((sigma2_dbm-30)/10);
P0vec = [10 20 30 40];
Kvec = [2  6 11 16];
for ip = 1: length(Kvec)
    K = Kvec(ip);
    P0=10^((P0vec(3)-30)/10);% 1.5;
    N=K;
    
    for ict = 1 : 5000
        H = complex(sqrt(0.5)*randn(N,K),sqrt(0.5)*randn(N,K));
        H = H;
        Px = H*inv(H'*H);
        P = Px/sqrt(trace(Px*Px'))*sqrt(10^((20-30)/10));%P_sdma = 20 dBm
        locm = sign(rand(K,2)-0.5) .* disc.*rand(K,2); %location  of users
        dist = max(1, sqrt(sum(locm.^2,2))); %distance to the base station
        H = H* diag(1./dist.^(al/2));
        g = complex(sqrt(0.5)*randn(N,1),sqrt(0.5)*randn(N,1));;%H(:,1)+0.0001*complex(sqrt(0.5)*randn(N,1),sqrt(0.5)*randn(N,1)); 

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
        
        heff = real(diag(H'*H));
        [mx,mxid] = max(heff);
        for ii = mxid
            cvx_begin quiet sdp 
               cvx_precision best
               variable z(1,1) nonnegative
               variable W(N,N) hermitian
               dual variables c{K+4}
               maximize( z )
               subject to
                  c{1} : real(trace(W*G)) >=a0*z;
                  c{2} : real(trace(H(:,ii)*H(:,ii)'*W)) >=z*ai(ii);
                  jj=1;
                  for j = 1 : K                         
                      if ~mod(j,ii)
                          continue
                      end
                      jj = jj+1;
                      real(trace(W*H(:,j)*H(:,j)')) <=tau(j) : c{jj+1};
                  end
                  c{K+3} : real(trace(W)) <=P0;
                  c{K+4} : W>=0;
            cvx_end
            svd_case21 = svd(W);
            [U,di,V]=svd(W);
            w_case21 = U(:,1)*sqrt(di(1,1));
            rate_case21 = log2(z+1);  
        end

        rate_typex1(ict)=max(real([rate_case1 rate_case21]));

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

        rate_typex2(ict)=max(max(real([rx_case1; rx_case2])));
    end
    
    rate_type1(ip) = mean(rate_typex1);
    rate_type2(ip) = mean(rate_typex2);
end

%plot( Kvec,rate_type2  )
plot(Kvec,rate_type1,Kvec,rate_type2  )

 


