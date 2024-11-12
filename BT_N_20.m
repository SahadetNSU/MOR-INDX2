
errb=0;
%compute the SVD 
  
[Uc,sig,Vc]=svd((L{iter1})'*E1* (R{iter1}));   
tol1=10^(-2); %MOR tolerance
for l0=1:size(sig',1) 
   if (sig(l0,l0)>tol1) 
        rank_0c=l0;
    end;
end
UC0=Uc(:,1:rank_0c);
 VC0=Vc(:,1:rank_0c);
  ss0=sig(1:rank_0c,1:rank_0c);
  ss_0=sig(rank_0c+1:end,rank_0c+1:end);
  errb=errb+2*sum(diag(ss_0));
  
  
  S0=[real(L{iter1})*UC0*inv(sqrt(ss0))]; 
  T0=[real(R{iter1})*VC0*inv(sqrt(ss0))];
   
   Er1 = S0'*E1*T0;
   Ar1 = S0'*A1*T0;
   Br1 = S0'*B1;
   Cr1 = C1*T0;

% main_eig = eig(full(A),full(E));
% proj_eig = eig(full(A_n),full(E_n));
%red_eig = eig(full(Ar1),full(Er1));

% figure(1); clf;
% plot(main_eig,'ko','linewidth',1)
% hold on
% % plot(proj_eig,'r*')
% % hold on
 %plot(red_eig,'*')
% hold off
% % legend('Original System','Projected System','Reduced System')
% legend('Original System','Reduced System')

figure(2); clf;
semilogy(1:size(ss0,1),diag(ss0),'b*');
xlabel('j');
ylabel('Proper Hankel singular values')
legend('Hankel Singular Values');
  
pHsv=sort(diag(blkdiag(ss0))); %causal Hankel singular values
pHsv_trunc=sort(diag(blkdiag(ss0,ss_0)));

figure(3); clf;
semilogy(size(pHsv,1):-1:1,pHsv,'b*','linewidth',1); 
hold on
semilogy(size(pHsv_trunc,1):-1:1,pHsv_trunc,'r-o');
hold off
xlabel('j');
ylabel('Proper Hankel singular values')
legend('Truncated Hankel Singular Values', 'All HSV');

% clear A1 A2 A3 A4 ans B1 B2 C1 C2 E1 E2 E3 E4;
% clear Uc sig Vc UC0 VC0 ss0 ss_0 S0 T0 pHsv pHsv_trunc main_eig proj_eig red_eig;
%    