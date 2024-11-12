clc; clear all; 
close all;


%%%%% initialling the data %%%%%

n=20; m_k=2; p_k=3;

c1 = 0.1; c2 = 0.6; c3 = 0.9; 
s1 = 0.5; s2 = 0.3; s3 = 0.8;
t1 = 0.7; t2 = 0.9; t3 = 0.9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Evaluation of parameters and matrices for different time-step

E1 = eye(18,18);

n_f=size(E1,1);
E = blkdiag(E1, zeros(n-n_f,n-n_f));

A11 =  0.1*[     
                 -c1   c2   c3   0   0   0   0   0   0;...
                 0    0    0   s2  -s3  s1   0   0   0;...
                 0    0    0   0   0   0   t3 t1  t2;... 
                  s1   -s2   s3   0   0   0   0   0   0;...
                  0    0    0   -t2  t3  t1   0   0   0;...
                  0    0    0   0   0   0   c3 -c1  c2;...
                  t1   t2  -t3   0   0   0   0   0   0;...
                  0    0    0   c2  c3  -c1   0   0   0;...
                  0    0    0   0   0   0  -s3  s1 -s2 ];
                  
 A1 =   [A11 A11'; A11' A11];
 
 A12 = [c1  -s2  -t3  s1  -t2  -c3  t1  -c2  -s3;...
        -s3  t1  c2   s2  -t3  -c1  -s3  t1   c2]';
   
A2 = [A12 ;  A12];
A3 = A2';
A4 = zeros(n-n_f,n-n_f);
A = [A1, A2; A3, A4];




B = [4 -1 s3+1 1 0 -2 0 1 s1+1 -s1+2 ...
     1 0 s1+1 -2 1 -1 0 -13 0 0; ...
     1 -1 s3 1 0 -2 0 1 s1 2 ...
     1 -1 s2 -2 1 -2 0 13 0 0]';
 
  B1 = B(1:18,:) ;   
  B2 = B(19:20,:) ;  
%  Bk(:,:,i)=[4 -1 s3+1 1 0 -2 0 1 0 0;...
%             1 0 s1+1 -2 1 -1 0 -13 0 0]';
     
     
 C = [ 1 0 0 1 0 0 s1+1 -s3+1 0 0 ...
       0 1 0 0 1 0 s2+1 -s1+2 0 0; ...
       0 0 1 0 0 .05+c1 s3+1 -s2+2 0 0 ...
       1 0 0 1 1 0 -s1 -s3-1 0 0; ...
       0 1 0 0 1 0 -s2+1 -s1-2 0 0 ...
       0 0 1 0 -1 .05+c1 -s3 -s2 0 0];
   
 C1 = C(:,1:18);   
 C2 = C(:,19:20);
%  Ck(:,:,i)=[0 0 1 0 0 1 0 0 0 0;...
%             2 0 0 0 0 0 1 0 0 0;...
%             0 0 0 0 0 0 0 .05+c1 0 0];
     
%Pl = (E1 - (A2*inv(A2'*inv(E1)*A2)*A2'*inv(E1)));
Pl=eye(size(E1,1)) - A2*((A2'*(E1\A2))\(A2'/E1));

  [QL,sig,QR]=svd(Pl,0);
  
  Ql = QL(:,1:16);
  Qr = QR(:,1:16);
  

E_n = Qr'*E1*Qr;

A_n = Qr'*A1*Qr;

B_n = Qr'*B1;

C_n = C1*Qr;


E_m = Pl*E1*Pl';

A_m = Pl*A1*Pl';

B_m = Pl*B1;

C_m = C1*Pl';


%  Eig_lift=eig(A_n,E_n)
%  Eig_lift_ori=eig(A,E)    
%  
%  figure(1)
% plot(Eig_lift,'*');
% hold on 
% plot(Eig_lift_ori,'ro');
% xlabel('J')
% ylabel('Magnitude of eigenvalues')
% hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


     