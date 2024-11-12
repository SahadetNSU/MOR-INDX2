
points = 101;
freqlow = 10^(-6);
freqhig = 10^(6);
%linear space 
linst = (freqhig-freqlow)/(points-1);
freq = freqlow:linst:freqhig;
s =freq;

freqpts =max(size(s));
for q = 1:freqpts,
    
        HL_ori = C*((exp(sqrt(-1)*s(q))*E-A)\B);
      Hl = C_r*((exp(sqrt(-1)*s(q))*E_r-A_r)\B_r);
    
        nrmHL_ori(q)=norm(full(HL_ori));
        nrmHl(q)=norm(full(Hl));
       err(q)=norm(HL_ori-Hl);
         
end
% 
figure(6);  % norms in linear sclae
semilogx(s,nrmHL_ori,'r:o','linewidth',2);
hold on
semilogx(s, nrmHl,'b-*','linewidth',1);
hold off

xlabel('Frequency (rad/sec)')
ylabel('Magnitude')
legend('Original System','Reduced System')
%         
figure(7); 
semilogx(s,err,'b', s, ones(1,freqpts)*errb,'r--');
% semilogy(s,err,'b', s, errb,'r--');
xlabel('Frequency (rad/sec)')
ylabel('Absolute error')
legend('Error of System','Error Bound')


