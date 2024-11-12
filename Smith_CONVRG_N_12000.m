
%%%%%%%%%%%%%%Plot ADI convergence History%%%%%%%%%%%%%%%%%%%%%%%

figure(3);
%ns=niter_1
ns=20
J=1:ns;
semilogy(J,res1_cc(1:ns),'--rs','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',10); 
xlabel('normalized residual norm')
ylabel('number of iterations')
legend('Controllability Gramian' )                   
 title('Smith convergence history')
hold on
%figure(2);
%nr=niter_2
nr=20
J=1:nr;
semilogy(J,res1_oc(1:nr),'--rs','LineWidth',2,...
                       'MarkerEdgeColor','r',...
                       'MarkerFaceColor','m',...
                       'MarkerSize',10); 
ylabel('normalized residual norm')
xlabel('number of iterations')
legend('Observability Gramian')                   
% title('ADI convergence history: type C ')
hold off