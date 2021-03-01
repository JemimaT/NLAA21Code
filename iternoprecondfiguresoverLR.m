figure
hold on
for hno = 1:4
    subplot1 = subplot(2,2,hno);
    [itno,stepvec,LR] = iternoprecondoverLR(hno,1,1,0.01,1e-10); % (hno,SOARvalB,SOARvalR,inc,tol)

    plot((stepvec),(itno(1,:)),'Color', 'k','LineStyle','-','LineWidth',2,'DisplayName','$$L_B = 0.1$$')%,'DisplayName',sprintf('L_R = %f', LRvec(1)) );
    plot((stepvec),(itno(2,:)),'Color', 'r','LineStyle',':','LineWidth',2,'DisplayName','$$L_B = 0.4$$')%,'DisplayName',sprintf('L_R = %f', LRvec(2)));
    plot((stepvec),(itno(3,:)),'Color', 'b','LineStyle','-','LineWidth',2,'DisplayName','$$L_B = 0.7$$')%,'DisplayName',sprintf('L_R = %f', LRvec(3)));
    plot((stepvec),(itno(4,:)),'Color', 'c','LineStyle','--','LineWidth',2,'DisplayName','$$L_B = 1$$')%,'DisplayName',sprintf('L_R = %f', LRvec(3)));
  
    ylabel('Iterations','Interpreter','latex')
    set(subplot1,'FontSize',14,'TickLabelInterpreter','latex');    
   
end
subplot1 = subplot(2,2,1);
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.185610363361342 0.601921025404417 0.146859192451508 0.224464824629851],...
    'Interpreter','latex');
title('$$H_1$$','Interpreter','latex')
xlabel({'$$L_R$$';'$$(a)$$'},'Interpreter','latex')
subplot(2,2,2);
title('$$H_2$$','Interpreter','latex')
xlabel({'$$L_R$$';'$$(b)$$'},'Interpreter','latex')
subplot(2,2,3);
title('$$H_3$$','Interpreter','latex')
xlabel({'$$L_R$$';'$$(c)$$'},'Interpreter','latex')
subplot(2,2,4);
title('$$H_4$$','Interpreter','latex')
xlabel({'$$L_R$$';'$$(d)$$'},'Interpreter','latex')