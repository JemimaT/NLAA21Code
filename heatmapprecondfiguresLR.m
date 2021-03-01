figure
axes1 = axes;
hold(axes1,'on');
for hno = 1:4
    [stepvec,LRvec,condST] = heatmapprecond(100,200,1,1,hno,1,0.01); % (p,N,SOARvalB,SOARvalR,hno,condtype,inc)
    % where SOARval == 1 is SOAR, ==2 is Laplacian
    % condtype == 1 is just save condition number ==2 also saves bounds
    subplot1 = subplot(2,2,hno);
    hold(subplot1,'on');
    condS = condST.';
    logcondS=log10(condS);
    
    hold on
    h=pcolor(LRvec,stepvec,log10(condS));
    set(h,'edgecolor','none')
    [C,h]=contour(LRvec,stepvec,log10(condS),[0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.25 4.5 4.75 5],'k');
    axis('xy')
    caxis([1,5])
    c=colorbar;
    ylabel(c,'$$\log_{10}(\kappa(S))$$','Interpreter','latex')
   
    ylabel('$$L_B$$','Interpreter','latex')
    set(subplot1,'FontSize',14,'TickLabelInterpreter','latex');    
end



%plot([1,90],[1,90],'r')
subplot1 = subplot(2,2,1);
title('$$H_1$$','Interpreter','latex')
xlabel({'$$L_R$$';'$$(a)$$'},'Interpreter','latex')%';'there'})
subplot2 =subplot(2,2,2);
title('$$H_2$$','Interpreter','latex')
xlabel({'$$L_R$$';'$$(b)$$'},'Interpreter','latex')
subplot3 =subplot(2,2,3);
title('$$H_3$$','Interpreter','latex')
xlabel({'$$L_R$$';'$$(c)$$'},'Interpreter','latex')
subplot4 =subplot(2,2,4);
title('$$H_4$$','Interpreter','latex')
xlabel({'$$L_R$$';'$$(d)$$'},'Interpreter','latex')
