figure
hold on

LBvec = [0.1,0.5,1];
for index = 0:2
    for hno = 1:4
    subplot1 = subplot(3,4,hno+index*4);
    precondboundtestLBinc(hno,1,1,LBvec(index+1),100,0.01) %SOARval = 1 gives SOAR corr, = 2 gives Laplacian
    set(subplot1,'FontSize',12,'TickLabelInterpreter','latex');    
   
    end
end
subplot(3,4,1)
ylabel({'$$L_B = 0.1$$';'$$\kappa(S)$$'},'Interpreter','latex')
xlabel('$$(a)$$','Interpreter','latex')
title('$$H_1$$','Interpreter','latex')
subplot(3,4,2)
title('$$H_2$$','Interpreter','latex')
xlabel('$$(b)$$','Interpreter','latex')
subplot(3,4,3)
title('$$H_3$$','Interpreter','latex')
xlabel('$$(c)$$','Interpreter','latex')
subplot(3,4,4)
title('$$H_4$$','Interpreter','latex')
xlabel('$$(d)$$','Interpreter','latex')
subplot(3,4,5)
ylabel({'$$L_B = 0.5$$';'$$\kappa(S)$$'},'Interpreter','latex')
xlabel('$$(e)$$','Interpreter','latex')
subplot(3,4,6)
xlabel('$$(f)$$','Interpreter','latex')
subplot(3,4,7)
xlabel('$$(g)$$','Interpreter','latex')

subplot(3,4,8)
xlabel('$$(h)$$','Interpreter','latex')
subplot(3,4,9)
ylabel({'$$L_B = 1$$';'$$\kappa(S)$$'},'Interpreter','latex')
xlabel({'$$L_R$$';'$$(i)$$'},'Interpreter','latex')

subplot(3,4,10)
xlabel({'$$L_R$$';'$$(j)$$'},'Interpreter','latex')
subplot(3,4,11)
xlabel({'$$L_R$$';'$$(k)$$'},'Interpreter','latex')
subplot(3,4,12)
xlabel({'$$L_R$$';'$$(l)$$'},'Interpreter','latex')