%What is going on with structure?
% Fix LB=0.5 and vary LR. Why is convergence unchanged for LB >0.4 when the condition number changes 

p=100; N=2*p;hno=2 ;LB=0.5;H1 = eye([p,N]);
Htestunsort=zeros([p,N]);
H2=zeros([p,N]);H3=zeros([p,N]);H4=zeros([p,N]);
for k=1:p
    H2(k,2*k-1)=1;
end

for k=1:p-2
    H3(k+1, 2*k-1:2*k+3)=1/5;
end
H3(1,1:3)=1/5; H3(1,N-1:N)=1/5;
Ht3(p,1)=1/5; H3(p,N-3:N)=1/5;

%data=(1:N);
%vec = datasample(data,p,'Replace',false);
%vec = sort(vec);
H4=HtestRANDOM;
%     vec=randperm(N);
%     vecsort = sort(vec(1:p));
%     for k=1:p
%         Htest(k,vecsort(k))=1;
%         Htestunsort(k,vec(k))=1;
%     end


[R1,R1inv,R1vecs,R1vals]= SOARinv(p,0.1,1); [R2,R2inv,R2vecs,R2vals] = SOARinv(p,0.4,1); 
[R3,R3inv,R3vecs,R3vals]= SOARinv(p,0.7,1);
[B1,B1inv,B1vecs,B1vals] = SOARinv(N,LB,1);
[B2,B2inv,B2vecs,B2vals] = SOARinv(N,0.5,1);
[B3,B3inv,B3vecs,B3vals] = SOARinv(N,0.7,1);
[R4,R4inv,R4vecs,R4vals] = SOARinv(p,1,1);
%Bsq = sqrtm(B);




prodvalslow = real(eig(H3*B2*H3.'*R1inv));
prodvals1 = real(eig(H1*B1*H1.'*R1inv));
prodvals2 = real(eig(H1*B1*H1.'*R2inv));
prodvals3 = real(eig(H1*B1*H1.'*R3inv));
prodvals4 = real(eig(H1*B1*H1.'*R4inv));
figure
subplot(2,2,1)
semilogx((prodvals1),ones([1,p]),'bx')
hold on
semilogx(prodvals2,1.2*ones([1,p]),'ro')
semilogx(prodvals3,1.4*ones([1,p]),'kd')
semilogx(prodvals4,1.6*ones([1,p]),'cs')
xlim([min(prodvalslow),max(prodvals4)])
ylim([0.8,1.8])
yticks([1 1.2 1.4 1.6])
yticklabels({'L_R = 0.1','L_R = 0.4','L_R = 0.7','L_R = 1'})
title('H_1')
subplot(2,2,2)

prodvals1_2 = real(eig(H2*B1*H2.'*R1inv));
prodvals2_2 = real(eig(H2*B1*H2.'*R2inv));
prodvals3_2 = real(eig(H2*B1*H2.'*R3inv));
prodvals4_2 = real(eig(H2*B1*H2.'*R4inv));


semilogx((prodvals1_2),ones([1,p]),'bx')
hold on
semilogx(prodvals2_2,1.2*ones([1,p]),'ro')
semilogx(prodvals3_2,1.4*ones([1,p]),'kd')
semilogx(prodvals4_2,1.6*ones([1,p]),'cs')
xlim([min(prodvalslow),max(prodvals4)])
yticks([1 1.2 1.4 1.6])
yticklabels({'L_R = 0.1','L_R = 0.4','L_R = 0.7','L_R = 1'})
ylim([0.8,1.8])
title('H_2')
subplot(2,2,3)

prodvals1_3 = real(eig(H3*B1*H3.'*R1inv));
prodvals2_3 = real(eig(H3*B1*H3.'*R2inv));
prodvals3_3 = real(eig(H3*B1*H3.'*R3inv));
prodvals4_3 = real(eig(H3*B1*H3.'*R4inv));


semilogx((prodvals1_3),ones([1,p]),'bx')
hold on
semilogx(prodvals2_3,1.2*ones([1,p]),'ro')
semilogx(prodvals3_3,1.4*ones([1,p]),'kd')
semilogx(prodvals4_3,1.6*ones([1,p]),'cs')
xlim([min(prodvalslow),max(prodvals4)])
yticks([1 1.2 1.4 1.6])
yticklabels({'L_R = 0.1','L_R = 0.4','L_R = 0.7','L_R = 1'})
ylim([0.8,1.8])
title('H_3')
subplot(2,2,4)

prodvals1_4 = real(eig(H4*B1*H4.'*R1inv));
prodvals2_4 = real(eig(H4*B1*H4.'*R2inv));
prodvals3_4 = real(eig(H4*B1*H4.'*R3inv));
prodvals4_4 = real(eig(H4*B1*H4.'*R4inv));


semilogx((prodvals1_4),ones([1,p]),'bx')
hold on
semilogx(prodvals2_4,1.2*ones([1,p]),'ro')
semilogx(prodvals3_4,1.4*ones([1,p]),'kd')
semilogx(prodvals4_4,1.6*ones([1,p]),'cs')
yticks([1 1.2 1.4 1.6])
yticklabels({'L_R = 0.1','L_R = 0.4','L_R = 0.7','L_R = 1'})
ylim([0.8,1.8])
xlim([min(prodvalslow),max(prodvals4)])
title('H_4')
%%
% figure
% semilogx((prodvals1),ones([1,p]),'bx')
% hold on
% semilogx(prodvals2,1.2*ones([1,p]),'rx')
% semilogx(prodvals3,1.4*ones([1,p]),'kx')
% semilogx(prodvals4,1.6*ones([1,p]),'cx')
% 
% semilogx((prodvals1_2),1.05*ones([1,p]),'bo')
% 
% semilogx(prodvals2_2,1.25*ones([1,p]),'ro')
% semilogx(prodvals3_2,1.45*ones([1,p]),'ko')
% semilogx(prodvals4_2,1.65*ones([1,p]),'co')
% 
% semilogx((prodvals1_3),1.1*ones([1,p]),'bd')
% 
% semilogx(prodvals2_3,1.3*ones([1,p]),'rd')
% semilogx(prodvals3_3,1.5*ones([1,p]),'kd')
% semilogx(prodvals4_3,1.7*ones([1,p]),'cd')
% 
% semilogx((prodvals1_4),1.15*ones([1,p]),'bs')
% semilogx(prodvals2_4,1.35*ones([1,p]),'rs')
% semilogx(prodvals3_4,1.55*ones([1,p]),'ks')
% semilogx(prodvals4_4,1.75*ones([1,p]),'cs')