function precondboundtestLBinc(hno,SOARvalB,SOARvalR,LB,p,inc)
N=2*p;
a=1; %dimension of system
%Choose H
Htest=zeros([p,N]);Htestunsort=zeros([p,N]);
if hno==1
    Htest = eye([p,N]);
elseif hno==2
    for k=1:p
        Htest(k,2*k-1)=1;
    end
elseif hno==3
    for k=1:p-2
        Htest(k+1, 2*k-1:2*k+3)=1/5;
    end
    Htest(1,1:3)=1/5; Htest(1,N-1:N)=1/5;
    Htest(p,1)=1/5; Htest(p,N-3:N)=1/5;
elseif hno==4
    %data=(1:N);
    %vec = datasample(data,p,'Replace',false);
    %vec = sort(vec);
%     vec=randperm(N);
%     vecsort = sort(vec(1:p));
%     for k=1:p
%         Htest(k,vecsort(k))=1;
%         Htestunsort(k,vec(k))=1;
%     end
    Htest = HtestRANDOM;
end
hhmax = eigs(Htest.'*Htest,1);%hterm in bounds
hhmin = eigs(Htest*Htest.',1,'sm');

stepvec=[];condS=[];
lower=[];upper=[];upperhaben=[];
lowerhaben=[];
%figure
%hold on


%change LR as a number (easier to compare)
%2*a*sin(angle/2);

%generate R matrix and other relevant terms
if SOARvalB ==1
    [Btest,Binv,Bvecs,Bvals]=SOARinv(N,LB,a);
elseif SOARvalR ==2
    [Btest,Binv,Bvecs,Bvals]=Laplacian(N,LB,a);
end
%matrix square root of R
Bsq = sqrtm(Btest);%Rvecs*diag(1./sqrt(Rvals))*ctranspose(Rvecs);
m=1;

for step=0.01:inc:1
    if SOARvalR ==1
        [R,Rinv,Rvecs,Rvals]=SOARinv(p,step,a);%Blengthscale,a);
    elseif SOARvalB ==2
        [R,Rinv,Revecs,Rvals] = Laplacian(p,step,a);
    end
    
    %LBvec(k)=step;
    %matrix square root of B
    Rsq = sqrtm(Rinv);%Bvecs*diag(sqrt(Bvals))*ctranspose(Bvecs);
    %pxp version
    Rfirst = Rsq*Htest*Btest*Htest.'*Rsq;
    %NxN version
    Bfirst = Bsq*Htest.'*Rinv*Htest*Bsq;
    condS(m)=cond(eye(N)+Bfirst);%eigs(eye(p)+Rfirst,1);
    
    %calculate condition number using cond
    upperhaben(m)=1+norm(Rfirst,inf);
    %infinity norm in NxN space
    %upperhabenBfirst(k,m)=1+norm(Bfirst,inf);
    % my lower bound
    lower(m) = 1+min(Bvals)/min(Rvals)*hhmin;
    lower2(m) = 1+min(Bvals)*hhmax/max(Rvals);
    %my upper bound
    upper(m) = 1+max(Bvals)/min(Rvals)*hhmax;
    %haben's lower bound - sum of entries
    lowerhaben(m) = 1+1/p*sum(sum(Rfirst));
   
    stepvec(m)=step;

    m=m+1;
end

hold on
plot(stepvec,condS,'k')
plot(stepvec,lowerhaben,'b--')
plot(stepvec,upperhaben,'b--')
plot(stepvec,lower,'r-.')
%plot(stepvec,lower2,'c')
plot(stepvec,upper,'r-.')
set(gca, 'YScale', 'log')


end