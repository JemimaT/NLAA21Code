%Preconditioning tests - R and B are both SOAR, usual choices of H1-H4
function [stepvec,LRvec,condS]=heatmapprecond(p,N,SOARvalB,SOARvalR,hno,condtype,inc)
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
%     %data=(1:N);
%     %vec = datasample(data,p,'Replace',false);
%     %vec = sort(vec);
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

i=1;out=[];stepvec=[];condS=[];
k=1;LRvec=[];lower=[];upper=[];upperhaben=[];
lowerhaben=[];
%figure
%hold on
for Rstep =0.01:inc:1%:noRlines
    m=1;
    %change LR as a number (easier to compare)
    %2*a*sin(angle/2);
    LRvec(k)=Rstep;
    %generate R matrix and other relevant terms
    if SOARvalR ==1
        [Rtest,Rinv,Rvecs,Rvals]=SOARinv(p,Rstep,a);
    elseif SOARvalR ==2
        [Rtest,Rinv,Rvecs,Rvals]=Laplacian(p,Rstep,a);
    end
    %matrix square root of R
    Rsq = sqrtm(Rinv);%Rvecs*diag(1./sqrt(Rvals))*ctranspose(Rvecs);
  
    
    for step=0.01:inc:1
        if SOARvalB ==1
            [C,Cinv,Cvecs,Cvals]=SOARinv(N,step,a);%Blengthscale,a);
        elseif SOARvalB ==2
            [C,Cinv,evecs,Cvals] = Laplacian(N,step,a);
        end
        lmax=max(diag(Cvals));lmin=min(diag(Cvals));

        %LBvec(k)=step;
        %matrix square root of B
        Bsq = sqrtm(C);%Bvecs*diag(sqrt(Bvals))*ctranspose(Bvecs);
        %pxp version
        Rfirst = Rsq*Htest*C*Htest.'*Rsq;
        %NxN version
        Bfirst = Bsq*Htest.'*Rinv*Htest*Bsq;
        condS(k,m)=cond(eye(N)+Bfirst);%eigs(eye(p)+Rfirst,1);
        if condtype ==2
            %calculate condition number using cond
            upperhaben(k,m)=1+norm(Rfirst,inf);
            %infinity norm in NxN space
            %upperhabenBfirst(k,m)=1+norm(Bfirst,inf);
            % my lower bound
            lower(k,m) = 1+min(Cvals)/min(Rvals)*hhmin;
            lower2(k,m) = 1+min(Cvals)*hhmax/max(Rvals);
            %my upper bound
            upper(k,m) = 1+max(Cvals)/min(Rvals)*hhmax;
            %haben's lower bound - sum of entries
            lowerhaben(k,m) = 1+1/p*sum(sum(Rfirst));
        end
        stepvec(m)=step;

        i=i+1;
        m=m+1;
    end
    k=k+1;
end

end