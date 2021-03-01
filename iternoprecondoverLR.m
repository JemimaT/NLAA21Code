%Script to make Figure 3 b, d and f

%Want to produce no iterations to convergence for incremented Lr and Lb 

%[iterno,condS,diff]=conjgrad(p,Lr,Lb,hno,tol,maxit)
%Both B and R are SOAR matrices - want to see how the matrix condition
%number and bounds change as we change the lengthscale of one and fix the
%other

%Most simple case: let background error variance sig_b^2=1
%H ones on diagonal zero everywhere else
%R and of fixed lengthscale and radius
%Vary B lengthscale

%Rtest=SOAR(50,0.1,1);
function [itno,stepvec,LBvec] = iternoprecondoverLR(hno,SOARvalB,SOARvalR,inc,tol)
%Rtest=eye(100);
hold on
%Set variables
p=100;N=2*p;

% 
% if hno==1 || hno==4
%     tol=1e-6;
% else
%     tol=1e-6;
% end
%tol=1e-6;

maxit=3000;

%a=0.1/(2*sin(pi/500));
a=1;%a=l/(2*pi);
noRlines = 3;
%angle = 3*pi/(4);
%lengthscale=2*a*sin(angle/2);%3*(a/pi*sin(pi/p));
%SELECT H OBSERVATION OPERATOR:

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
%Construct truth, w
delx = 0.1;
q = 2*pi/(N*delx);
w= zeros([N,1]);
s=[];
for k=1:N
    s(k)=(k-1)*delx;
    w(k,1)= 2*sin(q*s(k))+ cos(3*q*s(k))-0.3*sin(125*q*s(k));
end

HtH=Htest.'*Htest;
HHt=Htest*Htest.';
hhmax=eigs(Htest*Htest.',1);
hhmin=eigs(Htest*Htest.',1,'sm');


stepvec=[];itno=[];
relresvec=[];outdiffvec=[];flagvec=[];k=1;yminvec=[];ymaxvec=[];bmaxvec=[];
%figure
%hold on
for Bstep =0.1:0.3:1
    %delxp=(a*sin(pi/p));
    %lengthscale = Rstep*pi/6;
    %% Think we want to change lengthscale in terms of angle - set L_R = 2asin(phi/2) 
    % where phi is the angle of the circle we want to cover with our
    % correlation lengthscale.
    
    LBvec(k)=Bstep;
    if SOARvalB ==1
        [Btest,Binv,Bvecs,Bvals]=SOARinv(N,Bstep,a);
    elseif SOARvalB ==2
        [Btest,Binv,Bvecs,Bvals]=Laplacian(N,Bstep,a);
    end    %ymax=real(ymax);ymin=real(ymin);
    %Rinv=Rtest\eye(p);
    i=1;
    for step=0.05:inc:1
        %     %If B is SOAR matrix
         if SOARvalR ==1
            [Rtest,Rinv,Rvecs,Rvals]=SOARinv(p,step,a);%Blengthscale,a);
        elseif SOARvalB ==2
            [Rtest,Rinv,Rvecs,Rvals] = Laplacian(p,step,a);
        end
        %LBvec(k)=step;
        %matrix square root of B
        Bsq = sqrtm(Btest);%Bvecs*diag(sqrt(Bvals))*ctranspose(Bvecs);
        %pxp version
        %NxN version
        Bfirst = Bsq*Htest.'*Rinv*Htest*Bsq;
        S=(eye(N)+Bfirst);%eigs(eye(p)+Rfirst,1);
        %Sunsort = Cinv+productunsort;
        b=S*w;
        
        [xout,flagno,relres,iterno,xvec] = pcg(S,b,tol,maxit);
        relresvec(k,i)=relres;
        flagvec(k,i)=flagno;
        %plot(resvec)
        xoutM(i,:)=xout;
        outdiffvec(k,i) = norm(w-xout);
        itno(k,i)=iterno;
        stepvec(i)=step;
        i=i+1;
    end
    k=k+1;
end
%figure
%plot(stepvec,itno,'Color', 'b','LineStyle','-','LineWidth',3)
%plot(stepvec,outdiffvec,'Color', 'k','LineStyle','-','LineWidth',3)
%figure

end