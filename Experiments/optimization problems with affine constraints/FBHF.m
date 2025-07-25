clear all
close all
clc


N=2000;
m=N/2;
k=100;
ub=1;
lb=0;
Ncomp=20;
itlim = 100000;


tol=1e-6;

rng('default') 




kappav = [0.8 0.9 0.999];
factorv = [0.9 0.95];
flambdav = [0.9 0.8 0.7];


for i=1:Ncomp
    i
    A=randn(m,N);
    D=randn(k,N);
    e=randn(k,1);
    b=randn(m,1);
    
    beta=1/(norm(A)^2);
    L=norm(D);
    % 
    AA = A'*A;
    Ab = A'*b;
    DT = D';
    MD = -D;


    %------ FBHF
    x0=zeros(N,1);
    u0=zeros(k,1);
    kappa1 = 0.999;
    ep = 2/(1+sqrt(1+16*beta^2*L^2));
    chi = 2*beta*ep;
    gam=kappa1*chi;
    error=1;
    
    r=0;

    tic
    while(error>tol) && r < itlim
        B1x=AA*x0-Ab;
        B2x=DT*u0;
        B2u=MD*x0;
        gx=x0-gam*(B1x+B2x);
        y=max(lb,min(ub,gx));
        gu=u0-gam*B2u;
        w=max(0,gu);
        B2y=DT*w;
        B2w=MD*y;
        x1=y-gam*(B2y-B2x);
        u1=w-gam*(B2w-B2u);
        error=sqrt((norm(x1-x0)^2+norm(u1-u0)^2)/(norm(x0)^2+norm(u0)^2));
        x0=x1;
        u0=u1;
        r=r+1;
    end
    time_FBHF(i)=toc;
    iter_FBHF(i)=r;
    f_FBHF(i)=norm(A*x0-b)^2/2+ones(k,1)'*abs(D*x0);


    % %------ FBHF inertial
    factor = [0.8 0.9 0.999];
    for j=1:3

        kappa1 = 1;
        
        ep = factor(j)*2/(1+sqrt(1+16*beta^2*L^2));
        chi = 2*beta*ep;
        gam = kappa1*chi; 
        psi = (2-ep)/(1+(L*gam)^2);
        lambda = 1;
        alpha = 0.9999*2*(psi/lambda-1)/((2*psi/lambda-1)+sqrt(8*psi/lambda-7));
        alpha1 = alpha+1;
        %(1-alpha)^2*(2-ep)-((1-alpha)^2+alpha*(1+alpha))*lambda*(1+L^2*gam^2);
        error=1;
        r=0;
        p = zeros(N,1);
        q = zeros(k,1);
        x0 = zeros(N,1);
        u0 = zeros(k,1);
        tic
        while(error>tol) && r < itlim
            r=r+1;
            B1x=AA*p-Ab;
            B2x=DT*q;
            B2u=MD*p;
            gx=p-gam*(B1x+B2x);
            y=max(lb,min(ub,gx));
            gu=q-gam*B2u;
            w=max(0,gu);
            B2y=DT*w;
            B2w=MD*y;
            x1=y-gam*(B2y-B2x);
            u1=w-gam*(B2w-B2u);
            x10 = x1-x0;
            u10 = u1-u0;
            error=sqrt((norm(x10)^2+norm(u10)^2)/(norm(x0)^2+norm(u0)^2)); 
            p = x1 + alpha*x10;
            q = u1 + alpha*u10;
            x0=x1;
            u0=u1; 
            
        end
        time_FBHFI(i,j)=toc;
        iter_FBHFI(i,j)=r;
        f_FBHFI(i,j)=norm(A*x0-b)^2/2+ones(k,1)'*abs(D*x0);
    end

    % %------ FBHF inertial D1
    factor = 0.999;
    kappa1 = 1;
    
    ep = factor*2/(1+sqrt(1+16*beta^2*L^2));
    chi = 2*beta*ep;
    gam = kappa1*chi; %arriba 3.4553e-04
    psi = (2-ep)/(1+(L*gam)^2);
    lambda = 1;
    alphalim = 0.9999*2*(psi/lambda-1)/((2*psi/lambda-1)+sqrt(8*psi/lambda-7));

    %(1-alpha)^2*(2-ep)-((1-alpha)^2+alpha*(1+alpha))*lambda*(1+L^2*gam^2);
    error=1;
    r=0;
        p = zeros(N,1);
        q = zeros(k,1);
        x0 = zeros(N,1);
        u0 = zeros(k,1);
    tic
    while(error>tol) && r < itlim
        r=r+1;
        B1x=AA*p-Ab;
        B2x=DT*q;
        B2u=MD*p;
        gx=p-gam*(B1x+B2x);
        y=max(lb,min(ub,gx));
        gu=q-gam*B2u;
        w=max(0,gu);
        B2y=DT*w;
        B2w=MD*y;
        x1=y-gam*(B2y-B2x);
        u1=w-gam*(B2w-B2u);
        x10 = x1-x0;
        u10 = u1-u0;
        error=sqrt((norm(x10)^2+norm(u10)^2)/(norm(x0)^2+norm(u0)^2)); 
        alpha = 1/(1+0.001*r*(log(r)).^1.001); 
        p = x1 + alpha*x10;
        q = u1 + alpha*u10;
        x0=x1;
        u0=u1; 
    end
    time_FBHFDI1(i)=toc;
    iter_FBHFDI1(i)=r;
    f_FBHFDI1(i)=norm(A*x0-b)^2/2+ones(k,1)'*abs(D*x0);


   % %------ FBHF inertial D2
    factor = 0.999;
    kappa1 = 1;
    
    ep = factor*2/(1+sqrt(1+16*beta^2*L^2));
    chi = 2*beta*ep;
    gam = kappa1*chi; %arriba 3.4553e-04
    psi = (2-ep)/(1+(L*gam)^2);
    lambda = 1;
    %alphalim = 0.9999*2*(psi/lambda-1)/((2*psi/lambda-1)+sqrt(8*psi/lambda-7));

    error=1;
    r=0;
        p = zeros(N,1);
        q = zeros(k,1);
        x0 = zeros(N,1);
        u0 = zeros(k,1);
    tic
    while(error>tol) && r < itlim
        r=r+1;
        B1x=AA*p-Ab;
        B2x=DT*q;
        B2u=MD*p;
        gx=p-gam*(B1x+B2x);
        y=max(lb,min(ub,gx));
        gu=q-gam*B2u;
        w=max(0,gu);
        B2y=DT*w;
        B2w=MD*y;
        x1=y-gam*(B2y-B2x);
        u1=w-gam*(B2w-B2u);
        x10 = x1-x0;
        u10 = u1-u0;
        error=sqrt((norm(x10)^2+norm(u10)^2)/(norm(x0)^2+norm(u0)^2)); 
        alpha = 1/(3+0.00001*r*(log(r)).^1.00001);
        p = x1 + alpha*x10;
        q = u1 + alpha*u10;
        x0=x1;
        u0=u1; 
    end
    time_FBHFDI2(i)=toc;
    iter_FBHFDI2(i)=r;
    f_FBHFDI2(i)=norm(A*x0-b)^2/2+ones(k,1)'*abs(D*x0);

    % %------ FBHF inertial D3
    factor = 0.999;
    kappa1 = 1;
    
    ep = factor*2/(1+sqrt(1+16*beta^2*L^2));
    chi = 2*beta*ep;
    gam = kappa1*chi; %arriba 3.4553e-04
    psi = (2-ep)/(1+(L*gam)^2);
    lambda = 1;
    %alphalim = 0.9999*2*(psi/lambda-1)/((2*psi/lambda-1)+sqrt(8*psi/lambda-7));

    error=1;
    r=0;
        p = zeros(N,1);
        q = zeros(k,1);
        x0 = zeros(N,1);
        u0 = zeros(k,1);
    tic
    while(error>tol) && r < itlim
        r=r+1;
        B1x=AA*p-Ab;
        B2x=DT*q;
        B2u=MD*p;
        gx=p-gam*(B1x+B2x);
        y=max(lb,min(ub,gx));
        gu=q-gam*B2u;
        w=max(0,gu);
        B2y=DT*w;
        B2w=MD*y;
        x1=y-gam*(B2y-B2x);
        u1=w-gam*(B2w-B2u);
        x10 = x1-x0;
        u10 = u1-u0;
        error=sqrt((norm(x10)^2+norm(u10)^2)/(norm(x0)^2+norm(u0)^2)); 
        alpha = 1/(9+0.00001*r*(log(r)).^1.00001);
        p = x1 + alpha*x10;
        q = u1 + alpha*u10;
        x0=x1;
        u0=u1; 
    end
    time_FBHFDI3(i)=toc;
    iter_FBHFDI3(i)=r;
    f_FBHFDI3(i)=norm(A*x0-b)^2/2+ones(k,1)'*abs(D*x0);

    %------ FBHF relaxed inertial

    kappalam = [0.95 0.999];
    factor = 0.9;
    ep = factor*2/(1+sqrt(1+16*beta^2*L^2));
    chi = 2*beta*ep;
    gam=0.999*chi;
    psi = (2-ep)/(1+(L*gam)^2);
    for j=1:2
        r0 = zeros(N,1);
        s0 = zeros(k,1);
        p0 = r0;
        q0 = s0;

        lambda = kappalam(j)*psi;
        lambda1 = 1-lambda;
        alpha = 0.9999*2*(psi/lambda-1)/((2*psi/lambda-1)+sqrt(8*psi/lambda-7));
        (1-alpha)^2*(2-ep)-((1-alpha)^2+alpha*(1+alpha))*lambda*(1+L^2*gam^2);

        error=1;
        r=0;

        tic
        while(error>tol) && r < itlim
            
            B1x=AA*r0-Ab;
            B2x=DT*s0;
            B2u=MD*r0;
            gx=r0-gam*(B1x+B2x);
            y=max(lb,min(ub,gx));
            gu=s0-gam*B2u;
            w=max(0,gu);
            B2y=DT*w;
            B2w=MD*y;
            x1=y-gam*(B2y-B2x);
            u1=w-gam*(B2w-B2u);
            p = lambda*x1 + lambda1*r0;
            q = lambda*u1 + lambda1*s0;
            pp0 = p-p0;
            qq0 = q-q0;
            error=sqrt((norm(pp0)^2+norm(qq0)^2)/(norm(p0)^2+norm(q0)^2));
            r0 = p + alpha*(pp0);
            s0 = q + alpha*(qq0);
            p0 = p;
            q0 = q;

            r=r+1;
        end
        time_FBHFRI(i,j)=toc;
        iter_FBHFRI(i,j)=r;
        f_FBHFRI(i,j)=norm(A*x0-b)^2/2+ones(k,1)'*abs(D*x0);
    end
end

