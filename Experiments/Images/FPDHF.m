clc; 
clear all;
close all;

%Original Image
% z=imread('peppers2.tif'); %512
% z=imread('bird.tif'); %256
z=imread('saturn.tif'); % 128
[m,n]=size(z);
z=double(z)./max(double(z(:)));
zz=z;


rng('default');

%%parameters
lam1 = 1e-2; %1e-2;
lam2 = 1e-3; %1e-4;
lam3 = 1;


%Linear operator
% forward finite differences (with Neumann boundary conditions)
L1= @(x) [x(:,2:end,:)-x(:,1:end-1,:), zeros(size(x,1),1,size(x,3))]; % horizontal
L2 = @(x) [x(2:end,:,:)-x(1:end-1,:,:); zeros(1,size(x,2),size(x,3))]; % vertical
% backward finite differences (with Neumann boundary conditions)
L1T = @(x) [-x(:,1,:), x(:,1:end-2,:)-x(:,2:end-1,:), x(:,end-1,:)];    % horizontal
L2T = @(x) [-x(1,:,:); x(1:end-2,:,:)-x(2:end-1,:,:); x(end-1,:,:)];    % vertical

L=@(X) {L1(X),L2(X)};
 
LT=@(X) L1T(X{1})+L2T(X{2}) ;

%D = nabla d
psf = fspecial('average', 3); % Average Kernel 9x9
%psf = fspecial('average', 9); % Average Kernel 9x9
%psf = fspecial('gaussian', 3 , 1 ); % Gaussian 3x3
K = @(x) imfilter(x, psf,'symmetric');
KA = @(x) imfilter(x, rot90(psf,2),'symmetric'); 
    
numcases = 20; % Number of cases

% Blurred image
for j=1:numcases 
    b{j}=K(zz);
    % Noise
    b{j} = imnoise(b{j},'gaussian',0,1e-3);
end
%


%C wavelets
delta = 1e-2;
wav=@(X)wavelet('2D haar',3,X);
iwav=@(X)wavelet('2D haar',-3,X);
C = @(x) lam2*iwav(grad_Huber(wav(x),delta));


% parameters
nL=sqrt(8);
zeta = lam2/delta;
beta = 1/(lam3* sum(abs(psf(:))) ); % norma de K^2 ;


% iteration number
niter= 5000; %5000
tol = 1e-6;


factor = [0.5 0.75 0.999];
kappa1 = 0.29;
kappa2 = 0.99;

for j = 1 : numcases
    j
    Kb = KA(b{j});
    D = @(X)  lam3*(KA(K(X))-Kb);
    %% %%%% FPDHF
    % step-sizes
    ep = 0.999*2/(1+sqrt(1+16*beta^2*zeta^2));
    chi = 2*beta*ep; 
    t = kappa1*chi;
    s = kappa2*(1-t/chi)/(t*nL^2);
    
    lam1s=lam1/s;
    x = b{j};
    u{1}=zeros(m,n);
    u{2}=zeros(m,n);
    
    error=1;
    nit=0;
    r=0;
    
    errorFPDHF = 0.1;
    
    tic
    while errorFPDHF > tol && r < niter
        r=r+1;
            x_o=x;
            u_o=u;       
            p = C(x_o);
    
            z_= x_o - t*( LT(u_o) + p + D(x_o) );
            z = min(max(0,z_),1);
        
            q = t*(C(z)-p);
    
            uu = L( 2*z-x_o-q );
            
            u1 = (u_o{1}+s*uu{1})./s;
            prox = sign(u1) .* max(0, abs(u1) - lam1s);
            u{1} = (u1-prox)*s;
            
    
            u2 = (u_o{2}+s*uu{2})./s;
            prox = sign(u2) .* max(0, abs(u2) - lam1s);
            u{2} = (u2-prox)*s;
            
            x = z - q;
    
            errorFPDHF = sqrt((norm(x-x_o,'fro')^2+norm(u{1}-u_o{1},'fro')^2+norm(u{2}-u_o{2},'fro')^2) /(norm(x_o,'fro')^2 + norm(u{1},'fro')^2+norm(u{2},'fro')^2));
    
    end
    time_FPDHF(j) = toc;
    iter_FPDHF(j) = r;
    x_FPDHF{j}=x;
    xx = L(x);
    fo_FPDHF(j) = lam3*norm(K(x)-b{j},2)^2/2+lam1*sum(sum(abs(xx{1})))+lam1*sum(sum(abs(xx{2})))+huberf(wav(x),lam2,delta);
    SNR_FPDHF(j) = - 20*log10(norm(zz - x,2)/norm(zz,2));
    
    
    %% %%%% Inertial FPDHF
    for i = 1:3
        ep = factor(i)*2/(1+sqrt(1+16*beta^2*zeta^2));
        chi = 2*beta*ep; 
        t = kappa1*chi;
        s = kappa2*(1-t/chi)/(t*nL^2);
        zetat = t*zeta/(sqrt(1-s*t*nL^2));
        betat = beta*(1-s*t*nL^2);
        psi = (2-ep+2*zetat)/(1+zetat^2+2*zetat);
        lambda = 1;
        alpha = 0.9999*2*(psi/lambda-1)/((2*psi/lambda-1)+sqrt(8*psi/lambda-7));
        alpha1 = alpha+1;
        
       
        error=1;
        
        lam1s=lam1/s;
        
        x = b{j};
        u{1}=zeros(m,n);
        u{2}=zeros(m,n);
        x_o=x;
        u_o=u;
        p = alpha1*x-alpha*x_o;
        q={};
        q{1} = alpha1*u{1}-alpha*u_o{1};
        q{2} = alpha1*u{2}-alpha*u_o{2};
    
        
        r=0;
        
        errorFPDHFI = 1;
        tic
        while errorFPDHFI > tol && r < niter;
                r=r+1;
                x_o=x;
                u_o=u;
                
                Cp = C(p);
        
                z_= p - t*( LT(q) + Cp + D(p) );
                z = min(max(0,z_),1);
            
                x = z-t*(C(z)-Cp);
        
                uu = L( z + x -p );
                
                u1 = (q{1}+s*uu{1})./s;
                prox = sign(u1) .* max(0, abs(u1) - lam1s);
                u{1} = (u1-prox)*s;
                
        
                u2 = (q{2}+s*uu{2})./s;
                prox = sign(u2) .* max(0, abs(u2) - lam1s);
                u{2} = (u2-prox)*s;
                
    
                x2 =    x-x_o;
                u11 = u{1}-u_o{1};
                u12 = u{2}-u_o{2};
    
                errorFPDHFI = sqrt((norm(x2,'fro')^2+norm(u11,'fro')^2+norm(u12,'fro')^2) /(norm(x,'fro')^2 + norm(u_o{1},'fro')^2+norm(u_o{2},'fro')^2));
                
                p    =  x   +  alpha*x2;
                q{1} = u{1} +  alpha*u11;
                q{2} = u{2} +  alpha*u12;
    
        end
        time_FPDHFI(j,i) = toc;
        iter_FPDHFI(j,i) = r;
        x_FPDHFI{j,i}=x;
        xx = L(x);
        fo_FPDHFI(j,i) = lam3*norm(K(x)-b{j},2)^2/2+lam1*sum(sum(abs(xx{1})))+lam1*sum(sum(abs(xx{2})))+huberf(wav(x),lam2,delta);
        SNR_FPDHFI(j,i) = - 20*log10(norm(zz - x,2)/norm(zz,2));
    end
    
    
    %% %%%% Inertial FPDHFn 1
    for i = 1:1
        ep = 0.999*2/(1+sqrt(1+16*beta^2*zeta^2));
        chi = 2*beta*ep; 
        t = kappa1*chi;
        s = kappa2*(1-t/chi)/(t*nL^2);
        zetat = t*zeta/(sqrt(1-s*t*nL^2));
        betat = beta*(1-s*t*nL^2);
        psi = (2-ep+2*zetat)/(1+zetat^2+2*zetat);
        lambda = 1;
        alpha = 0.9999*2*(psi/lambda-1)/((2*psi/lambda-1)+sqrt(8*psi/lambda-7));
        alpha1 = alpha+1;
        
       
        error=1;
        
        lam1s=lam1/s;
        
        x = b{j};
        u{1}=zeros(m,n);
        u{2}=zeros(m,n);
        x_o=x;
        u_o=u;
        p = alpha1*x-alpha*x_o;
        q={};
        q{1} = alpha1*u{1}-alpha*u_o{1};
        q{2} = alpha1*u{2}-alpha*u_o{2};
    
        
        r=0;
        
        errorFPDHFI = 1;
        tic
        while errorFPDHFI > tol && r < niter;
                r=r+1;
                x_o=x;
                u_o=u;
                
                Cp = C(p);
        
                z_= p - t*( LT(q) + Cp + D(p) );
                z = min(max(0,z_),1);
            
                x = z-t*(C(z)-Cp);
        
                uu = L( z + x -p );
                
                u1 = (q{1}+s*uu{1})./s;
                prox = sign(u1) .* max(0, abs(u1) - lam1s);
                u{1} = (u1-prox)*s;
                
        
                u2 = (q{2}+s*uu{2})./s;
                prox = sign(u2) .* max(0, abs(u2) - lam1s);
                u{2} = (u2-prox)*s;
                
    
                x2 =    x-x_o;
                u11 = u{1}-u_o{1};
                u12 = u{2}-u_o{2};
    
                errorFPDHFI = sqrt((norm(x2,'fro')^2+norm(u11,'fro')^2+norm(u12,'fro')^2) /(norm(x,'fro')^2 + norm(u_o{1},'fro')^2+norm(u_o{2},'fro')^2));
                
                alpha = 1/(1+0.001*r*(log(r)).^1.001);
    
                p    =  x   +  alpha*x2;
                q{1} = u{1} +  alpha*u11;
                q{2} = u{2} +  alpha*u12;
    
        end
        time_FPDHFID1(j,i) = toc;
        iter_FPDHFID1(j,i) = r;
        x_FPDHFID1{j,i}=x;
        xx = L(x);
        fo_FPDHFID1(j,i) = lam3*norm(K(x)-b{j},2)^2/2+lam1*sum(sum(abs(xx{1})))+lam1*sum(sum(abs(xx{2})))+huberf(wav(x),lam2,delta);
        SNR_FPDHFID1(j,i) = - 20*log10(norm(zz - x,2)/norm(zz,2));
    end
    
    
    %% %%%% Inertial FPDHFn 2
    for i = 1:1
        ep = 0.999*2/(1+sqrt(1+16*beta^2*zeta^2));
        chi = 2*beta*ep; 
        t = kappa1*chi;
        s = kappa2*(1-t/chi)/(t*nL^2);
        zetat = t*zeta/(sqrt(1-s*t*nL^2));
        betat = beta*(1-s*t*nL^2);
        psi = (2-ep+2*zetat)/(1+zetat^2+2*zetat);
        lambda = 1;
        alpha = 0.9999*2*(psi/lambda-1)/((2*psi/lambda-1)+sqrt(8*psi/lambda-7));
        alpha1 = alpha+1;
        
       
        error=1;
        
        lam1s=lam1/s;
        
        x = b{j};
        u{1}=zeros(m,n);
        u{2}=zeros(m,n);
        x_o=x;
        u_o=u;
        p = alpha1*x-alpha*x_o;
        q={};
        q{1} = alpha1*u{1}-alpha*u_o{1};
        q{2} = alpha1*u{2}-alpha*u_o{2};
    
        
        r=0;
        
        errorFPDHFI = 1;
        tic
        while errorFPDHFI > tol && r < niter;
                r=r+1;
                x_o=x;
                u_o=u;
                
                Cp = C(p);
        
                z_= p - t*( LT(q) + Cp + D(p) );
                z = min(max(0,z_),1);
            
                x = z-t*(C(z)-Cp);
        
                uu = L( z + x -p );
                
                u1 = (q{1}+s*uu{1})./s;
                prox = sign(u1) .* max(0, abs(u1) - lam1s);
                u{1} = (u1-prox)*s;
                
        
                u2 = (q{2}+s*uu{2})./s;
                prox = sign(u2) .* max(0, abs(u2) - lam1s);
                u{2} = (u2-prox)*s;
                
    
                x2 =    x-x_o;
                u11 = u{1}-u_o{1};
                u12 = u{2}-u_o{2};
    
                errorFPDHFI = sqrt((norm(x2,'fro')^2+norm(u11,'fro')^2+norm(u12,'fro')^2) /(norm(x,'fro')^2 + norm(u_o{1},'fro')^2+norm(u_o{2},'fro')^2));
                
                alpha =  1/(3+0.00001*r*(log(r)).^1.00001);
    
                p    =  x   +  alpha*x2;
                q{1} = u{1} +  alpha*u11;
                q{2} = u{2} +  alpha*u12;
    
        end
        time_FPDHFID2(j,i) = toc;
        iter_FPDHFID2(j,i) = r;
        x_FPDHFID2{j,i}=x;
        xx = L(x);
        fo_FPDHFID2(j,i) = lam3*norm(K(x)-b{j},2)^2/2+lam1*sum(sum(abs(xx{1})))+lam1*sum(sum(abs(xx{2})))+huberf(wav(x),lam2,delta);
        SNR_FPDHFID2(j,i) = - 20*log10(norm(zz - x,2)/norm(zz,2));
    end
    
    %% %%%% Inertial FPDHFn 3
    for i = 1:1
        ep = 0.999*2/(1+sqrt(1+16*beta^2*zeta^2));
        chi = 2*beta*ep; 
        t = kappa1*chi;
        s = kappa2*(1-t/chi)/(t*nL^2);
        zetat = t*zeta/(sqrt(1-s*t*nL^2));
        betat = beta*(1-s*t*nL^2);
        psi = (2-ep+2*zetat)/(1+zetat^2+2*zetat);
        lambda = 1;
        alpha = 0.9999*2*(psi/lambda-1)/((2*psi/lambda-1)+sqrt(8*psi/lambda-7));
        alpha1 = alpha+1;
        
       
        error=1;
        
        lam1s=lam1/s;
        
        x = b{j};
        u{1}=zeros(m,n);
        u{2}=zeros(m,n);
        x_o=x;
        u_o=u;
        p = alpha1*x-alpha*x_o;
        q={};
        q{1} = alpha1*u{1}-alpha*u_o{1};
        q{2} = alpha1*u{2}-alpha*u_o{2};
    
        
        r=0;
        
        errorFPDHFI = 1;
        tic
        while errorFPDHFI > tol && r < niter;
                r=r+1;
                x_o=x;
                u_o=u;
                
                Cp = C(p);
        
                z_= p - t*( LT(q) + Cp + D(p) );
                z = min(max(0,z_),1);
            
                x = z-t*(C(z)-Cp);
        
                uu = L( z + x -p );
                
                u1 = (q{1}+s*uu{1})./s;
                prox = sign(u1) .* max(0, abs(u1) - lam1s);
                u{1} = (u1-prox)*s;
                
        
                u2 = (q{2}+s*uu{2})./s;
                prox = sign(u2) .* max(0, abs(u2) - lam1s);
                u{2} = (u2-prox)*s;
                
    
                x2 =    x-x_o;
                u11 = u{1}-u_o{1};
                u12 = u{2}-u_o{2};
    
                errorFPDHFI = sqrt((norm(x2,'fro')^2+norm(u11,'fro')^2+norm(u12,'fro')^2) /(norm(x,'fro')^2 + norm(u_o{1},'fro')^2+norm(u_o{2},'fro')^2));
                
                alpha =  1/(9+0.00001*r*(log(r)).^1.00001);
    
                p    =  x   +  alpha*x2;
                q{1} = u{1} +  alpha*u11;
                q{2} = u{2} +  alpha*u12;
    
        end
        time_FPDHFID3(j,i) = toc;
        iter_FPDHFID3(j,i) = r;
        x_FPDHFID3{j,i}=x;
        xx = L(x);
        fo_FPDHFID3(j,i) = lam3*norm(K(x)-b{j},2)^2/2+lam1*sum(sum(abs(xx{1})))+lam1*sum(sum(abs(xx{2})))+huberf(wav(x),lam2,delta);
        SNR_FPDHFID3(j,i) = - 20*log10(norm(zz - x,2)/norm(zz,2));
    end
    
    
    %% %%%% Relaxed Inertial FPDHF
    factorlam = [0.95 0.99];
    for i = 1:2
        ep = factor(3)*2/(1+sqrt(1+16*beta^2*zeta^2));
        chi = 2*beta*ep; 
        t = kappa1*chi;
        s = kappa2*(1-t/chi)/(t*nL^2);
        zetat = t*zeta/(sqrt(1-s*t*nL^2));
        betat = beta*(1-s*t*nL^2);
        psi = (2-ep+2*zetat)/(1+zetat^2+2*zetat);
        lambda = factorlam(i)*psi;
        lambda1 = 1-lambda;
        alpha = 0.9999*2*(psi/lambda-1)/((2*psi/lambda-1)+sqrt(8*psi/lambda-7));
        alpha1 = alpha+1;
        
     
     
        error=1;
        
        lam1s=lam1/s;
        
        x = b{j};
        u{1}=zeros(m,n);
        u{2}=zeros(m,n);
        x_o=x;
        u_o=u;
        p = alpha1*x-alpha*x_o;
        q={};
        q{1} = alpha1*u{1}-alpha*u_o{1};
        q{2} = alpha1*u{2}-alpha*u_o{2};
        w1_o = zeros(m,n);
        w21_o = zeros(m,n);
        w22_o = zeros(m,n);
        
    
        r=0;
        
        errorFPDHFRI = 1;
        tic
        while errorFPDHFRI > tol && r < niter;
                r=r+1;
                x_o=x;
                u_o=u;
                
                Cp = C(p);
        
                z_= p - t*( LT(q) + Cp + D(p) );
                z = min(max(0,z_),1);
            
                x = z-t*(C(z)-Cp);
        
                uu = L( z + x -p );
                
                u1 = (q{1}+s*uu{1})./s;
                prox = sign(u1) .* max(0, abs(u1) - lam1s);
                u{1} = (u1-prox)*s;
                
        
                u2 = (q{2}+s*uu{2})./s;
                prox = sign(u2) .* max(0, abs(u2) - lam1s);
                u{2} = (u2-prox)*s;
                
                
          
                w1 = lambda*x+lambda1*p;
                w21 = lambda*u{1}+lambda1*q{1};
                w22 = lambda*u{2}+lambda1*q{2};
                
               errorFPDHFRI = sqrt((norm(w1-w1_o,'fro')^2+norm(w21-w21_o,'fro')^2+norm(w22-w22_o,'fro')^2) /(norm(w1,'fro')^2 + norm(w21,'fro')^2+norm(w22,'fro')^2));
        
                
                p = alpha1*w1-alpha*w1_o;
                q{1} = alpha1*w21-alpha*w21_o;
                q{2} = alpha1*w22-alpha*w22_o;
                
                w1_o = w1;
                w21_o = w21;
                w22_o = w22;
        end
        time_FPDHFRI(j,i) = toc;
        iter_FPDHFRI(j,i) = r;
        x_FPDHFI{j,i}=x;
        xx = L(x);
        fo_FPDHFRI(j,i) = lam3*norm(K(x)-b{j},2)^2/2+lam1*sum(sum(abs(xx{1})))+lam1*sum(sum(abs(xx{2})))+huberf(wav(x),lam2,delta);
        SNR_FPDHFRI(j,i) = - 20*log10(norm(zz - x,2)/norm(zz,2));
    
    end
end


