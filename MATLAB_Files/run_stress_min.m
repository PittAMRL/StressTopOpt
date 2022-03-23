clc
clear
close('all');
warning('off','all')
nelx=200;
nely=60;
nelz=1;
x=0.3*ones(nely,nelx,nelz);
[Hs,H]=prepare_filter(2.5,nelx,nely,nelz);
m =1;
epsimin = 0.0000001;
n=length(x(:));
xval=x(:);
xold1   = xval;
xold2   = xval;
xlb = 1e-3*ones(n,1);
xub = 1*ones(n,1);
xmin    = xlb;
xmax    = xub;
low     = xlb;
upp     = xub;
c       = [1e4]';
d       = [0]';
a0      = 0;
a       = [0]';
raa0    = 0.0001;
raa     = 0.0001;
raa0eps = 0.0000001;
raaeps  = 0.0000001;
outeriter = 0;
maxoutit  = 120;
kkttol  = 0;
x_his=zeros(nelx*nely*nelz,maxoutit);
if outeriter < 0.5
[f0val,df0dx,fval,dfdx]=stress_minimize(xval,Hs,H);
innerit=0;
outvector1 = [outeriter innerit xval'];
outvector2 = [f0val fval'];
end
kktnorm = kkttol+1;
outit = 0;
while  outit < maxoutit
outit   = outit+1;
outeriter = outeriter+1;
%%%% The parameters low, upp, raa0 and raa are calculated:
[low,upp,raa0,raa] = ...
asymp(outeriter,n,xval,xold1,xold2,xmin,xmax,low,upp, ...
raa0,raa,raa0eps,raaeps,df0dx,dfdx);
[xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = ...
gcmmasub(m,n,outeriter,epsimin,xval,xmin,xmax,low,upp, ...
raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d);
xold2 = xold1;
xold1 = xval;
xval  = xmma;
[f0val,df0dx,fval,dfdx]=stress_minimize(xval,Hs,H);
% PRINT RESULTS
fprintf(' It.:%5i      P-norm Stress.:%11.4f   Vol.:%7.3f \n',outit,f0val, ...
    mean(xval(:)));
%%%% The residual vector of the KKT conditions is calculated:
[residu,kktnorm,residumax] = ...
kktcheck(m,n,xmma,ymma,zmma,lam,xsi,eta,mu,zet,s, ...
xmin,xmax,df0dx,fval,dfdx,a0,a,c,d);
outvector1 = [outeriter innerit xval'];
outvector2 = [f0val fval'];
x_his(:,outit)=xmma;
end