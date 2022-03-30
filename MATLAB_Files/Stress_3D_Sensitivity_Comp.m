function [pnorm,pnorm_sen,MISES]=Stress_3D_Sensitivity_Comp(x,nelx,nely,nelz,pl,q,p)
[KE,B,D]=brick_stiffnessMatrix();
%  MATERIAL PROPERTIES
E0 = 1;           % Young's modulus of solid material
Emin = 1e-9;      % Young's modulus of void-like material
% USER-DEFINED LOAD DOFs
[il,jl,kl] = meshgrid(nelx, 0, 0:nelz);                 % Coordinates
loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1-jl); % Load Node IDs
loaddof = 3*loadnid(:) - 1;                             % Load DOFs
% USER-DEFINED SUPPORT FIXED DOFs
[iif,jf,kf] = meshgrid(0,0:nely,0:nelz);                  % Coordinates
fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1-jf); % Fixed Node IDs
fixeddof = [3*fixednid(:); 3*fixednid(:)-1; 3*fixednid(:)-2]; % Fixed DOFs
nele = nelx*nely*nelz;
ndof = 3*(nelx+1)*(nely+1)*(nelz+1);
F = sparse(loaddof,1,-1,ndof,1);     %External force
U = zeros(ndof,1);
freedofs = setdiff(1:ndof,fixeddof);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1);
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);
sK = reshape(KE(:)*(Emin+x(:)'.^pl*(E0-Emin)),24*24*nele,1);
K = sparse(iK,jK,sK); K = (K+K')/2; %global stiffness matrix assembly
U(freedofs,:) = K(freedofs,freedofs)\F(freedofs,:);  
MISES=zeros(nele,1); %von Mises stress vector
S=zeros(nele,6); 
for i=1:nele
    temp=x(i)^q*(D*B*U(edofMat(i,:)))';
    S(i,:)=temp;
    MISES(i)=sqrt(0.5*((temp(1)-temp(2))^2+(temp(1)-temp(3))^2....
    +(temp(2)-temp(3))^2+6*sum(temp(4:6).^2)));
end
nele=size(x(:),1); %total element number
[ndof,~]=size(U);
DvmDs=zeros(nele,6);
dpn_dvms=(sum(MISES.^p))^(1/p-1);
index_matrix=edofMat';
pnorm=(sum(MISES.^p))^(1/p);
for i=1:nele
    DvmDs(i,1)=1/2/MISES(i)*(2*S(i,1)-S(i,2)-S(i,3));
    DvmDs(i,2)=1/2/MISES(i)*(2*S(i,2)-S(i,1)-S(i,3));
    DvmDs(i,3)=1/2/MISES(i)*(2*S(i,3)-S(i,1)-S(i,2));
    DvmDs(i,4)=3/MISES(i)*S(i,4);
    DvmDs(i,5)=3/MISES(i)*S(i,5);
    DvmDs(i,6)=3/MISES(i)*S(i,6);  
end
beta=zeros(nele,1);
for i=1:nele
    u=reshape(U(edofMat(i,:),:)',[],1);
    beta(i)=q*(x(i))^(q-1)*MISES(i)^(p-1)*DvmDs(i,:)*D*B*u;
end
T1=dpn_dvms*beta;
gama=zeros(ndof,1);
for i=1:nele
    index=index_matrix(:,i);
    gama(index)=gama(index)+x(i)^q*dpn_dvms*B'*D'*DvmDs(i,:)'*MISES(i).^(p-1);
end
lamda=zeros(ndof,1);
lamda(freedofs,:)=K(freedofs,freedofs)\gama(freedofs,:);
T2=zeros(nele,1);
for i=1:nele
    index=index_matrix(:,i);
    T2(i)=-lamda(index)'*pl*x(i)^(pl-1)*KE*U(index);
end
pnorm_sen=T1+T2;


function [KE,B,D]=brick_stiffnessMatrix()
% elastic matrix formulation
nu=0.3;
D = 1./((1+nu)*(1-2*nu))*[1-nu nu nu 0 0 0; nu 1-nu nu 0 0 0;...
    nu nu 1-nu 0 0 0; 0 0 0 (1-2*nu)/2 0 0; 0 0 0 0 (1-2*nu)/2 0;...
    0 0 0 0 0 (1-2*nu)/2];
%stiffness matrix formulation
A = [32 6 -8 6 -6 4 3 -6 -10 3 -3 -3 -4 -8;
    -48 0 0 -24 24 0 0 0 12 -12 0 12 12 12];
k = 1/144*A'*[1; nu];
K1 = [k(1) k(2) k(2) k(3) k(5) k(5);
    k(2) k(1) k(2) k(4) k(6) k(7);
    k(2) k(2) k(1) k(4) k(7) k(6);
    k(3) k(4) k(4) k(1) k(8) k(8);
    k(5) k(6) k(7) k(8) k(1) k(2);
    k(5) k(7) k(6) k(8) k(2) k(1)];
K2 = [k(9)  k(8)  k(12) k(6)  k(4)  k(7);
    k(8)  k(9)  k(12) k(5)  k(3)  k(5);
    k(10) k(10) k(13) k(7)  k(4)  k(6);
    k(6)  k(5)  k(11) k(9)  k(2)  k(10);
    k(4)  k(3)  k(5)  k(2)  k(9)  k(12)
    k(11) k(4)  k(6)  k(12) k(10) k(13)];
K3 = [k(6)  k(7)  k(4)  k(9)  k(12) k(8);
    k(7)  k(6)  k(4)  k(10) k(13) k(10);
    k(5)  k(5)  k(3)  k(8)  k(12) k(9);
    k(9)  k(10) k(2)  k(6)  k(11) k(5);
    k(12) k(13) k(10) k(11) k(6)  k(4);
    k(2)  k(12) k(9)  k(4)  k(5)  k(3)];
K4 = [k(14) k(11) k(11) k(13) k(10) k(10);
    k(11) k(14) k(11) k(12) k(9)  k(8);
    k(11) k(11) k(14) k(12) k(8)  k(9);
    k(13) k(12) k(12) k(14) k(7)  k(7);
    k(10) k(9)  k(8)  k(7)  k(14) k(11);
    k(10) k(8)  k(9)  k(7)  k(11) k(14)];
K5 = [k(1) k(2)  k(8)  k(3) k(5)  k(4);
    k(2) k(1)  k(8)  k(4) k(6)  k(11);
    k(8) k(8)  k(1)  k(5) k(11) k(6);
    k(3) k(4)  k(5)  k(1) k(8)  k(2);
    k(5) k(6)  k(11) k(8) k(1)  k(8);
    k(4) k(11) k(6)  k(2) k(8)  k(1)];
K6 = [k(14) k(11) k(7)  k(13) k(10) k(12);
    k(11) k(14) k(7)  k(12) k(9)  k(2);
    k(7)  k(7)  k(14) k(10) k(2)  k(9);
    k(13) k(12) k(10) k(14) k(7)  k(11);
    k(10) k(9)  k(2)  k(7)  k(14) k(7);
    k(12) k(2)  k(9)  k(11) k(7)  k(14)];
KE = 1/((nu+1)*(1-2*nu))*...
    [ K1  K2  K3  K4;
    K2'  K5  K6  K3';
    K3' K6  K5' K2';
    K4  K3  K2  K1'];
% strain matrix formulation
B_1=[-0.044658,0,0,0.044658,0,0,0.16667,0
0,-0.044658,0,0,-0.16667,0,0,0.16667
0,0,-0.044658,0,0,-0.16667,0,0
-0.044658,-0.044658,0,-0.16667,0.044658,0,0.16667,0.16667
0,-0.044658,-0.044658,0,-0.16667,-0.16667,0,-0.62201
-0.044658,0,-0.044658,-0.16667,0,0.044658,-0.62201,0];
B_2=[0,-0.16667,0,0,-0.16667,0,0,0.16667
0,0,0.044658,0,0,-0.16667,0,0
-0.62201,0,0,-0.16667,0,0,0.044658,0
0,0.044658,-0.16667,0,-0.16667,-0.16667,0,-0.62201
0.16667,0,-0.16667,0.044658,0,0.044658,-0.16667,0
0.16667,-0.16667,0,-0.16667,0.044658,0,-0.16667,0.16667];
B_3=[0,0,0.62201,0,0,-0.62201,0,0
-0.62201,0,0,0.62201,0,0,0.16667,0
0,0.16667,0,0,0.62201,0,0,0.16667
0.16667,0,0.62201,0.62201,0,0.16667,-0.62201,0
0.16667,-0.62201,0,0.62201,0.62201,0,0.16667,0.16667
0,0.16667,0.62201,0,0.62201,0.16667,0,-0.62201];
B=[B_1,B_2,B_3];
end
end
% =========================================================================
% The present code is part of the journal paper by Deng et al. 2021 and is
% derived from the code which was part of the paper by Liu et al. 2014.
% -------------------------------------------------------------------------
% Please send your suggestions and comments to: albertto@pitt.edu
% -------------------------------------------------------------------------
% The code is intended for educations purposes, and the details and
% extensions can be found in the paper:
% Deng, H., Vulimiri, P.S. & To, A.C. An efficient 146-line 3D sensitivity
% analysis code of stress-based topology optimization written in MATLAB.
% Optim Eng (2021). https://doi.org/10.1007/s11081-021-09675-3
% -------------------------------------------------------------------------
% Details of the finite element formulation from Liu et al. can be found in
% the paper:
% Liu, K., Tovar, A. An efficient 3D topology optimization code written in
% Matlab. Struct Multidisc Optim 50, 1175–1196 (2014).
% https://doi.org/10.1007/s00158-014-1107-x
% -------------------------------------------------------------------------
% The code can be downloaded from the website:
% https://github.com/PittAMRL/StressTopOpt
% -------------------------------------------------------------------------
% The code from Liu et al. from which this is derived from, can be
% downloaded from: http://www.top3dapp.com/
% -------------------------------------------------------------------------
% Disclaimer:
% The code may be distributed and used for educational purposes.
% The authors do not guarantee that the code is free from errors.
