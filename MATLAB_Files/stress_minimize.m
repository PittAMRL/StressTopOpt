function [f0val,df0dx,fval,dfdx]=stress_minimize(x,Hs,H)
nelx=200;nely=60;nelz=1;
pl=3;q=0.5;p=10;
x(:)=(H*x(:))./Hs;
[pnorm,pnorm_sen,MISES]=Stress_3D_Sensitivity_Comp(x,nelx,nely,nelz,pl,q,p);
figure(1)
x_plot=reshape(x,nely,nelx,nelz);
figure(1); contourf(flipud(x_plot(:,:,1)),[0.5,0.5]);
colormap([0,0,0]); set(gcf,'color','w'); axis equal; axis off;title('material layout');drawnow;
figure(2)
imagesc(reshape(MISES.*(0.5*sign(x-0.5)+0.5),nely,nelx));axis equal;axis off; colormap('jet');title('Von-Mises Stress');colorbar;drawnow;
dv = ones(nely,nelx,nelz)/(nelx*nely*nelz); 
sen(:) = H*(pnorm_sen(:)./Hs);
dv(:) = H*(dv(:)./Hs);
fval=[mean(x(:))-0.3];
dfdx=[dv(:)'];
df0dx=sen';
f0val=pnorm;