%-------------------------------------------------------------
%
%    Copyright (C) 2008 Krister Svanberg
%
%    This file, gcmmasub.m.m, is part of GCMMA-MMA-code.
%    
%    GCMMA-MMA-code is free software; you can redistribute it and/or
%    modify it under the terms of the GNU General Public License as 
%    published by the Free Software Foundation; either version 3 of 
%    the License, or (at your option) any later version.
%    
%    This code is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%    
%    You should have received a copy of the GNU General Public License
%    (file COPYING) along with this file.  If not, see 
%    <http://www.gnu.org/licenses/>.
%    
%    You should have received a file README along with this file,
%    containing contact information.  If not, see
%    <http://www.smoptit.se/> or e-mail mmainfo@smoptit.se or krille@math.kth.se.
%
%
%------
%
%  Version Feb 2008.
%
function [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,f0app,fapp] = ...
gcmmasub(m,n,iter,epsimin,xval,xmin,xmax,low,upp, ...
         raa0,raa,f0val,df0dx,fval,dfdx,a0,a,c,d);
%
eeen = ones(n,1);
zeron = zeros(n,1);
%
% Calculations of the bounds alfa and beta.
albefa = 0.01;
move =0.2;
%
zzz1 = low + albefa*(xval-low);
zzz2 = xval - move*(xmax-xmin);
zzz  = max(zzz1,zzz2);
alfa = max(zzz,xmin);
zzz1 = upp - albefa*(upp-xval);
zzz2 = xval + move*(xmax-xmin);
zzz  = min(zzz1,zzz2);
beta = min(zzz,xmax);
%
% Calculations of p0, q0, r0, P, Q, r and b.
xmami = xmax-xmin;
xmamieps = 0.00001*eeen;
xmami = max(xmami,xmamieps);
xmamiinv = eeen./xmami;
ux1 = upp-xval;
ux2 = ux1.*ux1;
xl1 = xval-low;
xl2 = xl1.*xl1;
uxinv = eeen./ux1;
xlinv = eeen./xl1;
%
p0 = zeron;
q0 = zeron;
p0 = max(df0dx,0);
q0 = max(-df0dx,0);
pq0 = p0 + q0;
p0 = p0 + 0.001*pq0;
q0 = q0 + 0.001*pq0;
p0 = p0 + raa0*xmamiinv;
q0 = q0 + raa0*xmamiinv;
p0 = p0.*ux2;
q0 = q0.*xl2;
r0 = f0val - p0'*uxinv - q0'*xlinv;
%
P = sparse(m,n);
Q = sparse(m,n);
P = max(dfdx,0);
Q = max(-dfdx,0);
PQ = P + Q;
P = P + 0.001*PQ;
Q = Q + 0.001*PQ;
P = P + raa*xmamiinv';
Q = Q + raa*xmamiinv';
P = P * spdiags(ux2,0,n,n);
Q = Q * spdiags(xl2,0,n,n);
r = fval - P*uxinv - Q*xlinv;
b = -r;
%
% Solving the subproblem by a primal-dual Newton method
[xmma,ymma,zmma,lam,xsi,eta,mu,zet,s] = ...
subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);
%
% Calculations of f0app and fapp.
ux1 = upp-xmma;
xl1 = xmma-low;
uxinv = eeen./ux1;
xlinv = eeen./xl1;
f0app = r0 + p0'*uxinv + q0'*xlinv;
fapp  =  r +   P*uxinv +   Q*xlinv;
%
%---------------------------------------------------------------------
