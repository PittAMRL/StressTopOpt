%-------------------------------------------------------------
%
%    Copyright (C) 2007 Krister Svanberg
%
%    This file, asymp.m, is part of GCMMA-MMA-code.
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
%------
%
%  Values on the parameters raa0, raa, low and upp are
%  calculated in the beginning of each outer iteration.
%
function [low,upp,raa0,raa] = ...
asymp(outeriter,n,xval,xold1,xold2,xmin,xmax,low,upp, ...
      raa0,raa,raa0eps,raaeps,df0dx,dfdx);
%
eeen=ones(n,1);
xmami = xmax - xmin;
xmamieps = 0.00001*eeen;
xmami = max(xmami,xmamieps);
raa0 = abs(df0dx)'*xmami;
raa0 = max(raa0eps,(0.1/n)*raa0);
raa  = abs(dfdx)*xmami;
raa  = max(raaeps,(0.1/n)*raa);
if outeriter < 2.5
  low = xval - 0.5*xmami;
  upp = xval + 0.5*xmami;
else
  xxx = (xval-xold1).*(xold1-xold2);
  factor = eeen;
  factor(find(xxx > 0)) = 1.2;
  factor(find(xxx < 0)) = 0.7;
  low = xval - factor.*(xold1 - low);
  upp = xval + factor.*(upp - xold1);
  lowmin = xval - 10*xmami;
  lowmax = xval - 0.01*xmami;
  uppmin = xval + 0.01*xmami;
  uppmax = xval + 10*xmami;
  low = max(low,lowmin);
  low = min(low,lowmax);
  upp = min(upp,uppmax);
  upp = max(upp,uppmin);
end
%---------------------------------------------------------------------
