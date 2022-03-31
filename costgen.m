function [H,h]=costgen(predmodel,weighte,dime)

% Function computing the quadratic and the linear term of the optimization:
%   H: quadratic term
%   h: linear term

Qbar = blkdiag(kron(eye(dime.N),weighte.Q),weighte.P);
Rbar = kron(eye(dime.N),weighte.R);

% quadratic term
H = predmodel.S'*Qbar*predmodel.S+Rbar;   

% linear term
hx0 = predmodel.S'*Qbar*predmodel.T;
hxref = -predmodel.S'*Qbar*kron(ones(dime.N+1,1),eye(dime.nx));
huref = -Rbar*kron(ones(dime.N,1),eye(dime.nu));

h = [hx0 hxref huref];
 
end