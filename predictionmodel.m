function predmodel=predictionmodel(LTIe,dime)

% Function computing the matrices of the state solution function:
%   x_k = Tx_0 + Su_k

%Prediction matrix from initial state: T=[A^k]
T = zeros(dime.nx*(dime.N+1),dime.nx);
T(1:dime.nx,:) = eye(dime.nx,dime.nx);
for k=0:dime.N-1
    T( (k+1)*dime.nx+1:(k+2)*dime.nx , : ) = T( k*dime.nx+1:(k+1)*dime.nx, : ) * LTIe.A;
end

%Prediction matrix from input: S=[B AB A^2B ... A^kB]
S=zeros(dime.nx*(dime.N+1),dime.nu*(dime.N));
for k=1:dime.N
    S( k*dime.nx+1:(k+1)*dime.nx , : ) = LTIe.A*S( (k-1)*dime.nx+1:k*dime.nx , : );
    S( k*dime.nx+1:(k+1)*dime.nx , (k-1)*dime.nu+1:k*dime.nu ) = LTIe.B;
end

predmodel.T = T;
predmodel.S = S;



