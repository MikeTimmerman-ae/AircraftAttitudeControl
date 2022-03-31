function[xr,ur]=optimalss(LTI,dim,weight,constraints,eqconstraints)

% Compute optimal steady-state solution
%   xr: state reference
%   ur: input reference 

H=blkdiag(zeros(dim.nx),eye(dim.nu));
h=zeros(dim.nx+dim.nu,1);

%set solver options
options = optimoptions(@quadprog); 
options.OptimalityTolerance=1e-20;
options.ConstraintTolerance=1.0000e-15;
options.Display='off';

%solve optimization problem
xur=quadprog(H,h,[],[],eqconstraints.A,eqconstraints.b,[],[],[],options);

%get solutions
xr=xur(1:dim.nx);
ur=xur(dim.nx+1:end);

end