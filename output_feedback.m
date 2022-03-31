clear all
clc
close all

%% Data

% states = [angle of attack, pitch rate, pitch angle]

dt = 0.1;    %time step
T=15;        %simulation time
Nsim = T/dt; %number of simulation steps

%Continuous-time model dynamics

A = [ -0.5507  1      0;
      -9.7621 -0.9983 0;
       0       1      0 ];              
B = [ -0.0545
      -14.494
      0 ];               
C = [ 0 1 0;
      0 0 1 ];
D = [ 0;  0 ];

%Discrete-time model dynamics

sys = ss(A, B, C, D);                       %state-space model
sys = c2d(sys, dt);                         %discrete state-space model

%LTI system definition

LTI.A = sys.A;

LTI.B = sys.B;

LTI.C = sys.C;

LTI.D = sys.D;
              
LTI.x0 = [ 0.0; 0; 0.0 ];
LTI.d = [ 0.0; 0.0 ];
LTI.yref = [ 0.0; 0.0 ];

%Definition of system dimension

dim.nx = 3;     %state dimension
dim.nu = 1;     %input dimension
dim.ny = 2;     %output dimension
dim.nd = 2;     %disturbance dimension
dim.N =  5;     %horizon

%Definition of quadratic cost function

weight.Q = diag([10, 10, 500]);                  %weight on output
weight.R = eye(dim.nu);                          %weight on input
weight.P = dare(sys.A,sys.B,weight.Q,weight.R);  %terminal cost

%% Observer gain and disturbance matrices

K = place(LTI.A',(LTI.C*LTI.A)',[0.85; 0.90; 0.55])';      %tuning gain matrix
L = [K; eye(2)];                                           %observer gain
LTI.Bd = K;
LTI.Cd = eye(2) - C*K;

%% Check observability condition

rank = rank([eye(dim.nx)-LTI.A -LTI.Bd; LTI.C LTI.Cd]);
disp(rank)

%% Augmented system dynamics

LTIe.A = [ LTI.A LTI.Bd; zeros(dim.nd,dim.nx) eye(dim.nd) ];
LTIe.B = [ LTI.B; zeros(dim.nd,dim.nu) ];
LTIe.C = [ LTI.C LTI.Cd ];
LTIe.x0 = [ LTI.x0; LTI.d ];
LTIe.yref = LTI.yref;

%Definition of system dimension
dime.nx = 5;     %state dimension
dime.nu = 1;     %input dimension
dime.ny = 2;     %output dimension
dime.N  = dim.N; %horizon

%Definition of quadratic cost function
weighte.Q = blkdiag(weight.Q, zeros(dim.nd));            %weight on output
weighte.R = weight.R;                                    %weight on input
weighte.P = blkdiag(weight.P, zeros(dim.nd));            %terminal cost

%% Offset-free MPC using state-feedback

predmodel = predictionmodel(LTIe, dime);     %generation of state solution function
[He,he] = costgen(predmodel, weighte, dime); %writing cost function in quadratic form

%Receding horizon

xe = zeros(dime.nx, Nsim);                                  %extended state
y = zeros(dime.ny, Nsim);                                   %system output

xref = zeros(dime.nx, Nsim);                                %reference state
yref = zeros(dime.ny, Nsim);                                %reference outputs

u = zeros(dime.nu, Nsim);                                   %control input

xehat = zeros(dime.nx, Nsim);                               %state prediction


xe(:,1) = LTIe.x0;
y(:,1) = LTIe.C*LTIe.x0;
yref(:,1) = LTIe.yref;
xehat(:,1) = zeros(1,dime.nx);

% Set the constraint matrix

predmodel = predictionmodel(LTI, dim);
constraint.A = [ predmodel.S(dim.nx+1:end, :);
                 -predmodel.S(dim.nx+1:end, :);
                 eye(dime.N);
                -eye(dime.N)];

Xub = [35*pi/180; 100*pi/180; 30*pi/180];           %state upper bound constraint
Xlb = [-15*pi/180; -100*pi/180; -15*pi/180];        %state lower bound constraint
Uub = 25*pi/180;                                    %input upper bound constraint
Ulb = -25*pi/180;                                   %input lower bound constraint

%closed-loop simulation

for k=1:Nsim
    
    xe_0 = xe(:,k);                                   %current extended state
    dhat = xehat(end-dim.nd+1:end,k);                 %disturbance estimate

    %Set RHS of the constraints
    constraint.b = [ repmat(Xub, dime.N, 1) - predmodel.T(dim.nx+1:end, :)*xehat(1:3,k);
                    -repmat(Xlb, dime.N, 1) + predmodel.T(dim.nx+1:end, :)*xehat(1:3,k);
                     Uub*ones(1*dim.N,1);
                    -Ulb*ones(1*dim.N,1) ];
    
    %Compute optimal target selection
    eqconstraints = eqconstraintsgen(LTI,dim,dhat);
    [xr,ur] = optimalss(LTI,dim,weight,[],eqconstraints);
    xref(:,k) = [xr;dhat];

    useq = sdpvar(dime.nu*dime.N,1);                                         %define optimization variable
    Constraint = [constraint.A*useq<=constraint.b];                          %define constraints
    Objective = 0.5*useq'*He*useq+(he*[xehat(:,k); xref(:,k); ur])'*useq;    %define cost function
    options = sdpsettings('verbose',0,'solver','quadprog','quadprog.maxiter',100);
    optimize(Constraint,Objective, options);                                 %solve optimization problem
    
    %Get optimal control input
    uopt = value(useq);      
    u(:,k) = uopt(1:dim.nu);

    %Update state/output
    xe(:,k+1) = LTIe.A*xe_0 + LTIe.B*u(:,k);
    y(:,k+1) = LTIe.C*xe(:,k+1);
    clear u_uncon
    
    %Update extended state estimation
    xehat(:,k+1) = LTIe.A*xehat(:,k)+LTIe.B*u(:,k)+L*(y(:,k)-LTIe.C*xehat(:,k));
    
    %Vary distrubance and reference output signals
    if 6 > k*dt && k*dt > 3
        xe(4:5, k+1) = [-0.05;-0.05];
    end
    if 9 > k*dt && k*dt > 6
        xe(4:5, k+1) = [0.;0.];
    end    
    if 12 > k*dt && k*dt > 9
        xe(4:5, k+1) = [-0.05;-0.05];
        LTI.yref = [0.0; 0.3491];
        LTIe.yref = [0.0; 0.3491];
    end   
    if k*dt > 12
        xe(4:5, k+1) = [0.;0.];
    end
    
    %Save reference outputs
    yref(:, k+1) = LTIe.yref;

end

%% Plots

%State Trajectories

figure

subplot(3,1,1);
plot(dt*(0:Nsim-1), xe(1, 1:Nsim)*180/pi, dt*(0:Nsim-1), xe(3, 1:Nsim)*180/pi, dt*(0:Nsim-1), xref(3, 1:Nsim)*180/pi, dt*(0:Nsim-1), xe(4, 1:Nsim)*180/pi);
grid on
xticks([0 3 6 9 12 15]);
xlabel('Time [s]');
ylabel('Angle [deg]'); 
legend({'Angle of Attack', 'Pitch Angle', 'Reference Pitch', 'Disturbance'},'Location','northeast', 'FontSize', 6);

subplot(3,1,2);
plot(dt*(0:Nsim-1), xe(2, 1:Nsim)*180/pi)
grid on
xticks([0 3 6 9 12 15]);
xlabel('Time [s]');
ylabel('Anglular Velocity [deg/s]');
legend({'Pitch Rate'},'Location','northeast');

subplot(3,1,3);
plot(dt*(0:Nsim-1), u(:, 1:Nsim)*180/pi);
grid on
xticks([0 3 6 9 12 15]);
xlabel('Time [s]');
ylabel('Angle [deg]');
legend({'Elevator Deflection'},'Location','northeast');

% Observer error

figure
e = xe - xehat;

plot((0:T/dt)*dt,e),
grid on
xticks([0 3 6 9 12 15]);
xlabel('Time [s]') 
ylabel('Prediction error')
legend({'Angle of Attack', 'Pitch Rate', 'Pitch Angle', 'Disturbance 1', 'Disturbance 2'});

%Output evolution/tracking error

figure
e=(y(2,:) - yref(2,:))*180/pi;

subplot(2,1,1);
plot(dt*(0:Nsim), y(2,:)*180/pi, dt*(0:Nsim), yref(2,:)*180/pi)
grid on
xticks([0 3 6 9 12 15]);
xlabel('Time [s]') 
ylabel('Angel [deg]')

subplot(2,1,2)
plot((0:T/dt)*dt,e),
grid on
xticks([0 3 6 9 12 15]);
xlabel('Time [s]') 
ylabel('Tracking Error [deg]')

%Observer convergence

figure
e = xehat(4:5, 31:61) - xehat(4:5, 30:60);

loglog( (31:61), (1:31).^-1, (31:61), (abs(e(31:61))) );

grid on
xlabel('Time step k [-]') 
ylabel('Observer Error [-]')
legend({'-log(k)', 'log(ε)'});
