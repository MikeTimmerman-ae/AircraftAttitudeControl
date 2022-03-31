clear all
clc
close all

%% Data

% states = [angle of attack, pitch rate, pitch angle]

dt = 0.1;    %time step
T=10;        %simulation time
Nsim = T/dt; %number of simulation steps

%Continuous-time model dynamics

A = [ -0.5507  1      0;
      -9.7621 -0.9983 0;
       0       1      0 ];              
B = [ -0.0545
      -14.494
       0 ];               
C = [ 0 0 1 ];
D = [ 0 ];

%Discrete-time model dynamics

sys = ss(A, B, C, D);                       %state-space model
sys = c2d(sys, dt);                         %discrete state-space model

%LTI system definition

LTI.A = sys.A;

LTI.B = sys.B;

LTI.C = sys.C;

LTI.D = sys.D;
              
LTI.x0 = [ 0.0; 0; 0.0 ];
LTI.yref = [ 0.3491 ];

%Definition of system dimension

dim.nx = 3;     %state dimension
dim.nu = 1;     %input dimension
dim.ny = 1;     %output dimension
dim.N  = 5;     %horizon

%% Optimal Target Selection
 
A = [ eye(3)-sys.A, -sys.B;
      sys.C,   0 ];
 
B = [ 0;
      0;
      0;
      LTI.yref ];
 
X = linsolve(A,B);               %steady-state reference state and input

x_ref = X(1:dim.nx);             %reference state
u_ref = X(dim.nx+1);             %reference input
    
%Definition of quadratic cost function

weight.Q = diag([10, 10, 50]);                  %weight on state
weight.R = eye(dim.nu);                         %weight on input
weight.P = dare(sys.A,sys.B,weight.Q,weight.R); %terminal cost

%% Prediction model

predmodel = predictionmodel(LTI,dim);     %generation of state solution function
[H,h] = costgen(predmodel,weight,dim);    %writing cost function in quadratic form

%% Closed-loop simulation

%Receding horizon

x = zeros(dim.nx, Nsim); x(:,1) = LTI.x0;
u = zeros(dim.nu, Nsim); u(:,1) = 0;

% Set the constraint matrix

predmodel = predictionmodel(LTI, dim);
A = [ predmodel.S(dim.nx+1:end, :);
     -predmodel.S(dim.nx+1:end, :);
      eye(dim.N);
     -eye(dim.N)];

Xub = [35*pi/180; 100*pi/180; 30*pi/180];           %state upper bound constraint
Xlb = [-15*pi/180; -100*pi/180; -15*pi/180];        %state lower bound constraint
Uub = 25*pi/180;                                    %input upper bound constraint
Ulb = -25*pi/180;                                   %input lower bound constraint

%Closed-loop simulation loop

for k=1:Nsim 
    
    %Current state
    x_0 = x(:, k);
    
    %Set RHS of the constraints
    b = [ repmat(Xub, dim.N, 1) - predmodel.T(dim.nx+1:end, :)*x_0;
         -repmat(Xlb, dim.N, 1) + predmodel.T(dim.nx+1:end, :)*x_0;
          Uub*ones(1*dim.N,1);
         -Ulb*ones(1*dim.N,1) ];
    
    %Solve optimization problem   
    useq = sdpvar(dim.nu*dim.N,1);                                                 %define optimization variable
    Constraint=[A*useq<=b];                                                        %define constraints
    Objective = 0.5*useq'*H*useq+(h*[x_0; x_ref; u_ref])'*useq;                  %define cost function
    options = sdpsettings('verbose',0,'solver','quadprog','quadprog.maxiter',100);
    optimize(Constraint,Objective,options);                                         %solve optimization problem
    
    %Get optimal control input
    uopt = value(useq);
    u(:,k) = uopt(1:dim.nu);
    
    %Update system state
    x(:,k+1)=sys.A*x(:,k)+sys.B*u(:,k);

end

%% Plot data

subplot(3,1,1);
plot(dt*(1:Nsim), x(1, 1:Nsim)*180/pi, dt*(1:Nsim), x(3, 1:Nsim)*180/pi)
grid on
ylabel('Angle [deg]') 
legend({'Angle of Attack', 'Pitch Angle'},'Location','northeast')

subplot(3,1,2);
plot(dt*(1:Nsim), x(2, 1:Nsim)*180/pi)
grid on
ylabel('Anglular Velocity [deg/s]') 
legend({'Pitch Rate'},'Location','northeast')

subplot(3,1,3);
plot(dt*(1:Nsim), u(:, 1:Nsim)*180/pi)
grid on
xlabel('Time [s]') 
ylabel('Angle [deg]') 
legend({'Elevator Deflection'},'Location','northeast')
