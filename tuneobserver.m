clear all
clc
close all

%% Data

dt = 0.1;    % time step
T = 10;      % simulation horizon
Nsim = T/dt; % number of simulation steps

%Continuous-time model dynamics

A = [  -0.5507 1      0;
      -9.7621 -0.9983 0;
       0       1      0];              
B = [ -0.0545
      -14.494
      0];               
C = [ 0 1 0;
      0 0 1];
D = [ 0;  0];

%Discrete-time model dynamics

sys = ss(A, B, C, D);                       %state-space model
sys = c2d(sys, dt);                         %discrete state-space model

%LTI system definition

LTI.A = sys.A;

LTI.B = sys.B;

LTI.C = sys.C;

LTI.D = sys.D;

LTI.yref = [0; 0];

%% Observer gain and disturbance matrices

K = place(LTI.A',(LTI.C*LTI.A)',[0.85; 0.9; 0.55])';    %tuning gain matrix

L = [K; eye(2)];                                        %observer gain
LTI.Bd = K;
LTI.Cd = eye(2) - C*K;

%% Augmented system dynamics

LTIe.A = [LTI.A LTI.Bd; zeros(dim.nd,dim.nx) eye(dim.nd)];
LTIe.B = [LTI.B; zeros(dim.nd,dim.nu)];
LTIe.C = [LTI.C LTI.Cd];

eigs(LTIe.A-L*LTIe.C)                              %Observer error dynamics
