% PID control of quadrotor height dynamics


% close all
clear
clc
clf

% Mode
% 0 - Open-loop response
% 1 - Height control
% 2 - Velocity control
Mode = 1;

% Gain selection
Kp = 1;
Ki = 0;
Kd = 0;

% Choose animation framerate and end time
f = 0.01;
T = 2;

% Send gains to PID controller and height model
QuadrotorModel(Kp,Ki,Kd,f,T,Mode)


