%% 2DOF Helicopter Parameters
% distance from pivot to center of rotor/propeller (m)
Dt = 0.16743; 
% total aero body mass (kg)
Mb = 1.07; 
% distance from aero body plane to center of mass below (m)
Dm = 2.4e-3; 
% Moment of inertia about pitch-axis (kg-m^2)
Jp = 23188500.45 / 1000 * 0.001^2; % from CAD: 23188500.45 g-mm^2
% Moment of inertia about pitch-axis (kg-m^2)
Jy = 23810415.95 / 1000 * 0.001^2; % from CAD: 23810415.95 g-mm^2;
% Gravity (m/s^2)
g = 9.81;

%%
% These parameters were found experimentally.
% 
%% Pitch
% Stiffness (N-m/V)
Ksp = 0.00744;
% Damping (N-m/V)
Dp = 0.00199;
% Pitch thrust gain (N/V)
Kpp = 0.00321;
% Pitch from yaw thrust gain (N/V)
Kpy = 0.00137;
% 
%% Yaw
% Damping (N-m/V)
Dy = 0.00192;
% Thrust (N-m/V)
Kyy = 0.00610;
% Yaw from pitch thrust gain (N/V)
Kyp = -0.00319;

