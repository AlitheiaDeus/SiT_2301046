% Control Design Project
% Dr Arturo Molina-Cristobal
% Created May 2024
% Model of 2DOF Helicopter

clear
clc
FILE1_model_parameters

sim(['FILE3_s_2dof_Model_vis'], 30);
load('rsp_sim.mat');
rsp_data=rsp_sim;

t=rsp_data(1,:);
y_pitch=rsp_data(4,:); % Pitch position(deg)
y_pitch_r=rsp_data(3,:); % Pitch reference(deg)
u_p=rsp_data(2,:);      % Pitch Motor (V)

u_y=rsp_data(5,:)   ;   % Yaw Motor (V)
y_yaw_r=rsp_data(6,:) ; % Yaw reference(deg)
y_yaw  =rsp_data(7,:) ;% Yaw position(deg)

figure
subplot(2,2,1); plot(t,y_pitch,t,y_pitch_r)
title('Pitch position(deg)'); grid on
xlabel('time [sec]')
subplot(2,2,3); plot(t,u_p);
title('Pitch Motor (V)'); grid on
xlabel('time [sec]')
subplot(2,2,2); plot(t,y_yaw,t,y_yaw_r)
title('Yaw position(deg)'); grid on
xlabel('time [sec]')
subplot(2,2,4); plot(t,u_y)
title('Yaw Motor (V)'); grid on
xlabel('time [sec]')

%%
oldpath = path;
path(oldpath,'3DVisualisation\')
Visualise_model

