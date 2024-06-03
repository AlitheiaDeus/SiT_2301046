%%
% Author: Abdul Karim Bin Abdul Mathanni
% SIT-UofG Mechanical Engineering Student
% Project: Visualisation of Quanser 2DoF Helicopter
% University of Glasgow Singapore, June 2022
% Academic Supervisor: Dr Arturo Molina-Cristobal
%% Example script to visualize the aircraft simulation data
% Add the path of the aircraft_3d_animation function
%addpath('../src/');
% path of the *.mat file containing the 3d model information
model_info_file = 'QUANSER_BODY16.mat';

% Load the simulation data
% load('scissors_maneuver.mat')
% load('breakaway_maneuver.mat')
% load('split_s_maneuver.mat')
% define the reproduction speed factor
speedx = 1; 
% Do you want to save the animation in a mp4 file? (0.No, 1.Yes)
isave_movie = 0;
% Movie file name
movie_file_name = '';

% -------------------------------------------------------------------------
% The frame sample time shall be higher than 0.02 seconds to be able to 
% update the figure (CPU/GPU constraints)\

signals.time(:,1)=t';
signals.signals.values(:,2)=y_pitch';
signals.signals.values(:,4)=y_yaw';
%frame_sample_time = max(0.04, t(2)-t(1));

 times = signals.time(:,1);
 pit = signals.signals.values(:,2);
 yaw = signals.signals.values(:,4);
 frame_sample_time = max(0.04, times(2)-times(1));


% Resample the time vector to modify the reproduction speed
t_new   = times(1):frame_sample_time*(speedx):times(end);
% Resample the recorded data
pit= interp1(times, pit, t_new','linear');
% We have to be careful with angles with ranges
%y_new = atan2(interp1(times, sin(f), t_new','linear'), interp1(times, cos(signals.signals.values(:, 4)), t_new','linear')) * pi / 180;
y_new = interp1(times,yaw, t_new', 'linear');
% Assign the data
heading_deg           =  y_new;
pitch_deg             =  pit;
bank_deg              =  y_new;
angle_of_attack_deg   =  pit;
mach                  =  y_new;
altitude_ft           = -y_new;
nz_g                  =  y_new;
% Flight control surfaces
dr     = pit;
% Control array assignation
% (modify the order according to your particular 3D model)
controls_deflection_deg = [dr(:),dr(:), dr(:)];

%% Run aircraft_3d_animation function
% -------------------------------------------------------------------------
model_3d_animation(model_info_file,...
    heading_deg, ...            Heading angle [deg]
    pitch_deg, ...              Pitch angle [deg]
    bank_deg, ...               Roll angle [deg]
    controls_deflection_deg, ...Flight control deflection (each column is a control surface)
    frame_sample_time, ...      Sample time [sec]
    speedx, ...                 Reproduction speed
    t_new, ...
    isave_movie, ...            Save the movie? 0-1
    movie_file_name);           % Movie file name
disp('--------------------------------')
disp('--------------------------------')
disp('The Visualisation tool called the Simulink Model')
disp('Visualisation by Student: A. K. Bin Abdul Mathanni,2022')
disp('SIT-UofG Mechanical Engineering Student')
disp('Project: Visualisation of Quanser 2DoF Helicopter')
disp('University of Glasgow Singapore, June 2022')
disp('Academic Supervisor: Dr Arturo Molina-Cristobal')
disp('------End of Visualisation-----')
disp('--------------------------------')
% 
