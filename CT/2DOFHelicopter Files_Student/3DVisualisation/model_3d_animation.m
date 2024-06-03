function [] = model_3d_animation(...         
    mat_file_name, ...
    heading_deg, ...            
    pitch_deg, ...              
    bank_deg, ...                                        
    controls_deflection_deg, ...
    frame_time, ...             
    speedx, ...
    t_new, ...
    isave_movie, ...           
    movie_file_name)           
%% Function Name: aircraft_3d_animation
%
% Description: A visualization tool that animates a 3D model using 
% simulation data
% % File modified by Abdul Karim Bin Abdul Mathanni
% SIT-UofG Mechanical Engineering Student
% Project: Visualisation of Quanser 2DoF Helicopter
% University of Glasgow Singapore, June 2022
% Academic Supervisor: Dr Arturo Molina-Cristobal
%
% Assumptions: None
%
% Inputs:
%     model_mat_file			Mat file containing the Model3D structure. 
%                               (see "./import_stl_model/import_stl_model.m")
%                               This structure must have the following fields:
%                                 - Aircraft: a structure vector containing:
%                                      * stl_data.vertices (vertices information for the
%                                                           patch command)
%                                      * stl_data.faces  (faces information for the patch command)
%                                      * model: stl file (string)
%                                      * color: surface color
%                                      * alpha: transparency
%                                 - Controls: a structure vector containing:
%                                      * stl_data.vertices (vertices information for the
%                                                           patch command)
%                                      * stl_data.faces  (faces information for the patch command)
%                                      * model: stl file (string)
%                                      * label: string
%                                      * color: surface color
%                                      * rot_offset_deg: rotation offset (deg)
%                                      * rot_point: rotation point coordinates 
%                                      * rot_vect: rotation vector
%                                      * max_deflection: maximum deflection angles [min, max]
%     heading_deg:              Heading angle [deg]
%     pitch_deg:                Pitch angle [deg]
%     bank_deg:                 Roll angle [deg]
%     roll_command:             Roll  stick command [-1,+1] [-1 -> left,            +1 -> right]
%     pitch_command:            Pitch stick command [-1,+1] [-1 -> full-back stick, +1 -> full-fwd stick]
%     angle_of_attack_deg:      AoA [deg]
%     angle_of_sideslip_deg:    AoS [deg]
%     fligh_path_angle_deg:     Flight path angle [deg]
%     mach:                     Mach number
%     altitude_ft:              Altitude [ft]
%     nz_g:                     Vertical load factor [g]
%     controls_deflection_deg:  Flight control deflection (deg) (each column is a different control surface)
%     frame_time:               Sample time [sec]
%     speedx:                   Reproduction speed factor
%     isave_movie:              Save the animation in a movie file? 0-1
%     movie_file_name:          Movie file name
%
% Outputs:
%     none
%
% $Revision: R2018b$ 
% $Author: Rodney Rodriguez Robles$
% $Date: January 25, 2021$
%------------------------------------------------------------------------------------------------------------

% Check Matlab version
% This is a "just in case" protection, I have not checked the compatibility for
% earlier versions than 2018b
% MatVersion = version('-release');
% MatVersion = str2double(MatVersion(1:end-1));
% if MatVersion < 2018
%     error('MATLAB version not supported [ < 2018 ]');
% end
% Select and load 3D model (see "generate_mat_from_stl.m" for examples)
load(mat_file_name, 'Model3D');
% Open the video output if we are recording the movie
if isave_movie == 1
    aviobj = VideoWriter(movie_file_name, 'MPEG-4');
    aviobj.Quality = 100;  % movie quality
    aviobj.FrameRate = 1/frame_time;
    open(aviobj);
end
% Get maximum dimension including all the aircraft's parts
AC_DIMENSION = max(max(sqrt(sum(Model3D.Aircraft(1).stl_data.vertices.^2, 2))));
for i=1:length(Model3D.Control)
    AC_DIMENSION = max(AC_DIMENSION, max(max(sqrt(sum(Model3D.Control(i).stl_data.vertices.^2, 2)))));
end

%% Initialize the figure
hf = figure;
AX = axes('position',[0.0 0.0 1 1]);
axis off
scrsz = get(0, 'ScreenSize');
set(gcf, 'Position',[scrsz(3)/40 scrsz(4)/12 scrsz(3)/2*1.0 scrsz(3)/2.2*1.0], 'Visible', 'on');
set(AX, 'color', 'none');
axis('equal')
hold on;
cameratoolbar('Show')

% Initializate transformation group handles
% -------------------------------------------------------------------------
% Aircraft transformation group handle
AV_hg         = hgtransform;
% controls_deflection_deg transformation group handles
% zeros is a Matrix
CONT_hg       = zeros(1,length(Model3D.Control));
for i=1:length(Model3D.Control)
    CONT_hg(i) = hgtransform('Parent', AV_hg, 'tag', Model3D.Control(i).label);
end
% Circles around the aircraft transformation group handles
% euler_hgt(1)  = hgtransform('Parent',           AX, 'tag', 'OriginAxes');
% euler_hgt(2)  = hgtransform('Parent', euler_hgt(1), 'tag', 'roll_disc');
% euler_hgt(4)  = hgtransform('Parent', euler_hgt(1), 'tag', 'heading_disc');
% euler_hgt(5)  = hgtransform('Parent', euler_hgt(2), 'tag', 'roll_line');
% euler_hgt(7)  = hgtransform('Parent', euler_hgt(4), 'tag', 'heading_line');
% Plot objects
% -------------------------------------------------------------------------
% Plot airframe
AV = zeros(1, length(Model3D.Aircraft));
for i = 1:length(Model3D.Aircraft)
    AV(i) = patch(Model3D.Aircraft(i).stl_data,  ...
        'FaceColor',        Model3D.Aircraft(i).color, ...
        'EdgeColor',        'none',        ...
        'FaceLighting',     'gouraud',     ...
        'AmbientStrength',   0.15,          ...
        'LineSmoothing',    'on',...
        'Parent',            AV_hg, ...
        'LineSmoothing', 'on');
end
CONT = zeros(1, (length(Model3D.Control)));
% Plot controls_deflection_deg
for i=1:length(Model3D.Control)
    CONT(i) = patch(Model3D.Control(i).stl_data,  ...
        'FaceColor',        Model3D.Control(i).color, ...
        'EdgeColor',        'none',        ...
        'FaceLighting',     'gouraud',     ...
        'AmbientStrength',  0.15,          ...
        'LineSmoothing', 'on',...
        'Parent',           CONT_hg(i));
end
% Fixing the axes scaling and setting a nice view angle
axis('equal');
axis([-1, 1, -1, 1, -1, 1] * 2.0 * AC_DIMENSION)
set(gcf, 'Color', [1, 1, 1])
axis off
view([100, 20])
zoom(2.0);
% Add a camera light, and tone down the specular highlighting
camlight('headlight');
material('dull');

%% Plot Euler angles references 
% R = AC_DIMENSION;

% % Plot outer circles
% phi = (-pi:pi/36:pi)';
% D1 = [sin(phi) cos(phi) zeros(size(phi))];
% plot3(R * D1(:,1), R * D1(:,2), R * D1(:,3), 'Color', 'b', 'tag', 'Zplane', 'Parent', euler_hgt(4));
% plot3(R * D1(:,3), R * D1(:,1), R * D1(:,2), 'Color', 'r', 'tag', 'Xplane', 'Parent', euler_hgt(2));
% 
% % Plot +0,+90,+180,+270 Marks
% S = 0.95;
% phi = -pi+pi/2:pi/2:pi;
% D1 = [sin(phi); cos(phi); zeros(size(phi))];
% plot3([S * R * D1(1, :); R * D1(1, :)],[S * R * D1(2, :); R * D1(2, :)],[S * R * D1(3, :); R * D1(3, :)], 'Color', 'b', 'tag', 'Zplane', 'Parent',euler_hgt(4));
% plot3([S * R * D1(3, :); R * D1(3, :)],[S * R * D1(1, :); R * D1(1, :)],[S * R * D1(2, :); R * D1(2, :)], 'Color', 'r', 'tag', 'Xplane', 'Parent',euler_hgt(2));
% text(R * 1.05 * D1(1, :), R * 1.05 * D1(2, :), R * 1.05 * D1(3, :), {'N', 'E', 'S', 'W'}, 'Fontsize',9, 'color', [0 0 0], 'HorizontalAlign', 'center', 'VerticalAlign', 'middle');
% 
% % Plot +45,+135,+180,+225,+315 Marks
% S = 0.95;
% phi = -pi+pi/4:2*pi/4:pi;
% D1 = [sin(phi); cos(phi); zeros(size(phi))];
% plot3([S*R * D1(1, :); R * D1(1, :)],[S*R * D1(2, :); R * D1(2, :)],[S*R * D1(3, :); R * D1(3, :)], 'Color', 'b', 'tag', 'Zplane', 'Parent',euler_hgt(4));
% text(R * 1.05 * D1(1, :), R * 1.05 * D1(2, :), R * 1.05 * D1(3, :), {'NW', 'NE', 'SE', 'SW'}, 'Fontsize',8, 'color',[0 0 0], 'HorizontalAlign', 'center', 'VerticalAlign', 'middle');
% 
% % 10 deg sub-division marks
% S = 0.98;
% phi = -180:10:180;
% phi = phi*pi / 180;
% D1 = [sin(phi); cos(phi); zeros(size(phi))];
% plot3([S * R * D1(1, :); R * D1(1, :)],[S * R * D1(2, :); R * D1(2, :)],[S * R * D1(3, :); R * D1(3, :)], 'Color', 'b', 'tag', 'Zplane', 'Parent', euler_hgt(4));
% plot3([S * R * D1(3, :); R * D1(3, :)],[S * R * D1(1, :); R * D1(1, :)],[S * R * D1(2, :); R * D1(2, :)], 'Color', 'r', 'tag', 'Xplane', 'Parent', euler_hgt(2));
% 
% % Guide lines
% plot3([-R, R], [ 0, 0], [0, 0], 'b-', 'tag', 'heading_line', 'parent', euler_hgt(7));
% plot3([ 0, 0], [-R, R], [0, 0], 'r-', 'tag',    'roll_line', 'parent', euler_hgt(5));


% Initialize text handles
FontSize    = 13;
text_color  = [1, 0, 1];
text_color1  = [0, 0, 1];
%text_color2  = [0.3010 0.7450 0.9330];
font_name   = 'Consolas';
hdle_text_t                 = text(0.45 * AC_DIMENSION * 1.5, 0.55 * AC_DIMENSION * 1.5, 't =  0 sec', 'Color',text_color, 'FontSize',FontSize, 'FontName', font_name);
hdle_text_th                = text(0.45 * AC_DIMENSION * 1.5, 0.55 * AC_DIMENSION * 1.5, 0.20 * AC_DIMENSION * 1.5, '', 'Color',text_color, 'FontSize', FontSize, 'FontName', font_name);
hdle_text_phi               = text(0.45 * AC_DIMENSION * 1.5, 0.55 * AC_DIMENSION * 1.5, 0.43 * AC_DIMENSION * 1.5, '', 'Color',text_color, 'FontSize', FontSize, 'FontName', font_name);
hdle_text_UGS  = text(2.0 * AC_DIMENSION * 1.5, 0.03 * AC_DIMENSION * 1.5, 1.10 * AC_DIMENSION * 2, '', 'Color',text_color1, 'FontSize', FontSize, 'FontName', font_name);
%hdle_text_UGS2 = text(2.0 * AC_DIMENSION * 1.5, 0.1 * AC_DIMENSION * 1.5, 1.05 * AC_DIMENSION * 2, '', 'Color',text_color2, 'FontSize', FontSize, 'FontName', font_name);
% Plot axes on top of graph
x = t_new;
y = pitch_deg;
z = heading_deg;
h_stick     = axes('Position',[0.8, 0.05, 0.15, 0.15], 'FontSize', 10);
hold(h_stick, 'on')
hold on; box on;
xlabel('Time (s)'); ylabel('Pitch (deg)')

h_stick2    = axes('Position',[0.1, 0.05, 0.15, 0.15], 'FontSize', 10);
hold(h_stick2, 'on')
hold on; box on;
xlabel('Time (s)'); ylabel('Yaw (deg)')
%text(0.5,0.5,'1DoF Helicopter')

%% Animation Loop
% Maximum and minimum surfaces' deflection
% max_deflection = reshape([Model3D.Control(:).max_deflection], 2, length(Model3D.Control(:)));
% Refresh plot for flight visualization
tic;
for i=1:length(heading_deg)

    % AIRCRAFT BODY
    M1 = makehgtform('zrotate',  -heading_deg(i) * pi / 180);  % Heading rotation
    %M2 = makehgtform('yrotate', pitch_deg(i) * pi / 180);  % Pitch rotation
    %M3 = makehgtform('xrotate', -bank_deg(i) * pi / 180);  % bank_deg rotation
    set(AV_hg, 'Matrix',M1)
    
    % controls_deflection_deg
    for j=1:length(Model3D.Control)
        M1 = makehgtform('translate', -Model3D.Control(j).rot_point);   % Heading
        M2 = makehgtform('axisrotate', Model3D.Control(j).rot_vect, controls_deflection_deg(i, j) * pi / 180);  % Pitch
        M3 = makehgtform('translate', Model3D.Control(j).rot_point);  % bank_deg
        set(CONT_hg(j), 'Matrix', M3 * M2 * M1);
    end
   
    % Compute Aerodynamic Speed Vector
    set(hdle_text_t, 'String',sprintf('t= %3.2f sec',(i-1)*frame_time*speedx))
    set(hdle_text_th, 'String',strcat('\theta=',num2str(pitch_deg(i), '%2.1f'), ' deg'))
    set(hdle_text_phi, 'String',strcat('\phi=',num2str(bank_deg(i), '%2.1f'), ' deg'))
    set(hdle_text_UGS  , 'String','2DoF Helicopter, SIT-UGS Aerospace Control')
    %set(hdle_text_UGS2  , 'String','(2022)')
    % Plot Graph (Pitch)
%     plot(h_stick, x(i),y(i), 'Color', [.6 0 0])
%     hold on
    plot(h_stick, x(1:i), y(1:i), 'Color', [0 0 0])
    pause(0.01)

    % Plot Graph (Heading)
%     plot(h_stick2, x(i),z(i), 'Color', [.6 0 0])
%     hold on
    plot(h_stick2, x(1:i), z(1:i), 'Color', [0 0 0])
    pause(0.01)


    % Detect control surfaces saturations
%     idx_sat = controls_deflection_deg(i, :) >= max_deflection(2, :)*0.99 | controls_deflection_deg(i, :) <= max_deflection(1, :)*0.99;
%     idx_nosat = ~idx_sat;
%     set(CONT(idx_sat), 'FaceColor', 'y');
%     set(CONT(idx_nosat), 'FaceColor', Model3D.Control(1).color);
%     
    % Real-time
    drawnow;
    if frame_time * i - toc > 0
        pause(max(0, frame_time * i - toc))
    end
    
    if isave_movie == 1
        writeVideo(aviobj, getframe(hf));
    end
%      for k = 1:length(x)
%    plot(x(1:k), y(1:k))
%     pause(0.005)
% %     if k ~= length(x)
%      %   clf
%     end
    
end
toc
if isave_movie == 1
    close(aviobj);
end
end
% end

% function Lbh = Lbh(heading_deg, pitch_deg, phi)
% % Rotation matrix from NED axis to Aircraft's body axis
% sps = sind(heading_deg);
% cps = cosd(heading_deg);
% sth = sind(pitch_deg);
% cth = cosd(pitch_deg);
% sph = sind(phi);
% cph = cosd(phi);
% Lb1 = [...
%     1   0   0
%     0   cph sph
%     0  -sph cph];
% L12 = [...
%     cth 0   -sth
%     0   1   0
%     sth 0   cth];
% L2h = [...
%     cps sps 0
%     -sps cps 0
%     0   0   1];
% Lbh = Lb1 * L12 * L2h;
% end

