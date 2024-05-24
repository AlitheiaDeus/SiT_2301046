function QuadrotorVis(t,Y,dt,fg,r,Mode)

%Attitude = [0.1 0 0]'*ones(1,length(Time));

% Turn path history on/off
Path = 1;

% Theatre limits (m)
Hlim = 1;
Zmax = max(max(Y(1:3,:)));
Zmin = min(min(Y(1:3,:)));

% Load model
load('QuadrotorGeometry');

V0 = GeoModel.Vertices;
F = GeoModel.Faces;
C = GeoModel.Colours;
A = GeoModel.Alpha;

figure(fg)

a3d = subplot(2,3,[1 2 4 5]);
rotate3d;
hold on
axis equal
box on
h = patch('Vertices',V0,'Faces',F,'FaceColor','flat','FaceVertexCData',C,...
    'CDataMapping','direct','FaceVertexAlphaData',A,'AlphaDataMapping','none',...
    'FaceAlpha','flat');
xlabel('x')
ylabel('y')
zlabel('z')
xlim([-1 1]*Hlim)
ylim([-1 1]*Hlim)
zlim([Zmin Zmax+0.1])
grid on
view(3)
lighting phong
light('Position',[10 10 10])
colormap(gray)

if Path == 1;
    p1 = plot3(Y(1,1),Y(2,1),Y(3,1),'b','LineWidth',2);
end

% Height plot
zmax = max(Y(3,:));
zmin = min(Y(3,:));
ap = subplot(2,3,3);
pp = plot(t(1),Y(1,3),'LineWidth',2);
grid on, hold on, box on

% Velocity plot
vmax = max(Y(6,:));
vmin = min(Y(6,:));
av = subplot(2,3,6);
pv = plot(t(1),Y(1,6),'LineWidth',2);
grid on, hold on, box on

switch Mode
    case 1
        subplot(ap)
        plot([0 max(t)],[r r],'r--','LineWidth',2)
        zmax = max([zmax r]);        
    case 2
        subplot(av)
        plot([0 max(t)],[r r],'r--','LineWidth',2)
        vmax = max([vmax r]);
end

% Plot limits
zrng = zmax - zmin;
zlims = [zmin-0.1*zrng zmax+0.1*zrng];
vrng = vmax - vmin;
vlims = [vmin-0.1*vrng vmax+0.1*vrng];

% Labels
subplot(ap)
xlabel('Time (s)'), ylabel('Height (m)')
xlim([0 max(t)]), ylim(zlims)

subplot(av)
xlabel('Time (s)'), ylabel('Height velocity (m/s)')
xlim([0 max(t)]), ylim(vlims)

% Plot time
tp = t(1);

for i = 1:length(t)
    
    if tp-t(i) < 1e-5

        % Move model to current position
        x = Y(1,i);
        y = Y(2,i);
        z = Y(3,i);
        V(:,1) = V0(:,1) + x;
        V(:,2) = V0(:,2) + y;
        V(:,3) = V0(:,3) + z;
        
        subplot(a3d)
        set(h,'Vertices',V);
        set(fg,'Name',['Vehicle pose at t = ',num2str(t(i),'%.2f'),'s'])
        
        if Path == 1;
            set(p1,'xdata',Y(1,1:i),'ydata',Y(2,1:i),'zdata',Y(3,1:i));
        end
        
        set(pp,'xdata',t(1:i),'ydata',Y(3,1:i))
        set(pv,'xdata',t(1:i),'ydata',Y(6,1:i))

%         refreshdata
        drawnow
        
        tp = tp + dt;
    
    end

end