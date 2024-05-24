function QuadrotorModel(Kp,Ki,Kd,tsamp,tfin,Mode)

% Plant
g = 9.81;
m = 1.5;
Kt = 120;
Cd = 0.1;
f = @(X,U) [X(2)
            (Kt*U - Cd*X(2))/m - g];

% Simulation properties
dt = 0.001;
t = 0:dt:tfin;
X = [0 0]';
r = 1;
eint = 0;
epr = r;

% Simulation loop
for i = 1:length(t)
    
    % Output
    Y(:,i) = [0 0 X(1) 0 0 X(2)]';
    
    % Controller
    switch Mode
        case 1
            e = r - X(1);
        case 2
            e = r - X(2);
    end
    switch Mode
        case 0
            U = 1;
        case {1,2}
            eint = eint + e*dt;
            edot = (e - epr)/dt;
            epr = e;
            U = Kp*e + Ki*eint + Kd*edot;
    end
    
    % Plant
    Xdot = f(X,U);
    
    % Integrate
    k1 = f(X,U);
    k2 = f(X+dt*k1/2,U);
    k3 = f(X+dt*k2/2,U);
    k4 = f(X+dt*k3,U);
    X = X + dt*(k1 + 2*k2 + 2*k3 + k4)/6;
    
end

% Run animation
fg = gcf;
figure(fg)
% set(fg,'OuterPosition',CentreFig(1300,800));
QuadrotorVis(t,Y,tsamp,fg,r,Mode);