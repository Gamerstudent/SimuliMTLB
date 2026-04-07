clear; clc;

%% Geometry
L = 1.0;
W = 0.5;
H = 0.5;

Nx = 30; Ny = 30; Nz = 30;

dx = L/(Nx-1);
dy = W/(Ny-1);
dz = H/(Nz-1);

%% Material properties
k = 0.8;        % thermal conductivity
rho = 1.2;      % density (air approx)
cp = 1000;      % heat capacity

alpha = k/(rho*cp);   % thermal diffusivity

%% Time stepping
dt = 0.05;      % time step
t_end = 200;    % simulation time

%% Initial temperature
T = ones(Nx,Ny,Nz)*300;

%% Heater power region
Q = zeros(Nx,Ny,Nz);
Q(3,:,:) = 5e5;         % heaters on one side
Q(end-2,:,:) = 5e5;

%% Ambient boundary
T_ambient = 300;
h = 10;  % convective heat loss coefficient

%% Time loop
for t = 0:dt:t_end

    T_old = T;

    for i = 2:Nx-1
        for j = 2:Ny-1
            for k_idx = 2:Nz-1

                laplacian = ...
                    (T_old(i+1,j,k_idx) + T_old(i-1,j,k_idx) + ...
                     T_old(i,j+1,k_idx) + T_old(i,j-1,k_idx) + ...
                     T_old(i,j,k_idx+1) + T_old(i,j,k_idx-1) ...
                     - 6*T_old(i,j,k_idx))/dx^2;

                T(i,j,k_idx) = T_old(i,j,k_idx) + ...
                    alpha*dt*laplacian + ...
                    dt*Q(i,j,k_idx)/(rho*cp);

            end
        end
    end

    %% Convective cooling at boundaries

    T(1,:,:) = T(2,:,:) - dx*h/k*(T(2,:,:)-T_ambient);
    T(end,:,:) = T(end-1,:,:) - dx*h/k*(T(end-1,:,:)-T_ambient);

    T(:,1,:) = T(:,2,:) - dy*h/k*(T(:,2,:)-T_ambient);
    T(:,end,:) = T(:,end-1,:) - dy*h/k*(T(:,end-1,:)-T_ambient);

    T(:,:,1) = T(:,:,2) - dz*h/k*(T(:,:,2)-T_ambient);
    T(:,:,end) = T(:,:,end-1) - dz*h/k*(T(:,:,end-1)-T_ambient);

    %% Visualization every few steps

    if mod(round(t/dt),20)==0

        [X,Y,Z] = meshgrid(linspace(0,L,Nx),linspace(0,W,Ny),linspace(0,H,Nz));

        slice(X,Y,Z,T,[L/2],[W/2],[H/2])
        %colormap hot
        colorbar
        caxis([300 1200])

        title(['Furnace Temperature at t = ',num2str(t),' s'])
        xlabel('X')
        ylabel('Y')
        zlabel('Z')

        drawnow

    end

end