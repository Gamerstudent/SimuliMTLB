clear; clc;

%% Geometry
L = 1.0;   
W = 0.5;   
H = 0.5;

Nx = 40; Ny = 40; Nz = 40;

dx = L/(Nx-1); 
dy = W/(Ny-1); 
dz = H/(Nz-1);

%% Ambient
T_ambient = 300;

%% Initialize temperature
T = ones(Nx,Ny,Nz)*T_ambient;

%% Conductivity field
k_field = ones(Nx,Ny,Nz)*0.8; % chamber gas

% material properties
k_steel = 16;
k_insulation = 0.2;

%% Define layered walls

% outer steel shell
k_field(1,:,:) = k_steel;
k_field(end,:,:) = k_steel;
k_field(:,1,:) = k_steel;
k_field(:,end,:) = k_steel;
k_field(:,:,1) = k_steel;
k_field(:,:,end) = k_steel;

% insulation layer
k_field(2,:,:) = k_insulation;
k_field(end-1,:,:) = k_insulation;
k_field(:,2,:) = k_insulation;
k_field(:,end-1,:) = k_insulation;
k_field(:,:,2) = k_insulation;
k_field(:,:,end-1) = k_insulation;

%% Heater layer (fixed temperature)

heater_temp = 1200;

T(3,:,:) = heater_temp;
T(end-2,:,:) = heater_temp;

%T(:,3,:) = heater_temp;
%T(:,end-2,:) = heater_temp;

%T(:,:,3) = heater_temp;
%T(:,:,end-2) = heater_temp;

%% Door opening (cold boundary)

T(:,2,:) = 300;

%% Solver parameters
max_iter = 5000;
tol = 1e-6;

for iter = 1:max_iter

    T_old = T;

    for i = 4:Nx-3
        for j = 4:Ny-3
            for k_idx = 4:Nz-3

                T(i,j,k_idx) = ( ...
                    T(i+1,j,k_idx) + T(i-1,j,k_idx) + ...
                    T(i,j+1,k_idx) + T(i,j-1,k_idx) + ...
                    T(i,j,k_idx+1) + T(i,j,k_idx-1) )/6;

            end
        end
    end

    % keep heater nodes fixed
    T(3,:,:) = heater_temp;
    T(end-2,:,:) = heater_temp;

    %T(:,3,:) = heater_temp;
    %T(:,end-2,:) = heater_temp;

    %T(:,:,3) = heater_temp;
    %T(:,:,end-2) = heater_temp;

    % convergence
    if max(abs(T(:)-T_old(:))) < tol
        disp(['Converged at iteration ',num2str(iter)])
        break
    end

end

%% Visualization

[X,Y,Z] = meshgrid(linspace(0,L,Nx), linspace(0,W,Ny), linspace(0,H,Nz));

slice(X,Y,Z,T,[L/2],[W/2],[H/2])

%colormap hot
colorbar

xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')

title('3D Furnace Temperature Distribution')