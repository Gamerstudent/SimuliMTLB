clear; clc;

%% Geometry and Grid
L = 1.0;   % Length (m)
W = 0.5;   % Width (m)
H = 0.5;   % Height (m)

Nx = 40; Ny = 40; Nz = 40;   % Grid resolution
dx = L/(Nx-1); dy = W/(Ny-1); dz = H/(Nz-1);

%% Material Properties
k = 0.8;   % Thermal conductivity (W/mK)
T_ambient = 300; % Ambient temperature (K)

%% Wall Layers (multi-layer insulation)
t1 = 0.05;   k1 = 1.5;   % refractory
t2 = 0.10;   k2 = 0.2;   % insulation
t3 = 0.005;  k3 = 16;    % steel shell

R_wall = t1/k1 + t2/k2 + t3/k3;
h_wall = 1 / R_wall;     % effective heat transfer coefficient

%% Boundary Conditions
T_left   = 1200;
T_right  = 1200;
T_top    = 1200;    
T_bottom = 1200;
T_front  = 300;    % Door side (cold spot) &  Insulated wall
T_back   = 1200;

%% Initialize Temperature Field
T = ones(Nx,Ny,Nz) * T_ambient;

% Apply boundary conditions
T(1,:,:)   = T_left;
T(end,:,:) = T_right;
T(:,1,:)   = T_front;
T(:,end,:) = T_back;
T(:,:,1)   = T_bottom;
T(:,:,end) = T_top;

%% Iterative Solver (Gauss-Seidel)
max_iter = 5000;
tol = 1e-6;

for iter = 1:max_iter
    T_old = T;
    
    for i = 2:Nx-1
        for j = 2:Ny-1
            for k_idx = 2:Nz-1
                T(i,j,k_idx) = ( ...
                    T(i+1,j,k_idx) + T(i-1,j,k_idx) + ...
                    T(i,j+1,k_idx) + T(i,j-1,k_idx) + ...
                    T(i,j,k_idx+1) + T(i,j,k_idx-1) ) / 6;
            end
        end
    end
    
    % Convergence check
    if max(abs(T(:)-T_old(:))) < tol
        disp(['Converged at iteration ', num2str(iter)]);
        break;
    end
end

%% Visualization
[X,Y,Z] = meshgrid(linspace(0,L,Nx), linspace(0,W,Ny), linspace(0,H,Nz));
slice(X,Y,Z,T, [L/2], [W/2], [H/2]); % mid-plane slices
colorbar;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('3D Furnace Temperature Distribution');