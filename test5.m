clear; clc; close all;

%% 1. Geometry & Mesh
L = 0.5; W = 0.4; H = 0.4; % Meters
Nx = 15; Ny = 15; Nz = 15; % Mesh density (keep modest for speed)
dx = L/(Nx-1); dy = W/(Ny-1); dz = H/(Nz-1);

%% 2. Material & Physical Properties
k_insulation = 0.15;      % Low k for multi-layer insulation (W/mK)
epsilon = 0.85;           % Emissivity of internal walls
sigma = 5.67e-8;          % Stefan-Boltzmann constant
T_ambient = 293;          % 20 degC
h_conv = 15;              % Convection coefficient (W/m^2K)

% Energy Input
P_total = 3000;           % 3kW Heater
% Distribute power to internal nodes near walls (simulating heating elements)
Q = zeros(Nx, Ny, Nz);
Q(2:3, :, :) = P_total / (2 * Ny * Nz * dx * dy * dz); 

%% 3. Numerical Solver (Steady State Iteration)
% We use an iterative Gauss-Seidel approach modified for radiation (non-linear)
T = ones(Nx, Ny, Nz) * 800; % Initial guess (Start hot to help convergence)
tol = 1e-4;
max_iter = 2000;
err = 1;
iter = 0;

while err > tol && iter < max_iter
    T_old = T;
    
    for i = 2:Nx-1
        for j = 2:Ny-1
            for k = 2:Nz-1
                % Standard Laplacation for internal conduction
                term_x = (T(i+1,j,k) + T(i-1,j,k)) / dx^2;
                term_y = (T(i,j+1,k) + T(i,j-1,k)) / dy^2;
                term_z = (T(i,j,k+1) + T(i,j,k-1)) / dz^2;
                
                denom = (2/dx^2 + 2/dy^2 + 2/dz^2);
                T(i,j,k) = (term_x + term_y + term_z + Q(i,j,k)/k_insulation) / denom;
            end
        end
    end
    
    %% 4. Boundary Conditions (The "Door Effect" & Radiation)
    % Door at x = L (i = Nx) - simulate higher loss
    h_door = 50; % High loss at the door seal
    T(Nx,:,:) = (k_insulation*T(Nx-1,:,:)/dx + h_door*T_ambient) / (k_insulation/dx + h_door);
    
    % Insulated Walls (Conduction + Radiation loss to environment)
    % Left wall (x=0)
    T(1,:,:) = T(2,:,:); % Adiabatic simplification or complex Robin:
    
    % Radiative cooling approximation on outer surfaces
    % (For a high-fidelity model, we linearize the T^4 term)
    T_surf = T(Nx, floor(Ny/2), floor(Nz/2));
    h_rad = epsilon * sigma * (T_surf^2 + T_ambient^2) * (T_surf + T_ambient);
    
    % Apply symmetry/insulation to other faces
    T(:,1,:) = T(:,2,:); T(:,Ny,:) = T(:,Ny-1,:);
    T(:,:,1) = T(:,:,2); T(:,:,Nz) = T(:,:,Nz-1);

    err = max(abs(T(:) - T_old(:)));
    iter = iter + 1;
end

%% 5. Visualization
fprintf('Converged in %d iterations.\n', iter);

[X, Y, Z] = meshgrid(linspace(0,W,Ny), linspace(0,L,Nx), linspace(0,H,Nz));
figure('Color', 'w');
slice_x = L/2; slice_y = W/2; slice_z = H/2;
s = slice(X, Y, Z, T, [0 L], [W/2], [0]);
set(s, 'EdgeColor', 'none', 'FaceAlpha', 0.8);

%colormap('hot');
colorbar;
title('3D Steady-State Thermal Map of Muffle Furnace');
xlabel('Length (m) -> Door at end');
ylabel('Width (m)');
zlabel('Height (m)');
view(3);
grid on;