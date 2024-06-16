% Constants
% Skillet Properties
k_skillet = 10; % W/(m*K)
rho_skillet = 3500; % kg/m^3
c_p_skillet = 1300; % J/(kg*K)
T_i_skillet = 24; % deg C
P_skillet = 500; % W
V_skillet = 0.09*0.21*0.005; % m^3

% Crust Properties
k_crust = 0.5; % W/(m*K)
rho_crust = 700; % kg/m^3
c_p_crust = 2500; % J/(kg*K)
T_i_crust = 4; % deg C
epsilon_m_crust = 0.05; % unitless
V_crust = pi/2*(0.04^2*0.2 - 0.034641^2*0.19); % m^3

% Filling Properties
k_filling = 1; % W/(m*K)
rho_filling = 1200; % kg/m^3
c_p_filling = 4200; % J/(kg*K)
T_i_filling = 4; % deg C
epsilon_m_filling = 0.4; % unitless
V_filling = pi*0.034641^2*0.19/2; % m^3

% Global Properties
h = 200; % W/(m^2*K)
T_ideal = 75; % deg C
T_tol = 4; % deg C
T_inf = 120; % deg C
t_i = 0; % s
t_step = 1; % s
t_f = 60; % s
d_step = 1e-3; % m
error_tol = 1e-3; % unitless
error = 1; % unitless
iter = 0; % unitless

% Guesses for Maximum Tastiness Score & Corresponding Microwave Power
guess_vec = [1000, 0; 2000, 0; 3000, 0]; % [W, unitless]

% Generate Geometry from STL
crust = importGeometry("Crust.stl");
filling = importGeometry("Filling.stl");
loaf = addCell(crust, filling);

% Generate FEA Model
model = createpde("thermal", "transient");
model.Geometry = loaf;
scale(model.Geometry, d_step);

% Plot of Geometry
pdegplot(model, "CellLabels", "on", "FaceAlpha", 0.5);

% Assign Regional Properties
thermalProperties(model, ThermalConductivity=k_skillet, MassDensity=rho_skillet, SpecificHeat=c_p_skillet, Cell=2);
thermalProperties(model, ThermalConductivity=k_crust, MassDensity=rho_crust, SpecificHeat=c_p_crust, Cell=1);
thermalProperties(model, ThermalConductivity=k_filling, MassDensity=rho_filling, SpecificHeat=c_p_filling, Cell=3);

% Set Initial Conditions
thermalIC(model, T_i_skillet, "Cell", 2);
thermalIC(model, T_i_crust, "Cell", 1);
thermalIC(model, T_i_filling, "Cell", 3);

% Set Boundary Conditions
thermalBC(model, "Face", [15, 16, 18, 19, 20], "HeatFlux", 0); % Insulated Skillet Edges
thermalBC(model, "Face", [2, 12, 13, 9, 8, 7, 6, 14, 5, 4, 11, 1, 10], "ConvectionCoefficient", h, "AmbientTemperature", T_inf); % Convection Occurs on Loaf

% Parabolic Interpolation to Find Maximum
while error > error_tol
    iter = iter + 1;

    % Calculate New Max Power
    if iter == 1
        P_microwave = guess_vec(1, 1);
    elseif iter < 4
        P_microwave = guess_vec(iter, 1);
        guess_vec(iter - 1, 2) = tastiness;
    elseif iter == 4
        guess_vec(3, 2) = tastiness;
        P_microwave = (guess_vec(1, 2)*(guess_vec(2, 1)^2 - guess_vec(3, 1)^2) + guess_vec(2, 2)*(guess_vec(3, 1)^2 - guess_vec(1, 1)^2) + guess_vec(3, 2)*(guess_vec(1, 1)^2 - guess_vec(2, 1)^2))/(2*(guess_vec(1, 2)*(guess_vec(2, 1) - guess_vec(3, 1)) + guess_vec(2, 2)*(guess_vec(3, 1) - guess_vec(1, 1)) + guess_vec(3, 2)*(guess_vec(1, 1) - guess_vec(2, 1))));
    else
        P_microwave = (guess_vec(1, 2)*(guess_vec(2, 1)^2 - guess_vec(3, 1)^2) + guess_vec(2, 2)*(guess_vec(3, 1)^2 - guess_vec(1, 1)^2) + guess_vec(3, 2)*(guess_vec(1, 1)^2 - guess_vec(2, 1)^2))/(2*(guess_vec(1, 2)*(guess_vec(2, 1) - guess_vec(3, 1)) + guess_vec(2, 2)*(guess_vec(3, 1) - guess_vec(1, 1)) + guess_vec(3, 2)*(guess_vec(1, 1) - guess_vec(2, 1))));
    end

    % Set Volumetric Heat Generation
    internalHeatSource(model, P_skillet/V_skillet, "Cell", 2);
    internalHeatSource(model, P_microwave*epsilon_m_crust/(V_crust*epsilon_m_crust + V_filling*epsilon_m_filling), "Cell", 1);
    internalHeatSource(model, P_microwave*epsilon_m_filling/(V_crust*epsilon_m_crust + V_filling*epsilon_m_filling), "Cell", 3);
    
    % Generate Mesh for FEA
    generateMesh(model);
    
    % Solve with FEA
    solution = solve(model, t_i:t_step:t_f);
    disp("Solved Iteration: " + iter);
    
    % Obtain & Format Temperature Distribution
    [X, Y, Z] = meshgrid(-0.045:d_step:0.045, -0.105:d_step:0.105, 0:d_step:0.045);
    V = interpolateTemperature(solution, X, Y, Z, 61);
    V = reshape(V, size(X));
    
    % Calculate Tastiness
    T_integral = 0;
    [m, n, o] = size(V);
    
    for i = 1:m
        for j = 1:n
            for k = 1:o
                T = V(i, j, k);
                x = X(i, j, k);
                y = Y(i, j, k);
                z = Z(i, j, k);
    
                if z < 0.005 || isnan(T)
                    continue;
                end
    
                if (y >= -0.95 || y <= 0.95) && sqrt(x^2 + y^2) <= 0.034641 && z >= 0.005
                    T_integral = T_integral + rho_filling*((T - T_ideal)/T_tol)^4*d_step^3;
                else
                    T_integral = T_integral + rho_crust*((T - T_ideal)/T_tol)^4*d_step^3;
                end
            end
        end
    end
    
    tastiness = (1 + T_integral/(V_filling*rho_filling + V_crust*rho_crust))^(-1);
    disp("Tastiness Score for Iteration: " + tastiness);

    % Update Guess Vector & Error
    if iter > 3
        if P_microwave < guess_vec(2, 1)
            guess_vec(3, :) = guess_vec(2, :);
            error = (guess_vec(3, 1) - P_microwave)/P_microwave;
        else
            guess_vec(1, :) = guess_vec(2, :);
            error = (P_microwave - guess_vec(1, 1))/P_microwave;
        end

        guess_vec(2, :) = [P_microwave, tastiness];
    end
end

% Plot 3D Contour of Temperature for Maximum Tastiness
figure;
colormap jet;
contourslice(X, Y, Z, V, [], -0.105:0.01:0.105, [])
xlabel("x-pos (m)");
ylabel("y-pos (m)");
zlabel("z-pos (m)");
colorbar;
view(-62,34);
axis equal;

% Plot Heat Map of Outer Mesh for Maximum Tastiness
figure;
pdeplot3D(solution.Mesh, ColorMapData=solution.Temperature(:, 61))

% Output Maximum Tastiness & Corresponding Power
disp("Maximum tastiness score of" + tastiness + " found at microwave power of " + P_microwave);
disp("Error of " + 100*error + "%");
