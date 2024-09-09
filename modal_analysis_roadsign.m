clc;
clear;
% Defining the parameters
nodes = 5;% Number of nodes
elements = 4; % Number of elements
E = 210000000000;% Young's modulus

L =[4.25, 8.25,8.25,4.25];% Array for length of each element
I = [0.1808,0.0302,0.0302,0.1808];% Array for moment of inertia of each element
A = [1.152e-5,1.926e-6,1.926e-6,1.152e-5]; % Array for cross sectional area of each element
rho = 7850;%density
boun = [0,0,0,1,0,1,0,1,1,0,0,1,0,0,0];% Boundary conditions

total_dofs = nodes * 3; % Number of nodesxdegree of freedom at each node
% Initialising the global stiffness matrix with zeros
K_mat = zeros(total_dofs);

% Calculating the stiffness matrix
for i=1:length(L)
    k1=[(E*A(i))/L(i), 0, 0, -(E*A(i))/L(i), 0, 0;
          0, (12*E*I(i))/(L(i)^3), (6*E*I(i))/(L(i)^2), 0, -(12*E*I(i))/(L(i)^3),(6*E*I(i))/(L(i)^2);
          0, (6*E*I(i))/(L(i)^2), (4*E*I(i))/L(i), 0 , -(6*E*I(i))/(L(i)^2),(2*E*I(i))/L(i);
          -(E*A(i))/L(i), 0, 0 , (E*A(i))/L(i), 0, 0;
          0, -(12*E*I(i))/(L(i)^3), -(6*E*I(i))/(L(i)^2),0, (12*E*I(i))/(L(i)^3),-(6*E*I(i))/(L(i)^2);
          0, (6*E*I(i))/(L(i)^2), (2*E*I(i))/L(i),0,-(6*E*I(i))/(L(i)^2),(4*E*I(i))/L(i)];

    % Update the appropriate elements in K_mat
    start_index = (i - 1) * 3 + 1;
    end_index = start_index + 5;
    K_mat(start_index:end_index, start_index:end_index) = K_mat(start_index:end_index, start_index:end_index) + k1;
end
K_mat;

% Generating reduced global stiffness matrix
% Initialising the reduced global stiffness matrix
red_K_mat = K_mat;

% Variable to keep track of deleted rows and columns
deleted = 0;

% Loop through each boundary condition
for i = 1:length(boun)
    % Check if the boundary condition is zero
    if boun(i) == 0
        % Delete the row corresponding to the boundary condition
        red_K_mat(i - deleted, :) = [];
        % Delete the column corresponding to the boundary condition
        red_K_mat(:, i - deleted) = [];
        % Increment the number of deleted rows and columns
        deleted = deleted + 1;
    end
end
K = red_K_mat;

% Initialising the global mass matrix with zeros
M_mat = zeros(total_dofs);

% Calculating the mass matrix
for i = 1:length(L)
    m1= (rho*A(i)*L(i))/420*[140, 0, 0,70,0,0;
                             0,156,22*L(i),0,54,-13*L(i);
                             0,22*L(i),4*(L(i)^2),0,13*L(i),-3*(L(i)^2);
                             70, 0, 0, 140, 0, 0;
                             0, 54, 13*L(i), 0, 156, -22*L(i);
                             0, -13*L(i), -3*(L(i)^2), 0,-22*L(i),4*(L(i)^2)];

    % Update the appropriate elements in M_mat
    start_index = (i - 1) * 3 + 1;
    end_index = start_index + 5;
    M_mat(start_index:end_index, start_index:end_index) = M_mat(start_index:end_index, start_index:end_index) + m1;
end
M_mat;

% Generating reduced mass matrix
% Initialising the reduced mass matrix
red_M_mat = M_mat;

% Variable to keep track of deleted rows and columns
deleted = 0;

% Loop through each boundary condition
for i = 1:length(boun)
    % Check if the boundary condition is zero
    if boun(i) == 0
        % Delete the row corresponding to the boundary condition
        red_M_mat(i - deleted, :) = [];
        % Delete the column corresponding to the boundary condition
        red_M_mat(:, i - deleted) = [];
        % Increment the number of deleted rows and columns
        deleted = deleted + 1;
    end
end
 M = red_M_mat;

%Eigen value and Eigen Vector calculation
[e_vec,e_val] = eig(K,M);

%angular frequency omega(w)
w = sqrt(e_val);
w1 = w(1,1);
w2 = w(2,2);
w3 = w(3,3);
w4 = w(4,4);
w5 = w(5,5);
disp([w1;w2;w3;w4;w5])

%calculating the natural time period
T1 =(2*3.14)/(w1);
T2 =(2*3.14)/(w2);
T3 =(2*3.14)/(w3);
T4 =(2*3.14)/(w4);
T5 =(2*3.14)/(w5);
disp([T1;T2;T3;T4;T5]);

% Initialize array to store normalized eigenvectors
norm_e_vec = zeros(size(e_vec));

% Loop through each eigenvector
for i = 1:size(e_vec, 2)
    % Get the eigenvector
    eigenvector = e_vec(:, i);
    
    % Normalize the eigenvector
    norm_eigenvector = eigenvector / norm(eigenvector);
    
    % Store the normalized eigenvector
    norm_e_vec(:, i) = norm_eigenvector;
end
 %plotting the mode shapes
for i = 1:5
    subplot(1,5,i)
    plot([0;norm_e_vec(:,i);0]);
end

% Modal m matrix
M_modal = transpose(norm_e_vec)*M*norm_e_vec;
disp(M_modal);

% Modal k matrix
K_modal = transpose(norm_e_vec)*K*norm_e_vec;
disp(K_modal);

% Rayleigh's Quotient
R = K_modal/M_modal;
disp(R);



