%% Assignment 3: Starry Night Dataset
%%% 
close all
clear all

% Add paths
addpath('./lie_functions');

%% Problem Setup
global fu fv cu cv b

% Load starry night data set
load('dataset3.mat');

% Constants
show_plot = 0;
D = [eye(3), zeros(3,1)];
num_L = size(rho_i_pj_i, 2);

% Get Odom in se3 algebra (remove the influence of dt)
% NOTE: w is being negated. Due to eq 6.45 in Barfoot book.
odom_combined = get_combined_odom(v_vk_vk_i, w_vk_vk_i, t);

odom_cov = [ v_var; w_var ];
odom_se3 = get_se3_array(odom_combined);

% Get groudtruth is se3 
gt_combined = [r_i_vk_i; theta_vk_i]; 
gt_se3 = get_se3_array(gt_combined);

% Convert landmark vectors to be one ending vectors
P_i = append_1(rho_i_pj_i); %clear rho_i_pj_i;

% Create the camera transform matrix
T_c = [ C_c_v      , -C_c_v*rho_v_c_v;
        zeros(1,3) , 1                 ];
% clear rho_v_c_v C_c_v;
    
%% Plot number of landmarks vs time
t_green = [];
count_green = [];
t_red = [];
count_red = [];

for k = 1:size(y_k_j,2)
    count = 0;
    for j = 1:size(y_k_j,3)
        if all(y_k_j(:,k,j) ~= -1)
            count = count + 1;
        end
    end
    if count >= 3
         count_green  = [count_green, count];
         t_green = [t_green, t(k)];
    else
         count_red  = [count_red, count];
         t_red = [t_red, t(k)];
    end
end
figure;
scatter(t_green, count_green, 'g.'); hold on;
scatter(t_red, count_red, 'r.');
ylabel("Number of Landmarks");
xlabel("Time [s]")

%% SANITY
% for k = 1000:1010
%     for j = 1:20
%         if all(y_k_j(:, k, j) ~= -1)
%             g1 = measurement_model( D* T_c * gt_se3{k} * P_i(:,j));
%             g2 = measurement_model(C_c_v*( so3_exp(theta_vk_i(:, k))*(rho_i_pj_i(:,j) - r_i_vk_i(:,k)) - rho_v_c_v));
%             diff1 = mean(g1 - y_k_j(:, k, j));
%             diff2 = mean(g2 - y_k_j(:, k, j));
%             disp(mat2str( [diff1; diff2]))
%         end 
%     end
% end

% Plot dead-reckoning
% T_prior = {};
% points = []; % Plot points
% T_prior{1} = gt_se3{1};
% for k = 2:length(gt_se3)
%     % Calculate prior through motion model
%     T_prior{k} = odom_se3{k-1} * T_prior{k-1};
%     C = T_prior{k}(1:3, 1:3);
%     r = -C'*T_prior{k}(1:3, 4);
%     points = [points, r];
% end
% figure;
% title("Dead Reckoning plot")
% plot3(points(1, :), points(2, :), points(3, :)); hold on;
% scatter3(P_i(1,:), P_i(2,:), P_i(3,:), 'r'); hold off;
% grid on;
% axis equal;
% xlabel('X [m]');
% ylabel('Y [m]');
% zlabel('Z [m]');


%% MAP setup

first_first = 1;
last_last = size(y_k_j, 2);
% last_last = 300;

%%%%% Plot Ground Truth points
if show_plot
    points = [];
    for k = first_first:last_last
        C = gt_se3{k}(1:3, 1:3);
        r = -C'*gt_se3{k}(1:3, 4);
        points = [points, r];
    end
    figure;
    plot3(points(1, :), points(2, :), points(3, :)); hold on;
    scatter3(P_i(1,:), P_i(2,:), P_i(3,:), 'r'); hold off;
    title("GT plot")
    grid on;
    axis equal;
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    %%%%% END: Plot Ground Truth points
end


T_final{first_first} = gt_se3{first_first};
P_final{first_first} = 0.0001*eye(6,6);
time = 0;
% init points
C = T_final{1}(1:3, 1:3);
r = -C'*T_final{1}(1:3, 4);
points = [r];
for k = first_first+1:last_last
    
    disp("Current Waypoint: " + num2str(k));

    tic;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% START OF MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Propogration
    T_predicted = odom_se3{k-1} * T_final{k - 1};
    F_k_ = se3_Ad(odom_se3{k-1});
    
    P_predicted = F_k_ * P_final{k-1} * F_k_' + diag(odom_cov*((t(k)-t(k-1))^2));
    
    % Estimation
    y_k = [];
    G_k = [];
    R_k = [];
    g_k = [];
    for j = 1:size(y_k_j,3)
        % Loop all valid measurements
        if all(y_k_j(:,k,j) ~= -1)
             y_k = [y_k; y_k_j(:,k,j)];
             
             g_j_k = measurement_model(D* T_c * T_predicted * P_i(:,j));
             g_k = [g_k; g_j_k];
             
             G_j_k = measurement_jacobian(D* T_c * T_predicted * P_i(:,j)) ...
                              * (D* T_c * se3_cd(T_predicted * P_i(:,j)));
             G_k = [G_k; G_j_k];
             
             R_k = [R_k; y_var];
        end
    end
    if length(R_k) > 0
        R_k = diag(R_k);

        temp = sparse(G_k * P_predicted * G_k' + R_k);
        K_k = (P_predicted * G_k') / temp;

        % Update
        temp1 = K_k*G_k;
        P_final{k} = (eye(size(temp1)) - temp1) * P_predicted;
        T_final{k} = se3_exp( K_k * (y_k - g_k)) * T_predicted;
    else
        P_final{k} = P_predicted;
        T_final{k} = T_predicted;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% END OF MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if show_plot
        %%%%%%%% Plot Resultant poses
        C = T_final{k}(1:3, 1:3);
        r = -C'*T_final{k}(1:3, 4);
        points = [points, r];
        figure(10);
        plot3(points(1, :), points(2, :), points(3, :)); hold on;
        scatter3(P_i(1,:), P_i(2,:), P_i(3,:), 'r'); hold off;
        grid on;
        title("Resultant Poses")
        axis equal;
        xlabel('X [m]');
        ylabel('Y [m]');
        zlabel('Z [m]');
        %%%%%%%% END: Plot Resultant poses
    end
    
    time = time + toc;
end
    
r_error = [];
C_error = [];
covs = [];
for k = first_first:last_last
    % Error Calculations
    C_op = T_final{k}(1:3, 1:3);
    r_op = -C_op'*T_final{k}(1:3, 4);
    C_gt = gt_se3{k}(1:3, 1:3);
    r_gt = -C_gt'*gt_se3{k}(1:3, 4);
    diff = (r_gt - r_op);
    r_error = [r_error, diff];
    C_error = [C_error, so3_hatinv(eye(3) - C_op*C_gt')];
    covs = [covs, diag(P_final{k})];
end
mean_error = mean(r_error, 2);
std_error = std(r_error, 0,2);
disp("Mean Positional error is: " + num2str(mean_error));
disp("Positional error std is: " + num2str(std_error));
disp("Total Time: " + num2str(time));

%%%%%%%% Plot Resultant poses
points = [];
for k = first_first:last_last
    C = T_final{k}(1:3, 1:3);
    r = -C'*T_final{k}(1:3, 4);
    points = [points, r];
end
figure;
plot3(points(1, :), points(2, :), points(3, :)); hold on;
scatter3(P_i(1,:), P_i(2,:), P_i(3,:), 'r'); hold off;
grid on;
title("Resultant Poses")
axis equal;
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
%%%%%%%% END: Plot Resultant poses

first = first_first;
last = last_last;

covs = sqrt(covs)*3;

% Plot Angle Error
figure;
subplot(3,1,1);
plot(t(first:last), covs(4,:), '--r', t(first:last), -covs(4,:), '--r', t(first:last), C_error(1,:));
xlabel("Time steps");
ylabel('\delta\theta_{x,k}');
legend('3\sigma bound');
grid on;
subplot(3,1,2);
plot(t(first:last), covs(5,:), '--r', t(first:last), -covs(5,:), '--r', t(first:last), C_error(2,:));
xlabel("Time steps");
ylabel('\delta\theta_{y,k}');
legend('3\sigma bound');
grid on;
subplot(3,1,3);
plot(t(first:last), covs(6,:), '--r', t(first:last), -covs(6,:), '--r', t(first:last), C_error(3,:));
xlabel("Time steps");
ylabel('\delta\theta_{z,k}');
legend('3\sigma bound');
grid on;

% Plot Position Errors
figure;
subplot(3,1,1);
plot(t(first:last), covs(1,:), '--r', t(first:last), -covs(1,:), '--r', t(first:last), r_error(1,:));
xlabel("Time steps");
ylabel('\deltar_{x,k}');
legend('3\sigma bound');
grid on;
subplot(3,1,2);
plot(t(first:last), covs(2,:), '--r', t(first:last), -covs(2,:), '--r', t(first:last), r_error(2,:));
xlabel("Time steps");
ylabel('\deltar_{y,k}');
legend('3\sigma bound');
grid on;
subplot(3,1,3);
plot(t(first:last), covs(3,:), '--r', t(first:last), -covs(3,:), '--r', t(first:last), r_error(3,:));
xlabel("Time [s]");
ylabel('\deltar_{z,k}');
legend('3\sigma bound');
grid on;



%% Helper Functions
function vec1 = append_1(X)
    assert( size(X,1)==3 || size(X,2)==3, "append_1: bad value " + mat2str(X));
    X = reshape(X, 3, []);
    one_array = ones(1, size(X,2));
    vec1 = [X; one_array];
end

function jac_jk = measurement_jacobian(X_op)
    global fu fv b;
    assert( numel(X_op) == 3, "measurement_jacobian: " + mat2str(X_op) );
    x = X_op(1);
    y = X_op(2);
    z = X_op(3);
    jac_jk = [  fu/z,    0,         -(fu*x)/(z^2);
                0,       fv/z,      -(fv*y)/(z^2);
                fu/z,    0,          (fu*(b - x))/(z^2);
                0,       fv/z,      -(fv*y)/(z^2)];
end

function g = measurement_model(X_jk)
    
    global fu fv cu cv b;
    assert( numel(X_jk) == 3, "measurement_jacobian: " + mat2str(X_jk) );
    
    x = X_jk(1);
    y = X_jk(2);
    z = X_jk(3);
    
    ul = fu * x / z + cu;
    vl = fv * y / z + cv;
    ur = fu * (x - b) / z + cu;
    vr = vl;

    g = [ul; vl; ur; vr];
end

function result = get_se3_array(phi)
    %%% Note: phi does not belong to Lie Algebra
    % Phi = [dR and dTheta];

    assert(size(phi,1)==6 || size(phi,2)==6, "get_se3_exp_array: argument error " + mat2str(phi))
    phi = reshape(phi, 6, []);
    result = {};
    for i = 1:size(phi,2)
        C = so3_exp(phi(4:6, i));
        r = phi(1:3, i);
        
        % Formula given in lecture 11 page 15
        result{i} = [C         ,  -C*r;
                     zeros(1,3),   1   ];
    end
end

function [result] = get_combined_odom(v, w, t)

    assert( size(v, 2) == size(w,2), 'v and w must be the same size');
    assert( size(t, 2) == size(v,2), 't, v and w must be the same size');
    
    tV = [v;
          w ];
    
    dt = t(:, 2:end) - t(:, 1:end-1);
    result = tV(:,1:end-1) .* dt;
end

function [result] = project_to_se3(T)
    
    assert(size(T,1) == 4 && size(T,2) == 4, "ERROR: project_to_se3: " + mat2str(T))
    
    result  = T;
    C = T(1:3, 1:3);
    result(1:3, 1:3) = ((C*C')^(-1/2)) * C;
    result(4, 1:3) = zeros(1,3);
    result(4,4) = 1;
end