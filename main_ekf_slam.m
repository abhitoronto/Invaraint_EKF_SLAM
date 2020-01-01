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
% pixel_means = [-1.067; -0.362; -1.163; -0.304];
pixel_means = [0;0;0;0];

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
         t_green = [t_green, k];
    else
         count_red  = [count_red, count];
         t_red = [t_red, k];
    end
end
figure;
scatter(t_green, count_green, 'g.'); hold on;
scatter(t_red, count_red, 'r.');
ylabel("Number of Landmarks");
xlabel("Time [s]");

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

%%% ALL
% first_first = 1;
% last_last = size(y_k_j, 2);

%%% Custom
first_first = 1;
last_last = 500;

%%%%% Plot Ground Truth points
points = [];
for k = first_first:last_last
    C = gt_se3{k}(1:3, 1:3);
    r = -C'*gt_se3{k}(1:3, 4);
    points = [points, r];
end
figure;
plot3(points(1, :), points(2, :), points(3, :)); hold on;
scatter3(P_i(1,:), P_i(2,:), P_i(3,:), 'g'); hold off;
title("GT plot")
grid on;
axis equal;
xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');
%%%%% END: Plot Ground Truth points

land_cov = 10^-6; %m

T_final{first_first} = gt_se3{first_first};
p_final{first_first} = P_i + randn(size(P_i))*sqrt(land_cov*400);
P_final{first_first} = land_cov*400*eye(6 + 3*(num_L),6 + 3*(num_L));
P_final{first_first}(1:6, 1:6) = P_final{first_first}(1:6, 1:6) / 400.0;
time = 0;

% init poses for plotting 
C = T_final{first_first}(1:3, 1:3);
r = -C'*T_final{first_first}(1:3, 4);
points = [r];
for k = first_first+1:last_last
    
    disp("Current Waypoint: " + num2str(k));

    tic;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% START OF MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Propogration
    T_predicted = odom_se3{k-1} * T_final{k - 1};
    p_predicted = p_final{k-1};
    F_k_ = eye(6 + 3*(num_L));
    F_k_(1:6, 1:6) = se3_Ad(odom_se3{k-1});
    Q_k = land_cov*eye(6 + 3*(num_L))*((t(k)-t(k-1))^2);
    Q_k(1:6, 1:6) = diag(odom_cov*((t(k)-t(k-1))^2));
    
    P_predicted = F_k_ * P_final{k-1} * F_k_' + Q_k;
    
    % Estimation
    y_k = [];
    G_k = [];
    R_k = [];
    g_k = [];
    landmarks = [];
    for j = 1:size(y_k_j,3)
        % Loop all valid measurements
        if all(y_k_j(:,k,j) ~= -1)
             landmarks = [landmarks; j];
             y_k = [y_k; y_k_j(:,k,j)];
             
             g_j_k = measurement_model(D* T_c * T_predicted * p_predicted(:,j));
             g_k = [g_k; g_j_k + pixel_means];
             
             G_1_j_k = measurement_jacobian(D* T_c * T_predicted * p_predicted(:,j)) ...
                              * (D* T_c * se3_cd(T_predicted * p_predicted(:,j)));
             G_2_j_k = measurement_jacobian(D* T_c * T_predicted * p_predicted(:,j)) ...
                              * (D* T_c * T_predicted * D');
             G_j_k = [G_1_j_k];
             G_j_k(:, size(G_1_j_k, 2) + size(G_2_j_k, 2)*num_L) = 0;
             G_j_k(:, size(G_1_j_k, 2) + size(G_2_j_k, 2)*(j-1)+1 : ...
                      size(G_1_j_k, 2) + size(G_2_j_k, 2)*j) = G_2_j_k;
             G_k = [G_k; G_j_k];
             
             R_k = [R_k; y_var];
        end
    end
    if length(R_k) > 0
        R_k = diag(R_k);

        temp = G_k * P_predicted * G_k' + R_k;
        K_k = (P_predicted * G_k') / temp;

        % Update
        temp1 = K_k*G_k;
        P_final{k} = (eye(size(temp1)) - temp1) * P_predicted;
        innovation = K_k * (y_k - g_k);
        
        [T_inn, p_inn] = extract_innovation(innovation);
        
        T_final{k} = se3_exp(T_inn) * T_predicted;
        p_final{k} = p_predicted + p_inn;
    else
        P_final{k} = P_predicted;
        T_final{k} = T_predicted;
        p_final{k} = p_predicted;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% END OF MAIN LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if show_plot
        %%%%%%%% Plot Resultant poses
        C = T_final{k}(1:3, 1:3);
        r = -C'*T_final{k}(1:3, 4);
        points = [points, r];
        figure(10);
        plot3(points(1, :), points(2, :), points(3, :)); hold on;
        scatter3(P_i(1,:), P_i(2,:), P_i(3,:), 'g'); hold on;
        scatter3(p_final{k}(1,:), p_final{k}(2,:), p_final{k}(3,:), 'r'); hold off;
        grid on;
        title("Resultant Poses")
        axis equal;
        xlabel('X [m]');
        ylabel('Y [m]');
        zlabel('Z [m]');
        %%%%%%%% END: plot resultant poses
        
%         %%%%%%%% VISUALIZE A
%         figure(11);
%         spy(P_final{k})
%         title('Covariance matrix')
    end
    
    time = time + toc;
end

%%%%% Calculate Errors
%%%%%%% Pose Errors
r_error = [];
C_error = [];
p_error = [];
state_nees = [];
r_nees = [];
C_nees = [];
p_nees = [];
covs = [];
for k = first_first:last_last
    % Error Calculations
    C_op = T_final{k}(1:3, 1:3);
    r_op = -C_op'*T_final{k}(1:3, 4);
    C_gt = gt_se3{k}(1:3, 1:3);
    r_gt = -C_gt'*gt_se3{k}(1:3, 4);
    diff = (r_gt - r_op);
    r_error = [r_error, diff];
    C_error = [C_error, so3_log(C_op*C_gt')];
    p_error = [p_error, rms(P_i - p_final{k}, 2)];
    p_error_l = P_i - p_final{k};
    state_error = [r_error(:, end); C_error(:, end); reshape(p_error_l(1:3,:), [], 1) ];
    state_nees = [state_nees, state_error' * inv(P_final{k}) * state_error];
    r_nees = [r_nees, r_error(:, end)' * inv(P_final{k}(1:3, 1:3)) * r_error(:, end)];
    C_nees = [C_nees, C_error(:, end)' * inv(P_final{k}(4:6, 4:6)) * C_error(:, end)];
    p_nees = [p_nees, state_error(7:end)' * inv(P_final{k}(7:end, 7:end)) * state_error(7:end)];
    [cov_pose, cov_landmarks] = extract_innovation(diag(P_final{k}));
    cov_landmarks = mean(cov_landmarks, 2);
    covs = [covs, [cov_pose; cov_landmarks]];
end
mean_error = mean(r_error, 2);
std_error = std(r_error, 0, 2);
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
scatter3(P_i(1,:), P_i(2,:), P_i(3,:), 'g'); hold on;
scatter3(p_final{end}(1,:), p_final{end}(2,:), p_final{end}(3,:), 'r'); hold off;
grid on;
title("Resultant Poses");
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

% Plot Landmark Errors
figure;
subplot(3,1,1);
plot(t(first:last), covs(7,:), '--r', t(first:last), -covs(7,:), '--r', t(first:last), p_error(1,:));
xlabel("Time steps");
ylabel('\deltap_{x,k}');
legend('3\sigma bound');
grid on;
subplot(3,1,2);
plot(t(first:last), covs(8,:), '--r', t(first:last), -covs(8,:), '--r', t(first:last), p_error(2,:));
xlabel("Time steps");
ylabel('\deltap_{y,k}');
legend('3\sigma bound');
grid on;
subplot(3,1,3);
plot(t(first:last), covs(9,:), '--r', t(first:last), -covs(9,:), '--r', t(first:last), p_error(3,:));
xlabel("Time [s]");
ylabel('\deltap_{z,k}');
legend('3\sigma bound');
grid on;

%%%%%% NEES plot
figure;
plot(t(first:last), state_nees, t(first:last), ones(size(t(first:last))), 'r');
xlabel("Time [s]");
ylabel('Normalized estimation error squared');
grid on;

figure;
subplot(3,1,1);
plot(t(first:last), r_nees, t(first:last), ones(size(t(first:last))), 'r');
xlabel("Time [s]");
ylabel('Position NEES');
grid on;
subplot(3,1,2);
plot(t(first:last), C_nees, t(first:last), ones(size(t(first:last))), 'r');
xlabel("Time [s]");
ylabel('Orientation NEES');
grid on;
subplot(3,1,3);
plot(t(first:last), p_nees, t(first:last), ones(size(t(first:last))), 'r');
xlabel("Time [s]");
ylabel('Landmark Position NEES');
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

function [pose_i, landmark_i] = extract_innovation(innovation)
    
    assert(length(innovation) == 66, "ERROR: extarct_innovation: " + mat2str(innovation))
    
    pose_i = innovation(1:6);
    pose_i = reshape(pose_i, 6, 1);
    landmark_i = innovation(7:end);
    
    landmark_i = reshape(landmark_i, 3, 20);
    landmark_i(4,1) = 0;
end