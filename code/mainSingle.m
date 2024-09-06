
%% set initial parameter
tau0 = 6667; 
patient = 'P1';
EZ=64
connectivity_path = strcat('../data/connectivity_', patient, '/');
K = load(strcat(connectivity_path, 'weights.txt')); % Load connectivity matrix
K = normal(K); % Normalize
[N, ~] = size(K); % Get the number of nodes
x0 = -2.3 + zeros(N, 1); 
x0(EZ)=-2.3

%% solve analytic parameter
 % a = 1;
 % b = 8 * x0 + 8 / 3;
 % c = (4 * x0 + 16 / 3).^2 + 629.6 / 27; 
 % D = b.^2 - 4 * a * c;
 % Z1 = NaN(N, 1); 
 % Z2 = NaN(N, 1); 
 % Z1_valid = false(N, 1);
 % Z2_valid = false(N, 1);
 % z_analytic_fixed_point = NaN(N, 1); % initialize zfixed_point
 % for i = 1:N
     % if D(i) >= 0
         % Z1(i) = (-b(i) + sqrt(D(i))) / (2 * a);
         % Z2(i) = (-b(i) - sqrt(D(i))) / (2 * a);
         % Z1_valid(i) = (8 * Z1(i) > 629.6 / 27) && (Z1(i) <= -4 * x0(i) - 16 / 3);
         % Z2_valid(i) = (8 * Z2(i) > 629.6 / 27) && (Z2(i) <= -4 * x0(i) - 16 / 3);    
         % if Z1_valid(i)
             % z_analytic_fixed_point(i) = Z1(i);
         % elseif Z2_valid(i)
             % z_analytic_fixed_point(i) = Z2(i);
         % end
     % end
     % fprintf('Root 1: %.4f (Valid: %d)\n', Z1(i), Z1_valid(i));
     % fprintf('Root 2: %.4f (Valid: %d)\n', Z2(i), Z2_valid(i));
     % fprintf('zfixed_point: %.4f\n', z_analytic_fixed_point(i));
     % fprintf('----------------------------------\n');
 % end


%% solve numerical-fixed-point
Z0 = 3 + zeros(N, 1); 
one_dim_epileptor_fun = @(z) oneDepileptor(z, x0, K, tau0); 
opt = optimset('TolFun', 1e-14, 'TolX', 1e-14);
[z_fixed_numerical, fval, exitflag, output] = fsolve(one_dim_epileptor_fun, Z0, opt);
z_fixed_numerical = real(z_fixed_numerical);
%- print result
% for i = 1:length(z_fixed_numerical)
    % disp(['z_fixed_numerical at N = ', num2str(i), ' is ', num2str(z_fixed_numerical(i))])
% end


%% Create plot for z_fixed_numerical
% figure; % Create a new figure
% plot(1:N, z_fixed_numerical, '-o'); % Plot the values
% title(['x0(EZ)=(', num2str(x0(EZ)), ')']);
% xlabel('Nodes');
% ylabel('z\_fixed\_numerical');
% grid on; % Add grid for better readability


%% Calculate Coupling Matrix
C = CouplingMatrix(z_fixed_numerical , K, tau0);


%% Calculate Synchronizability
[sortedLambdasCU, UcovarianceU, err] = covarianceUGaussianNet(C, false, 10000000, false, 1);
syncWidth = synchronizability(UcovarianceU); 
fprintf('Synchronization Width: %.3f\n', syncWidth);

%% Calculate Stability
[covariance,err] = con2cov(C,false,100000,10,true)
stability = Stability(covariance); 
fprintf('Stability: %.3f\n', stability);

%% Calculate 'how much is node i being pushed away from stability'
% BeAwayStability=NaN(1,N)
% trace=diag(covariance)
% for i=1:N
    % BeAwayStability(i)=trace(i)
% end


%% Create plot for BeAwayStability
% figure; % Create a new figure
% plot(1:N, BeAwayStability, '-o'); % Plot the values
% title(['x0(EZ)=(', num2str(x0(EZ)), ')']);
% xlabel('Nodes');
% ylabel('BeAwayStability');
% grid on; % Add grid for better readability


%% calculate 'how much node j is pushing all nodes of the network away' from stability
% PushAwayStability=NaN(1,N)
% C_transpose=C'
% [covariance_new,err] = con2cov(C_transpose,false,100000,10,1)
% trace_new=diag(covariance_new)
% for i=1:N
    % PushAwayStability(i)=trace_new(i)
% end

%% Create plot for PushAwayStability
% figure; % Create a new figure
% plot(1:N, PushAwayStability, '-o'); % Plot the values
% title(['x0(EZ)=(', num2str(x0(EZ)), ')']);
% xlabel('Nodes');
% ylabel('PushAwayStability');
% grid on; % Add grid for better readability

