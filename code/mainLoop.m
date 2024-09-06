diary('today.txt')

%% clear parameter
clear
clc

%% set basic parameter
tau0 = 6667; 
patient = 'P1';
EZ = 64;
connectivity_path = strcat('../data/connectivity_', patient, '/');
K = load(strcat(connectivity_path, 'weights.txt')); % Load connectivity matrix
K = normal(K); % Normalize
[N, ~] = size(K); % Get the number of nodes

%%%%%%%%%%%%%%%% this is for analytic fixed-point range %%%%%%%%%%%%%%%%\

%% initialize
%x0_range = -2.30;
%x0_EZ_range = -1.9:0.1:-1.67
%syncWidthMatrix = NaN(length(x0_range), length(x0_EZ_range)); % save syncWidth into matrix
%StabilityMatrix = NaN(length(x0_range), length(x0_EZ_range)); % save Stability into matrix
%LargestLambdaMatrix = NaN(length(x0_range), length(x0_EZ_range)); % save Largest eigenvalue into matrix
%ZiAnalytic = NaN(length(x0_range), length(x0_EZ_range)); % save z_i-analytic

%% loop x0 and x0(EZ)
%for ix0 = 1:length(x0_range)
    %for ix0_EZ = 1:length(x0_EZ_range)
        %x0 = x0_range(ix0) + zeros(N, 1); 
        %x0(EZ) = x0_EZ_range(ix0_EZ);
        
        %% solve fixed-point
        %a = 1;
        %b = 8 * x0 + 8 / 3;
        %c = (4 * x0 + 16 / 3).^2 + 629.6 / 27; 
        %D = b.^2 - 4 * a * c; 
        %Z1 = NaN(N, 1); % initialize as NaN vector
        %Z2 = NaN(N, 1)
        %Z1_valid = false(N, 1); % initialize as false vector
        %Z2_valid = false(N, 1); % initialize as false vector
        %z_analytic_fixed_point = NaN(N, 1); % initialize zfixed_point
        
        %for i = 1:N
            %if D(i) >= 0
                %Z1(i) = (-b(i) + sqrt(D(i))) / (2 * a);
                %Z2(i) = (-b(i) - sqrt(D(i))) / (2 * a);
                %% Check for additional constraints
                %Z1_valid(i) = (8 * Z1(i) > 629.6 / 27) && (Z1(i) <= -4 * x0(i) - 16 / 3);
                %Z2_valid(i) = (8 * Z2(i) > 629.6 / 27) && (Z2(i) <= -4 * x0(i) - 16 / 3);    
                %% Set Valid fixed-point
                %if Z1_valid(i)
                    %z_analytic_fixed_point(i) = Z1(i);
                %elseif Z2_valid(i)
                    %z_analytic_fixed_point(i) = Z2(i);
                %end
            %end
            
        %end
        
        %% calculate coupling matrix
        %C = CouplingMatrix(z_analytic_fixed_point, K, tau0);
        
        %% calculate synchronizability
        %[sortedLambdasCU, UcovarianceU, err] = covarianceUGaussianNet(C, false, 100000, false, 1);
        %syncWidth = synchronizability(UcovarianceU); 
                
        %% save syncWidth into matrix
        %syncWidthMatrix(ix0, ix0_EZ) = syncWidth;
        %fprintf('Synchronization Width: %.3f\\n', syncWidth);

        %% calculate Stability
        %[covariance,err] = con2cov(C,false,100000,10,1)
        %stability = Stability(covariance); 
       
        %% save Stability into matrix
        %StabilityMatrix(ix0, ix0_EZ) = stability;
        %fprintf('Stability: %.3f\\n', stability);
        
        %% save zi of analytic fixed-point
        %ZiAnalytic(ix0, ix0_EZ)=z_analytic_fixed_point
        %fprintf('ZiAnalytic: %.3f\\n', z_analytic_fixed_point);

        %% calculate and save largest eigenvalue
        %eigen=eig(C)
        %largestEigen=max(eigen)
        %fprintf('largest eigenvalue is: %.3f\\n', largestEigen);
        %LargestLambdaMatrix(ix0, ix0_EZ) = largestEigen

            
    %end
%end


%%%%%%%%%%%%%%%% this is for numerical fixed-point range %%%%%%%%%%%%%%%% \

%%- initialize
format long
x0_range1 = -2.30;
%x0_EZ_range1 = linspace(-1.672548, -1.672545, 10);
x0_EZ_range1 = -2.3:0.05:-1.6
disp(x0_EZ_range1) 
syncWidthMatrix1 = NaN(length(x0_range1), length(x0_EZ_range1)); 
StabilityMatrix1 = NaN(length(x0_range1), length(x0_EZ_range1)); 
TotalCouplingMatrix= NaN(length(x0_range1), length(x0_EZ_range1),N,N) 
LargestLambdaMatrix1 = NaN(length(x0_range1), length(x0_EZ_range1)); 
BeAwayStabilityMatrix= NaN(length(x0_range1), length(x0_EZ_range1),N) 
PushAwayStabilityMatrix= NaN(length(x0_range1), length(x0_EZ_range1),N) 
NumericFixedPointMatrix = NaN(length(x0_range1),length(x0_EZ_range1),N)

%%-dual loop x0 and x0(EZ)
for ix0 = 1:length(x0_range1)
    for ix0_EZ = 1:length(x0_EZ_range1)
        x0=x0_range1(ix0)+zeros(N,1);
        x0(EZ) = x0_EZ_range1(ix0_EZ);
        
        %%- solve fixed-point
        Z0 = 3 + zeros(N, 1); 
        one_dim_epileptor_fun = @(z) oneDepileptor(z, x0, K, tau0); 
        opt = optimset('TolFun', 1e-14, 'TolX', 1e-14);
        [z_fixed_numerical, fval, exitflag, output] = fsolve(one_dim_epileptor_fun, Z0, opt);
         

        %%- save zi for numerical fixed-point
        NumericFixedPointMatrix(ix0, ix0_EZ) = z_fixed_numerical(EZ);
        
        %%- calculate Coupling Matrix for numerical solution
        C1 = CouplingMatrix(z_fixed_numerical, K, tau0);
        TotalCouplingMatrix(ix0, ix0_EZ,:,:) = C1;
        
        %%- calcaulate sycnchronizability for numerical solution
        [sortedLambdasCU1, UcovarianceU1, err] = covarianceUGaussianNet(C1, false, 100000, false, 1);
        syncWidth1 = synchronizability(UcovarianceU1); 
                
        %%- save syncWidth into matrix
        syncWidthMatrix1(ix0, ix0_EZ) = syncWidth1;
        fprintf('Synchronization Width1: %.3f\\n', syncWidth1);
       

        %%- calculate Stability
        [covariance1,err] = con2cov(C1,false,100000,10,1)
        stability1 = Stability(covariance1); 
       
        %%- save Stability into matrix
        StabilityMatrix1(ix0, ix0_EZ) = stability1;
        fprintf('Stability1: %.3f\\n', stability1);


        %%- calculate and save largest eigenvalue
        eigen1=eig(C1)
        eigen1=real(eigen1)
        largestEigen1=max(eigen1)
        LargestLambdaMatrix1(ix0, ix0_EZ, :) = largestEigen1;
        fprintf('largest eigenvalue1 is: %.3f\\n', largestEigen1);
        
        
        %%- calculate BeAwayStability
        BeAwayStability = NaN(1, N);
        trace = diag(covariance1);
        for i = 1:N
            BeAwayStability(i) = trace(i);
        end
        BeAwayStabilityMatrix(ix0, ix0_EZ, :) = BeAwayStability;
        fprintf('BeAwayStability: %.3f\\n', BeAwayStability(EZ));

        %%- calculate PushAwayStability
        C_transpose = C1';
        [covariance_new, err] = con2cov(C_transpose, false, 100000, 10, 1);
        trace_new = diag(covariance_new);
        PushAwayStability = NaN(1, N);
        for i = 1:N
            PushAwayStability(i) = trace_new(i);
        end
        PushAwayStabilityMatrix(ix0, ix0_EZ, :) = PushAwayStability;
        fprintf('PushAwayStability: %.3f\\n', PushAwayStability(EZ));


    end
end


%% - plot out details for coupling matrix
% C1_values_EZ_EZ = NaN(length(x0_EZ_range1), 1);
% C1_values_other_EZ = NaN(length(x0_EZ_range1), N-1); 
% C1_sum_col_EZ = NaN(length(x0_EZ_range1), 1);

% for ix0 = 1:length(x0_range1)
    % for ix0_EZ = 1:length(x0_EZ_range1) 
	% C1 = squeeze(TotalCouplingMatrix(ix0, ix0_EZ, :, :)); 
	% C1_values_64_64(ix0_EZ) = C1(EZ, EZ); 
	% C1_values_other_EZ(ix0_EZ,:)=C1([1:EZ-1,EZ+1:N],EZ)
	% C1_sum_col_EZ(ix0_EZ) = sum(C1(:, EZ));
    % end
% end

% - plot all lines together
% figure; 
% hold on; 
% plot(x0_EZ_range1, real(C1_values_EZ_EZ), 'LineWidth', 2, 'DisplayName', 'C1(EZ, EZ)\'92); plot(x0_EZ_range1, real(C1_sum_col_EZ), 'LineWidth', 2, 'DisplayName', 'sum(C1(:, EZ))\'92); 
% for node = 1:(N-1) 
      % plot(x0_EZ_range1, real(C1_values_other_EZ(:, node)), 'LineWidth', 1, 'DisplayName', ['Node ' num2str(node)]); 
% end 
% hold off; 
% xlabel('x0(EZ)'); 
% ylabel('Values'); 
% title('Details for CouplingMatrix'); legend('show');

%- plot subgraph version

% figure;

% subplot(3, 1, 1);
% plot(x0_EZ_range1, real(C1_values_EZ_EZ), 'LineWidth', 2);
% xlabel('x0(EZ)');
% ylabel('C(EZ,EZ)');
% title('Self-Coupling');
% legend('C(EZ, EZ)');

% subplot(3, 1, 2);
% hold on;
% for node = 1:(N-1)
    % plot(x0_EZ_range1, real(C1_values_other_EZ(:, node)), 'LineWidth', 1);
% end
% hold off;
% xlabel('x0(EZ)');
% ylabel('C(other nodes, EZ)');
% title('Cross-Coupling');
% legend('C(other nodes, EZ)');

% subplot(3, 1, 3);
% plot(x0_EZ_range1, real(C1_sum_col_EZ), 'LineWidth', 2);
% xlabel('x0(EZ)');
% ylabel('sum(C(:, EZ))');
% title('Sum of Self and Cross Coupling to EZ');
% legend('sum(C(:, EZ))');



%%- turn off diary
diary off

function w = normal(W)
    [N, ~] = size(W);
    u = zeros(N * (N - 1) / 2, 1);
    k = 0;
    for i = 1:N
        for j = i + 1:N
            k = k + 1;
            u(k, 1) = W(i, j);
        end
    end
    u = sort(u);
    uu = u(floor(length(u) * 0.95));
    w = W - diag(diag(W));
    idx = w > uu;
    w(idx) = uu;
    w = w / (max(max(w)));
end
