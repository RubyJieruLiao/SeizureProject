%diary('today.txt')

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


%%- initialize
format long
x0_range = -2;

% x0_EZ_range = -1.76:0.005:-1.64 % for all C's change-big version
% x0_EZ_range = linspace(-1.68, -1.672, 30); % for all C's change-zoom in version, so that we can see transition at different time
% x0_EZ_range = -1.674:0.00002:-1.673 %for critical point "EZ-to-EZ from C"
% x0_EZ_range = -1.6726:0.00002:-1.672 %for critical point "EZ-to-All from C"

% x0_EZ_range = [-1.67254746492271, -1.67254746198135, -1.67254745903999, -1.67254745609863, -1.67254745315727] %for critical point "Stability"

x0_EZ_range = -1.68:0.0001:-1.676

disp(x0_EZ_range) 

NumericFixedPointMatrix = NaN(length(x0_range),length(x0_EZ_range),N)
syncWidthMatrix = NaN(length(x0_range), length(x0_EZ_range)); 
StabilityMatrix = NaN(length(x0_range), length(x0_EZ_range)); 
TotalCouplingMatrix= NaN(length(x0_range), length(x0_EZ_range),N,N);
MaxNonzeroLambdaMatrix = NaN(length(x0_range), length(x0_EZ_range)); 
MinNonzeroLambdaMatrix = NaN(length(x0_range), length(x0_EZ_range)); 
BeAwayStabilityMatrix= NaN(length(x0_range), length(x0_EZ_range),N) 
PushAwayStabilityMatrix= NaN(length(x0_range), length(x0_EZ_range),N) 


%%-dual loop x0 and x0(EZ)
for ix0 = 1:length(x0_range)
    for ix0_EZ = 1:length(x0_EZ_range)
        x0=x0_range(ix0)+zeros(N,1);
        x0(EZ) = x0_EZ_range(ix0_EZ);
%%*******************************************z_fixed_point******************************************
        %%- solve fixed-point
        Z0 = 3 + zeros(N, 1); 
        one_dim_epileptor_fun = @(z) oneDepileptor(z, x0, K, tau0); 
        opt = optimset('TolFun', 1e-14, 'TolX', 1e-14);
        [z_fixed_numerical, fval, exitflag, output] = fsolve(one_dim_epileptor_fun, Z0, opt); 

        %%- save zi for numerical fixed-point
        NumericFixedPointMatrix(ix0, ix0_EZ) = z_fixed_numerical(EZ);

        %%- Create plot for z_fixed_numerical
        %figure; 
        %plot(1:N, z_fixed_numerical, '-o'); 
        %title(['x0(EZ)=(', num2str(x0(EZ)), ')']);
        %xlabel('Nodes');
        %ylabel('z_fixed_numerical');
        %grid on; % Add grid for better readability

%%**********************************************Coupling Matrix****************************************
        %%- calculate Coupling Matrix for numerical solution
        C = CouplingMatrix(z_fixed_numerical, K, tau0);
        TotalCouplingMatrix(ix0, ix0_EZ,:,:) = C;

        %% - plot out details for coupling matrix
        C_values_EZ_EZ = NaN(length(x0_EZ_range), 1);
        C_values_other_EZ = NaN(length(x0_EZ_range), N-1); 
        C_sum_col_EZ = NaN(length(x0_EZ_range), 1);

	C_values_EZ_other= NaN(length(x0_EZ_range), N-1);
	C_sum_row_EZ = NaN(length(x0_EZ_range), 1);


        for ix0 = 1:length(x0_range)
            for ix0_EZ = 1:length(x0_EZ_range) 
	            C = squeeze(TotalCouplingMatrix(ix0, ix0_EZ, :, :)); 
	            C_values_EZ_EZ(ix0_EZ) = C(EZ, EZ); 
	            C_values_other_EZ(ix0_EZ,:)=C([1:EZ-1,EZ+1:N],EZ)
	            C_sum_col_EZ(ix0_EZ) = sum(C(:, EZ));

		    C_values_EZ_other(ix0_EZ,:) =C(EZ,[1:EZ-1,EZ+1:N])
		    C_sum_row_EZ(ix0_EZ) = sum(C(EZ, :));
			
            end
        end

	C_EZtoEZ_mod=abs(C_values_EZ_EZ)
	C_OthertoEZ_mod=abs(C_values_other_EZ)
	C_sum_EZcol_mod=abs(C_sum_col_EZ)
	
	C_EZtoOther_mod=abs(C_values_EZ_other)
	C_sum_EZrow_mod=abs(C_sum_row_EZ)

        %%- plot subgraph version
        %figure;

        %subplot(3, 1, 1);
        %plot(x0_EZ_range, C_EZtoEZ_mod, 'LineWidth', 2);
        %xlabel('x0(EZ)');
        %ylabel('C(EZ,EZ)');
        %title('Self-Coupling');
        %legend('C(EZ, EZ)');

        %subplot(3, 1, 2);
        %hold on;
        %for node = 1:(N-1)
            %plot(x0_EZ_range, C_OthertoEZ_mod(:, node), 'LineWidth', 1);
        %end
        %hold off;
        %xlabel('x0(EZ)');
        %ylabel('C(other nodes, EZ)');
        %title('Cross-Coupling');
        %legend('C(other nodes, EZ)');

        %subplot(3, 1, 3);
        %plot(x0_EZ_range, C_sum_EZcol_mod, 'LineWidth', 2);
        %xlabel('x0(EZ)');
        %ylabel('sum(C(:, EZ))');
        %title('Sum of Self and Cross Coupling to EZ');
        %legend('sum(C(:, EZ))');


%%******************************************Synchronizability*****************************************
        %%- calcaulate sycnchronizability for numerical solution
        %[sortedLambdasCU, UcovarianceU, err] = covarianceUGaussianNet(C, false, 10000000, false, 1);
        % syncWidth = synchronizability(UcovarianceU); 
                
        %%- save syncWidth into matrix
        %syncWidthMatrix(ix0, ix0_EZ) = syncWidth;
        %fprintf('Synchronization Width: %.3f\\n', syncWidth);

        %%- Create plot for synchronizability
	    %figure;
	    %plot(x0_EZ_range,syncWidthMatrix(ix0,:),'-o');
	    %xlabel('x0(EZ)')
	    %ylabel('Synchronizability')
	    %title('Synchronizability over x0(EZ)')
        %grid on
      
%%*******************************************Stability*************************************************
        %%- calculate Stability
        %[covariance,err] = con2cov(C,false,100000,10,true)
        %stability = Stability(covariance); 
       
        %%- save Stability into matrix
        %StabilityMatrix(ix0, ix0_EZ) = stability;
        %fprintf('Stability: %.3f\\n', stability);

	    %% Create plot for Stability
	    %figure;
	    %plot(x0_EZ_range,StabilityMatrix(ix0,:),'-o');
	    %xlabel('x0(EZ)')
	    %ylabel('Stability')
	    %title('Stability over x0(EZ)')
        %grid on

%%*******************************************Largest Eigenvalue*************************************************
        %%- calculate and save largest eigenvalue
        %eigen=eig(C)
        %eigen_modulus=abs(eigen) % calculate the modulus of the eigenvalues
        %nonzero_modulus = eigen_modulus(eigen_modulus ~= 0); %Exclude eigenvalues with modulus equal to 0
        %% Handle the special case where there are no non-zero eigenvalues
        %if isempty(nonzero_modulus)
            %error('All eigenvalue moduli are zero');
        %end
        %max_nonzero_modulus = max(nonzero_modulus);
        %min_nonzero_modulus = min(nonzero_modulus);
        %% save the result
        %MaxNonzeroLambdaMatrix(ix0, ix0_EZ) = max_nonzero_modulus;
        %MinNonzeroLambdaMatrix(ix0, ix0_EZ) = min_nonzero_modulus;
       
        %% Create plot for Largest Eigenvalue
	    %figure;
	    %plot(x0_EZ_range,MaxNonzeroLambdaMatrix(ix0,:),'-o');
	    %xlabel('x0(EZ)')
	    %ylabel('Largest Eigenvalue')
	    %title('Largest Eigenvalue')
        %grid on

%%*******************************************BeAwayStability*************************************************

        %%- calculate BeAwayStability
        %BeAwayStability = NaN(1, N);
        %trace = diag(covariance);
        %for i = 1:N
            %BeAwayStability(i) = trace(i);
        %end
        %BeAwayStabilityMatrix(ix0, ix0_EZ, :) = BeAwayStability;
        %fprintf('BeAwayStability: %.3f\\n', BeAwayStability(EZ));

        %%- calculate PushAwayStability
        %C_transpose = C';
        %[covariance_new, err] = con2cov(C_transpose, false, 100000, 10, 1);
        %trace_new = diag(covariance_new);
        %PushAwayStability = NaN(1, N);
        %for i = 1:N
            %PushAwayStability(i) = trace_new(i);
        %end
        %PushAwayStabilityMatrix(ix0, ix0_EZ, :) = PushAwayStability;
        %fprintf('PushAwayStability: %.3f\\n', PushAwayStability(EZ));


    end
end



%%- turn off diary
%diary off

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
