

%% clear parameter
clear
clc

%% set basic parameter
tau0 = 6667; 
patient = 'P1';
EZ = 48;
connectivity_path = strcat('../data/connectivity_', patient, '/');
K = load(strcat(connectivity_path, 'weights.txt')); % Load connectivity matrix
K = normal(K); % Normalize
[N, ~] = size(K); % Get the number of nodes


%%%%%%%%%%%%%%%% this is for numerical fixed-point range %%%%%%%%%%%%%%%% \

%%- initialize
format long
x0_range1 = -2.30;
%x0_EZ_range1 = linspace(-1.672548, -1.672545, 10);
x0_EZ_range1 = -2.3:0.1:-1
disp(x0_EZ_range1) 
StabilityMatrix1 = NaN(length(x0_range1), length(x0_EZ_range1)); 

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

     
        %%- calculate Coupling Matrix for numerical solution
        C1 = CouplingMatrix(z_fixed_numerical, K, tau0);
        
        
        %%- calculate Stability
        [covariance1,err] = con2cov(C1,false,100000,10,1)
        stability1 = Stability(covariance1); 
       
        %%- save Stability into matrix
        StabilityMatrix1(ix0, ix0_EZ) = stability1;
        fprintf('Stability1: %.3f\\n', stability1);


      

    end
end






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
