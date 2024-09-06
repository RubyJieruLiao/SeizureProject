

%% define start and end index as parameters
startIndex = 1;
endIndex = 2500;

% Adjust the size of stabilityValues to match the range from startIndex to endIndex
stabilityValues = zeros(1, endIndex - startIndex + 1); % Array to store stability values、

%% it's for the whole network
for kIndex = startIndex:endIndex
    k = motifLengthsToCheck(kIndex);
    % verbose = true to display iteration information
    [covarianceApprox, err] = con2cov(C, false, k, 10, true); % Set verbose to true
    if err == 0
        stability = Stability(covarianceApprox);
        stabilityValues(kIndex - startIndex + 1) = stability; % Store stability value
    else
        disp(['Error occurred for k = ', num2str(k)]);
    end
end