clc; clear; close all;

% read in files
files = { ...
    'Case 1 - load in center.csv', ... % center load
    'Case 2 data.txt', ...             % single off-center load
    'Case 3 data.txt'  ...             % two symmetric loads
};

caseNames = { ...
    'Case 1: Center Load', ...
    'Case 2: Single Off-Center Load', ...
    'Case 3: Two Symmetric Loads' ...
};

%% constants
lbf2N = 4.44822;
in2m  = 0.0254;

%% load the data
Cases = struct([]);
for k = 1:length(files)

    % Robust read for csv/txt (handles headers better than plain readmatrix)
    try
        opts = detectImportOptions(files{k}, 'FileType', 'text');
        D = readmatrix(files{k}, opts);
    catch
        D = readmatrix(files{k});
    end

    Cases(k).P    = D(:,1) * lbf2N;     % Applied load [N]
    Cases(k).F0   = D(:,2) * lbf2N;     % Reaction forces [N]
    Cases(k).F1   = D(:,3) * lbf2N;
    Cases(k).F2   = D(:,4) * lbf2N;
    Cases(k).Fint = D(:,5) * lbf2N;     % Internal force [N]
    Cases(k).LVDT = D(:,6) * in2m;      % Displacement [m]
end

%% Plot 1 - Reaction Forces
figure; hold on; grid on;
for k = 1:length(Cases)
    scatter(Cases(k).P, Cases(k).F0, 'filled');
    scatter(Cases(k).P, Cases(k).F1, 'filled');
    scatter(Cases(k).P, Cases(k).F2, 'filled');
end
xlabel('Applied Load P (N)');
ylabel('Reaction Forces (N)');
title('Reaction Forces vs Applied Load');
legend( ...
    'C1 F0','C1 F1','C1 F2', ...
    'C2 F0','C2 F1','C2 F2', ...
    'C3 F0','C3 F1','C3 F2', ...
    'Location','bestoutside');

%% Plot 2 - Displacement
figure; hold on; grid on;
for k = 1:length(Cases)
    scatter(Cases(k).P, Cases(k).LVDT, 'filled');
end
xlabel('Applied Load P (N)');
ylabel('Displacement (m)');
title('LVDT Displacement vs Applied Load');
legend(caseNames,'Location','best');

%% Plot 3 - Internal Force
figure; hold on; grid on;
for k = 1:length(Cases)
    scatter(Cases(k).P, Cases(k).Fint, 'filled');
end
xlabel('Applied Load P (N)');
ylabel('Internal Force (N)');
title('Internal Force vs Applied Load');
legend(caseNames,'Location','best');

%% linear regression using polyval and polyfit functions
for k = 1:length(Cases)

    % Example: LVDT regression
    x = Cases(k).P;
    y = Cases(k).LVDT;

    % Remove any NaNs just in case
    mask = isfinite(x) & isfinite(y);
    x = x(mask);
    y = y(mask);

    %linear regression: y = a*x + b 
    p = polyfit(x, y, 1);
    a = p(1);    % slope
    b = p(2);    % intercept

    % Regression curve for plotting
    xfit = linspace(min(x), max(x), 100);
    yfit = polyval(p, xfit);

    % Predicted values at data points + residuals
    yhat = polyval(p, x);
    residuals = y - yhat;

    %Linearity metrics 
    SS_res = sum(residuals.^2);
    SS_tot = sum((y - mean(y)).^2);
    R2 = 1 - SS_res/SS_tot;

    RMSE = sqrt(mean(residuals.^2));

    %Plot regression 
    figure; hold on; grid on;
    scatter(x, y, 'filled');
    plot(xfit, yfit, 'LineWidth', 1.5);
    xlabel('Applied Load P (N)');
    ylabel('Displacement (m)');
    title(['Linear Regression (polyfit): ' caseNames{k}]);
    legend('Data','Linear Fit','Location','best');

    % Residual plot (linearity check)
    figure; hold on; grid on;
    scatter(x, residuals, 'filled');
    yline(0,'-');
    xlabel('Applied Load P (N)');
    ylabel('Residual (m)');
    title(['Residuals vs Load: ' caseNames{k}]);

    %Print stats for your discussion
    fprintf('\n%s\n', caseNames{k});
    fprintf('Slope = %.4e m/N\n', a);
    fprintf('Intercept = %.4e m\n', b);
    fprintf('R^2 = %.4f\n', R2);
    fprintf('RMSE = %.4e m\n', RMSE);
end
