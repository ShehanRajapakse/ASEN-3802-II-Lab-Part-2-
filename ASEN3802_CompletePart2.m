%% ASEN 3802 | Lab 1 Part 2

clc; clear; close all;

% read in files
files = { ...
    'Case 1 - load in center.csv', ... % center load
    'Case 2 data.txt', ...             % single off-center load
    'Case 3 data.txt'  ...             % two symmetric loads
};

caseNames = {'Case 1: Center Load', 'Case 2: Single Off-Center Load', 'Case 3: Two Symmetric Loads'};

% constants
lbf2N = 4.44822;
in2m  = 0.0254;

% load the data
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

%% Part 1

% Plot 1 - Reaction Forces

% F0
figure(1);
grid on;
for k  = 1:length(Cases)
    % Scatter
    hold on;
    scatter(Cases(k).P, Cases(k).F0, 'filled', 'Marker', 'o', 'SizeData', 40);
    
    % Calculate line of best fit
    p0 = polyfit(Cases(k).P, Cases(k).F0, 1);
    fit_F0 = polyval(p0, Cases(k).P);

    % Plotting Best fit
    plot(Cases(k).P, fit_F0,'-');

end
xlabel('Applied Load [P] (N)');
ylabel('F0 | Reaction Forces (N)');
title('F0 Reaction Force v.s Applied Load');
legend('Case 1 Data', 'Case 1: Line of Best Fit', 'Case 2 Data', 'Case 2: Line of Best Fit', 'Case 3 Data', 'Case 3: Line of Best Fit', 'Location','bestoutside');

% F1
figure(2);
grid on;
for k  = 1:length(Cases)
    % Data Plot
    hold on;
    scatter(Cases(k).P, Cases(k).F1, 'filled', 'Marker', 'd', 'SizeData', 40);

    % Calculate line of best fit
    p1 = polyfit(Cases(k).P, Cases(k).F1, 1);
    fit_F1 = polyval(p1, Cases(k).P);

    % Plotting Best fit
    plot(Cases(k).P, fit_F1,'-');

end
xlabel('Applied Load [P] (N)');
ylabel('F1 | Reaction Forces (N)');
title('F1 Reaction Force v.s Applied Load');
hold off;
legend('Case 1 Data', 'Case 1: Line of Best Fit', 'Case 2 Data', 'Case 2: Line of Best Fit', 'Case 3 Data', 'Case 3: Line of Best Fit', 'Location','bestoutside');

% F2
figure(3);
grid on;
for k  = 1:length(Cases)
    % Data Plot
    hold on;
    scatter(Cases(k).P, Cases(k).F2, 'filled', 'Marker', 'h', 'SizeData', 80);

    % Calculate line of best fit
    p2 = polyfit(Cases(k).P, Cases(k).F2, 1);
    fit_F2 = polyval(p2, Cases(k).P);

    % Plot For Best fit
    plot(Cases(k).P, fit_F2,'-');
end
xlabel('Applied Load [P] (N)');
ylabel('F2 | Reaction Forces (N)');
title('F2 Reaction Forces v.s Applied Load');
legend('Case 1 Data', 'Case 1: Line of Best Fit', 'Case 2 Data', 'Case 2: Line of Best Fit', 'Case 3 Data', 'Case 3: Line of Best Fit', 'Location','bestoutside');
hold off;

% F3D
figure(4);
grid on;
for k  = 1:length(Cases)
    % Data Plot
    hold on;
    scatter(Cases(k).P, Cases(k).Fint, 'filled', 'Marker', '^', 'SizeData', 50);

    % Calculate line of best fit
    p3D = polyfit(Cases(k).P, Cases(k).Fint, 1);
    fit_F3D = polyval(p3D, Cases(k).P);

    % Plot For Best fit
    plot(Cases(k).P, fit_F3D,'-');
end
xlabel('Applied Load [P] (N)');
ylabel('F3D | Reaction Force (N)');
title('F3D Reaction Forces v.s Applied Load');
legend('Case 1 Data', 'Case 1: Line of Best Fit', 'Case 2 Data', 'Case 2: Line of Best Fit', 'Case 3 Data', 'Case 3: Line of Best Fit', 'Location','bestoutside');
hold off;

% LVDT
figure(5);
grid on;
for k  = 1:length(Cases)
    % Data Plot
    hold on;
    x = Cases(k).P;
    y = Cases(k).LVDT;
    
    scatter(x, y, 'filled', 'Marker', 's', 'SizeData', 50);

    % Remove any NaNs just in case
    mask = isfinite(x) & isfinite(y);
    x = x(mask);
    y = y(mask);

    % Calculate line of best fit
    pLVDT = polyfit(x, y, 1);
    xfit = linspace(min(Cases(k).P), max(Cases(k).LVDT), 1);
    yfit = polyval(pLVDT, xfit);

    % Plot For Best fit
    plot(xfit, yfit,'-', 'LineWidth', 2);
end
xlabel('Applied Load [P] (N)');
ylabel('LVDT | Displacement (m)');
title('LVDT Displacement v.s Applied Load');
legend('Case 1 Data', 'Case 1: Line of Best Fit', 'Case 2 Data', 'Case 2: Line of Best Fit', 'Case 3 Data', 'Case 3: Line of Best Fit', 'Location','bestoutside');

% linear regression: LVDT cs P (all Cases)
for k = 1:length(Cases)

    % Example: LVDT regression
    x = Cases(k).P;
    y = Cases(k).LVDT;

    % Remove any NaNs just in case
    mask = isfinite(x) & isfinite(y);
    x = x(mask);
    y = y(mask);

    %linear regression: y = a*x + b 
    scatter(x,y, 50, 'filled')
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

end
xlabel('Applied Load P (N)');
ylabel('Displacement (m)');
title('LVDT Displacement vs Load with Regression Fits');
legend('Case 1 Data','Case 1 Fit', 'Case 2 Data','Case 2 Fit', 'Case 3 Data','Case 3 Fit', 'Location','bestoutside');

% Regression Plot: F2 vs P (All Cases)

    %Plot regression 
    figure(7);
    hold on;
    grid on;
        for k = 1:length(Cases)
            x = Cases(k).P;
            y = Cases(k).F2;
            mask = isfinite(x) & isfinite(y);
            x = x(mask);
            y = y(mask);
            scatter(x,y, 50, 'Filled');
            p = polyfit(min(x), max(x), 100);
            %% 
            xfit = linspace(min(x), max(x), 100);
            yfit = polyval(p, xfit);
            plot(xfit, yfit, 'LineWidth', 1.5);
        end
    xlabel('Applied Load P (N)');
    ylabel('Reaction Force F2 (N)');
    title('F2 Reaction Force vs Load with Regression Fits');
    legend('Case 1 Data','Case 1 Fit', 'Case 2 Data','Case 2 Fit', 'Case 3 Data','Case 3 Fit', 'Location','bestoutside');
        
%     % Residual plot (linearity check)
%     figure;
%     hold on;
%     grid on;
%     scatter(x, residuals, 'filled');
%     yline(0,'-');
%     xlabel('Applied Load P (N)');
%     ylabel('Residual (m)');
%     title(['Residuals vs Load: ' caseNames{k}]);
% 
%     %Print stats for your discussion
%     fprintf('\n%s\n', caseNames{k});
%     fprintf('Slope = %.4e m/N\n', a);
%     fprintf('Intercept = %.4e m\n', b);
%     fprintf('R^2 = %.4f\n', R2);
%     fprintf('RMSE = %.4e m\n', RMSE);
% end
% 
% 
% % Prepare for analysis of internal forces
% for k = 1:length(Cases)
%     x = Cases(k).P;
%     y = Cases(k).Fint;
% 
%     % Remove any NaNs just in case
%     mask = isfinite(x) & isfinite(y);
%     x = x(mask);
%     y = y(mask);
% 
%     % Linear regression: y = a*x + b 
%     p = polyfit(x, y, 1);
%     a = p(1);    % slope
%     b = p(2);    % intercept
% 
%     % Regression curve for plotting
%     xfit = linspace(min(x), max(x), 100);
%     yfit = polyval(p, xfit);
% 
%     % Predicted values at data points + residuals
%     yhat = polyval(p, x);
%     residuals = y - yhat;
% 
%     % Linearity metrics 
%     SS_res = sum(residuals.^2);
%     SS_tot = sum((y - mean(y)).^2);
%     R2 = 1 - SS_res/SS_tot;
% 
%     RMSE = sqrt(mean(residuals.^2));
% 
%     % Plot regression 
%     figure; hold on; grid on;
%     scatter(x, y, 'filled');
%     plot(xfit, yfit, 'LineWidth', 1.5);
%     xlabel('Applied Load P (N)');
%     ylabel('Internal Force (N)');
%     title(['Linear Regression (polyfit): ' caseNames{k}]);
%     legend('Data','Linear Fit','Location','best');
% 
%     % Residual plot (linearity check)
%     figure; hold on; grid on;
%     scatter(x, residuals, 'filled');
%     yline(0,'-');
%     xlabel('Applied Load P (N)');
%     ylabel('Residual (N)');
%     title(['Residuals vs Load: ' caseNames{k}]);
% 
%     % Print stats for your discussion
%     fprintf('\n%s\n', caseNames{k});
%     fprintf('Slope = %.4e N/N\n', a);
%     fprintf('Intercept = %.4e N\n', b);
%     fprintf('R^2 = %.4f\n', R2);
%     fprintf('RMSE = %.4e N\n', RMSE);
% end