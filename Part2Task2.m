%% ASEN 3802 Lab 1 (Spring 2026) - Part 2 Task 2 (Case 2, OFF-CENTER LOAD)
% Script version 
% Reads a TXT file, estimates load position from reactions, and compares
% measured LVDT + F3D to equivalent beam model.
%

clear; clc; close all;

CASE2_TXT = "Case 2 data.txt";   % <-- TXT file path

% Geometry
L   = 4.0;        % [m] span (16 bays * 0.25 m)
bay = 0.250;      % [m] node spacing

% Section properties from Part 1
I_mm4 = 2.47e6;   % [mm^4]
A_mm2 = 39.57;    % [mm^2]
c_mm  = 125.0;    % [mm]

% Material
E = 69e9;         % [Pa] Aluminum 

% Unit conversions 
LB2N = 4.4482216152605;
IN2M = 0.0254;

outDir = "task2_outputs";
if ~exist(outDir, "dir"); mkdir(outDir); end

% Load the TXT file
opts = detectImportOptions(CASE2_TXT, "FileType","text");
% If your TXT is whitespace-delimited, make sure delimiter is whitespace:
opts.Delimiter = {'\t',' ',';','|',','};   % robust; will auto-handle most
T = readtable(CASE2_TXT, opts);
% Assign columns by position (TXT has no headers)
P_lbf   = T.Var1;   % Applied load
F0      = T.Var2;   % Reaction load cell 0
F1      = T.Var3;   % Reaction load cell 1
F2      = T.Var4;   % Reaction load cell 2
F3D     = T.Var5;   % Inline force
LVDT_in = T.Var6;   % Midspan deflection


% Reactions and regression (position estimate)
RA = F0 + F1;
RB = F2;

% Filter out near-zero loads for stable regression
mask = abs(P_lbf) > 0.5;
Pfit  = P_lbf(mask);
RAfit = RA(mask);
RBfit = RB(mask);

% Linear regression y = m*P + b
pA = polyfit(Pfit, RAfit, 1);  mA = pA(1); bA = pA(2);
pB = polyfit(Pfit, RBfit, 1);  mB = pB(1); bB = pB(2);

% Statics: RB/P ≈ x/L. Using slopes is robust to offsets.
% Use both slopes so scale issues cancel:
x_est  = (mB / (mA + mB)) * L;

% Snap to nearest node (loads only at nodes)
nodeIdx = round(x_est / bay);
x_node  = nodeIdx * bay;

fprintf("\n--- TASK 2 (CASE 2) ---\n");
fprintf("Fit R_A = mA*P + bA: mA=%.6g, bA=%.6g\n", mA, bA);
fprintf("Fit R_B = mB*P + bB: mB=%.6g, bB=%.6g\n", mB, bB);
fprintf("Estimated load location: x_est = %.6f m\n", x_est);
fprintf("Nearest node: x_node = %.6f m (node #%d, bay=%.3f m)\n\n", x_node, nodeIdx, bay);

% Equivalent beam model @ midspan for point load at x_node
I = I_mm4 * 1e-12;   % [m^4]
A = A_mm2 * 1e-6;    % [m^2]
c = c_mm  * 1e-3;    % [m]

P_N = P_lbf * LB2N;
x_mid = L/2;

% Deflection prediction at midspan
v_pred_m  = zeros(size(P_N));
for i = 1:numel(P_N)
    v_pred_m(i) = pointLoadDeflection_midspan(x_node, P_N(i), L, E, I);
end
v_pred_in = v_pred_m / IN2M;

% Inline force prediction using bending-stress equivalence at midspan:
% sigma = M*c/I, Fi = sigma*A
Fi_pred_N = zeros(size(P_N));
for i = 1:numel(P_N)
    Mmid = midspanMoment_pointLoad(x_node, P_N(i), L);
    sigma = Mmid * c / I;
    Fi_pred_N(i) = sigma * A;
end

% Preload correction for F3D (subtract value at nearest-zero-load row)
[~, idx0] = min(abs(P_lbf));
F3D0 = F3D(idx0);
Fi_meas = F3D - F3D0;

fprintf("F3D preload correction: using row %d as baseline -> F3D0 = %.6g\n\n", idx0, F3D0);

% PLOTS
Pline = linspace(min(P_lbf), max(P_lbf), 200);

% 1) Reactions vs P + fits
figure("Color","w"); hold on; grid on;
plot(P_lbf, RA, "o", "DisplayName","R_A = F0+F1");
plot(P_lbf, RB, "o", "DisplayName","R_B = F2");
plot(Pline, polyval(pA, Pline), "-", "DisplayName", sprintf("Fit R_A (m=%.3g)", mA));
plot(Pline, polyval(pB, Pline), "-", "DisplayName", sprintf("Fit R_B (m=%.3g)", mB));
xlabel("Applied Load P (lbf)");
ylabel("Reaction (file units)");
title(sprintf("Task 2 Case 2: Reactions vs P | x_{est}=%.3f m, nearest node=%.3f m", x_est, x_node));
legend("Location","best");
saveas(gcf, fullfile(outDir, "task2_reactions_txt.png"));

% 2) LVDT vs P + linear fit + beam prediction
pL = polyfit(Pfit, LVDT_in(mask), 1); mL = pL(1); bL = pL(2);

figure("Color","w"); hold on; grid on;
plot(P_lbf, LVDT_in, "o", "DisplayName","Measured LVDT (in)");
plot(Pline, polyval(pL, Pline), "-", "DisplayName", sprintf("Fit LVDT (m=%.3g)", mL));
plot(P_lbf, v_pred_in, "s-", "DisplayName","Beam model @ midspan (in)");
xlabel("Applied Load P (lbf)");
ylabel("Deflection (in)");
title("Task 2 Case 2: Midspan Deflection (Measured vs Beam Model)");
legend("Location","best");
saveas(gcf, fullfile(outDir, "task2_deflection_txt.png"));

% 3) Inline force trend (measured preload-corrected) vs prediction
figure("Color","w"); hold on; grid on;
plot(P_lbf, Fi_meas, "o", "DisplayName","Measured F3D (preload removed)");
plot(P_lbf, Fi_pred_N, "s-", "DisplayName","Beam-based Fi prediction (N)");
xlabel("Applied Load P (lbf)");
ylabel("Inline force (meas units) / Prediction (N)");
title("Task 2 Case 2: Inline Force Trend Comparison");
legend("Location","best");
saveas(gcf, fullfile(outDir, "task2_inline_force_txt.png"));

fprintf("Saved plots to: %s\n", outDir);

% Local helper functions
function colName = findColumn(vars, varsLower, candidates)
    if ~isstring(candidates); candidates = string(candidates); end
    cands = lower(regexprep(candidates, '\s+', ''));
    for k = 1:numel(cands)
        cand = cands(k);
        idx = find(varsLower == cand, 1);
        if ~isempty(idx)
            colName = vars(idx);
            return;
        end
        idx = find(contains(varsLower, cand), 1);
        if ~isempty(idx)
            colName = vars(idx);
            return;
        end
    end
    error("Could not find column for: %s. Available: %s", ...
        strjoin(candidates, ", "), strjoin(vars, ", "));
end

function v = pointLoadDeflection_midspan(a, P, L, E, I)
% Midspan deflection (x=L/2) of simply supported beam with point load at x=a
% Returns v in meters (negative downward).
    x = L/2;
    b = L - a;
    if x <= a
        v = -P * b * x / (6*L*E*I) * (L^2 - b^2 - x^2);
    else
        v = -P * a * (L - x) / (6*L*E*I) * (L^2 - a^2 - (L - x)^2);
    end
end

function Mmid = midspanMoment_pointLoad(a, P, L)
% Bending moment at midspan x=L/2 for point load at x=a
    x = L/2;
    RA = P * (L - a) / L;
    if a <= x
        Mmid = RA*x - P*(x - a);
    else
        Mmid = RA*x;
    end
end
