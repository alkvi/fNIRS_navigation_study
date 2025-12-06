%% Select dataset

clear;

% File locations
results_file = "../data/results_protocol_2.csv";

% Load
subjstats_file = "../data/mat_files/SubjStats_setup_2_cbsi.mat";
SubjStats = importdata(subjstats_file);

%% Set groups in demographics

demographics = nirs.createDemographicsTable(SubjStats);
disp(demographics);

ya_ids = readtable("../../Park-MOVE_fnirs_dataset_v2/identifiers_YA.csv");
oa_ids = readtable("../../Park-MOVE_fnirs_dataset_v2/identifiers_OA.csv");
pd_ids = readtable("../../Park-MOVE_fnirs_dataset_v2/identifiers_PD.csv");
subj_ids = cell2table(demographics.SubjectID, "VariableNames", ["id_nummer"]);

ya_idx = ismember(subj_ids, ya_ids);
oa_idx = ismember(subj_ids, oa_ids);
pd_idx = ismember(subj_ids, pd_ids);

demographics.group = repmat("NA", height(demographics), 1);
demographics.group(ya_idx) = "YA";
demographics.group(oa_idx) = "OA";
demographics.group(pd_idx) = "PD";

job=nirs.modules.AddDemographics;
job.varToMatch = 'UUID';
job.demoTable=demographics;
SubjStats=job.run(SubjStats);

demographics = nirs.createDemographicsTable(SubjStats);
disp(demographics);

fprintf('Sum YA: %d\n', sum(ya_idx));
fprintf('Sum OA: %d\n', sum(oa_idx));
fprintf('Sum PD: %d\n', sum(pd_idx));

%% Read clusters

oa_cluster = readtable("../data/clusters_oa.csv");
pd_cluster = readtable("../data/clusters_pd.csv");

combined_cluster = [oa_cluster; pd_cluster];

%% Add cluster to demographics

demographics.cluster = repmat("NA", height(demographics), 1);

% Add cluster data
for idx=1:height(combined_cluster)
    subj_id_seek = string(combined_cluster.subject(idx)); 
    match_idx = strcmp(demographics.SubjectID, subj_id_seek);
    if sum(match_idx) < 1
        continue
    end
    demographics.cluster(match_idx) = repmat(num2str(combined_cluster.cluster(idx)), sum(match_idx), 1);
end

job = nirs.modules.AddDemographics;
job.demoTable = demographics;
job.varToMatch = 'UUID';
SubjStats = job.run(SubjStats);

%% Select groups

demographics = nirs.createDemographicsTable(SubjStats);

ya_idx = strcmpi(demographics.group,'YA');
oa_idx = strcmpi(demographics.group,'OA');
pd_idx = strcmpi(demographics.group,'PD');

% Keep all subjects
nan_idx = logical(zeros(height(demographics),1));

% Only YA
selected_idx = zeros(size(SubjStats,1),1);
selected_idx(ya_idx,1) = 1;
selected_idx(nan_idx,1) = 0;
selected_idx = logical(selected_idx);
SubjStats_ya = SubjStats(selected_idx);
demographics_ya = nirs.createDemographicsTable(SubjStats_ya);

% Only OA
selected_idx = zeros(size(SubjStats,1),1);
selected_idx(oa_idx,1) = 1;
selected_idx(nan_idx,1) = 0;
selected_idx = logical(selected_idx);
SubjStats_oa = SubjStats(selected_idx);
demographics_oa = nirs.createDemographicsTable(SubjStats_oa);

% Only PD
selected_idx = zeros(size(SubjStats,1),1);
selected_idx(pd_idx,1) = 1;
selected_idx(nan_idx,1) = 0;
selected_idx = logical(selected_idx);
SubjStats_pd = SubjStats(selected_idx);
demographics_pd = nirs.createDemographicsTable(SubjStats_pd);

fprintf('Sum OA: %d\n', length(SubjStats_oa));
fprintf('Sum PD: %d\n', length(SubjStats_pd));

%% Check age

ages_cluster1 = demographics_oa.age(demographics_oa.cluster == '1');
ages_cluster2 = demographics_oa.age(demographics_oa.cluster == '2');
fprintf('Age OA cluster 1: %.3g\n', mean(ages_cluster1));
fprintf('Age OA cluster 2: %.3g\n', mean(ages_cluster2));

ages_cluster1 = demographics_pd.age(demographics_pd.cluster == '1');
ages_cluster2 = demographics_pd.age(demographics_pd.cluster == '2');
fprintf('Age PD cluster 1: %.3g\n', mean(ages_cluster1));
fprintf('Age PD cluster 2: %.3g\n', mean(ages_cluster2));

%% Set up ROI

% This selection excludes bad quality channels
% in areas where most subjects have data set to NaN
source = [8 6 5 3 2]';    
detector = [NaN NaN NaN NaN NaN]';
ROI_PFC = table(source,detector);

% note: S7-D7 / F4-F6 or S1-D1 / F3-F5 also included even with spec < 30

% ba46
source = [6 8 3 2]';
detector = [7 7 1 1]';
ROI_ba_46 = table(source,detector);
weights_ba46 = [0.47, 0.43, 0.47, 0.43]';
weights_ba46 = normalize(weights_ba46, 'norm', 1);
ROI_ba_46.weight = weights_ba46;

%% Add diagnostics

function [table_out] = AddDiagnostics(roi_result, formula, group, model_name, covar)

    % Diagnostics and validation
    [AIC, BIC] = PlotDiagnostics(roi_result.model{1}, model_name , covar);
    roi_result.formula = repmat(string(formula), size(roi_result,1),1);
    roi_result.AIC = repmat(AIC, size(roi_result,1),1);
    roi_result.BIC = repmat(AIC, size(roi_result,1),1);
    roi_result.group = repmat(group, size(roi_result,1),1);
    table_out = roi_result;

end

%% Cluster - OA

folder = '../figures/figures_p2_oa';
formula = 'beta ~ -1 + cond:cluster + (1|SubjectID)'; 

% Run group model
job = nirs.modules.MixedEffects();
job.formula = formula;
job.dummyCoding = 'full';
job.include_diagnostics = true;
GroupStats_oa = job.run(SubjStats_oa);

GroupStats_oa.probe.defaultdrawfcn='3D mesh (frontal)';
GroupStats_oa.probe = GroupStats_oa.probe.SetFiducials_Visibility(false);
GroupStats_oa.draw('tstat', [-10 10], 'q < 0.05');
GroupStats_oa.printAll('tstat', [-10 10], 'q < 0.05', folder, 'png')

% PFC
roi_result_oa = nirs.util.roiAverage(GroupStats_oa, ROI_PFC, "PFC");
disp(roi_result_oa);
roi_result_oa_pfc = AddDiagnostics(roi_result_oa, formula, "OA", "OA_model", "");

% BA46
roi_result_oa = nirs.util.roiAverage(GroupStats_oa, ROI_ba_46, "BA46");
disp(roi_result_oa);
roi_result_oa_ba46 = AddDiagnostics(roi_result_oa, formula, "OA", "OA_model", "");

c = [-1 1 0 0];
ContrastStats = GroupStats_oa.ttest(c);
ContrastStats.probe.defaultdrawfcn='3D mesh (frontal)';
ContrastStats.probe = ContrastStats.probe.SetFiducials_Visibility(false);
ContrastStats.draw('tstat', [-10 10], 'q < 0.05');
ContrastStats.printAll('tstat', [-10 10], 'q < 0.05', folder, 'png')

roi_result_oa_contrast_1_pfc = nirs.util.roiAverage(ContrastStats, ROI_PFC, "PFC");
roi_result_oa_contrast_1_pfc.group = repmat("OA", size(roi_result_oa_contrast_1_pfc,1),1);
disp(roi_result_oa_contrast_1_pfc);
roi_result_oa_contrast_1_ba46 = nirs.util.roiAverage(ContrastStats, ROI_ba_46, "BA46");
roi_result_oa_contrast_1_ba46.group = repmat("OA", size(roi_result_oa_contrast_1_ba46,1),1);
disp(roi_result_oa_contrast_1_ba46);

c = [0 0 -1 1];
ContrastStats = GroupStats_oa.ttest(c);
ContrastStats.probe.defaultdrawfcn='3D mesh (frontal)';
ContrastStats.probe = ContrastStats.probe.SetFiducials_Visibility(false);
ContrastStats.draw('tstat', [-10 10], 'q < 0.05');
ContrastStats.printAll('tstat', [-10 10], 'q < 0.05', folder, 'png')

roi_result_oa_contrast_2_pfc = nirs.util.roiAverage(ContrastStats, ROI_PFC, "PFC");
roi_result_oa_contrast_2_pfc.group = repmat("OA", size(roi_result_oa_contrast_2_pfc,1),1);
disp(roi_result_oa_contrast_2_pfc);
roi_result_oa_contrast_2_ba46 = nirs.util.roiAverage(ContrastStats, ROI_ba_46, "BA46");
roi_result_oa_contrast_2_ba46.group = repmat("OA", size(roi_result_oa_contrast_2_ba46,1),1);
disp(roi_result_oa_contrast_2_ba46);

%% Cluster - PD

folder = '../figures/figures_p2_pd';

% Run group model
job = nirs.modules.MixedEffects();
job.formula = 'beta ~ -1 + cond:cluster + (1|SubjectID)'; 
job.dummyCoding = 'full';
job.include_diagnostics = true;
GroupStats_pd = job.run(SubjStats_pd);

GroupStats_pd.probe.defaultdrawfcn='3D mesh (frontal)';
GroupStats_pd.probe = GroupStats_pd.probe.SetFiducials_Visibility(false);
GroupStats_pd.draw('tstat', [-10 10], 'q < 0.05');
GroupStats_pd.printAll('tstat', [-10 10], 'q < 0.05', folder, 'png')

% PFC
roi_result_pd = nirs.util.roiAverage(GroupStats_pd, ROI_PFC, "PFC");
disp(roi_result_pd);
roi_result_pd_pfc = AddDiagnostics(roi_result_pd, formula, "PD", "PD_model", "");

% BA46
roi_result_pd = nirs.util.roiAverage(GroupStats_pd, ROI_ba_46, "BA46");
disp(roi_result_pd);
roi_result_pd_ba46 = AddDiagnostics(roi_result_pd, formula, "PD", "PD_model", "");

c = [1 -1 0 0];
ContrastStats = GroupStats_pd.ttest(c);
ContrastStats.probe.defaultdrawfcn='3D mesh (frontal)';
ContrastStats.probe = ContrastStats.probe.SetFiducials_Visibility(false);
ContrastStats.draw('tstat', [-10 10], 'q < 0.05');
ContrastStats.printAll('tstat', [-10 10], 'q < 0.05', folder, 'png')

roi_result_pd_contrast_1_pfc = nirs.util.roiAverage(ContrastStats, ROI_PFC, "PFC");
roi_result_pd_contrast_1_pfc.group = repmat("PD", size(roi_result_pd_contrast_1_pfc,1),1);
disp(roi_result_pd_contrast_1_pfc);
roi_result_pd_contrast_1_ba46 = nirs.util.roiAverage(ContrastStats, ROI_ba_46, "BA46");
roi_result_pd_contrast_1_ba46.group = repmat("PD", size(roi_result_pd_contrast_1_ba46,1),1);
disp(roi_result_pd_contrast_1_ba46);

c = [0 0 1 -1];
ContrastStats = GroupStats_pd.ttest(c);
ContrastStats.probe.defaultdrawfcn='3D mesh (frontal)';
ContrastStats.probe = ContrastStats.probe.SetFiducials_Visibility(false);
ContrastStats.draw('tstat', [-10 10], 'q < 0.05');
ContrastStats.printAll('tstat', [-10 10], 'q < 0.05', folder, 'png')

roi_result_pd_contrast_2_pfc = nirs.util.roiAverage(ContrastStats, ROI_PFC, "PFC");
roi_result_pd_contrast_2_pfc.group = repmat("PD", size(roi_result_pd_contrast_2_pfc,1),1);
disp(roi_result_pd_contrast_2_pfc);
roi_result_pd_contrast_2_ba46 = nirs.util.roiAverage(ContrastStats, ROI_ba_46, "BA46");
roi_result_pd_contrast_2_ba46.group = repmat("PD", size(roi_result_pd_contrast_2_ba46,1),1);
disp(roi_result_pd_contrast_2_ba46);

%% Save result tables

% Compile ROI results
result_table = [ 
    roi_result_oa_pfc;
    roi_result_oa_ba46;
    roi_result_pd_pfc;
    roi_result_pd_ba46;
];

% Drop unneeded columns
result_table(:,[10,11]) = [];
disp(result_table);

% Compile contrast results
contrast_table = [
    roi_result_oa_contrast_1_pfc;
    roi_result_oa_contrast_1_ba46;
    roi_result_pd_contrast_1_pfc;
    roi_result_pd_contrast_1_ba46;
    roi_result_oa_contrast_2_pfc;
    roi_result_oa_contrast_2_ba46;
    roi_result_pd_contrast_2_pfc;
    roi_result_pd_contrast_2_ba46;
];
contrast_table(:,[10]) = [];
disp(contrast_table);

% Write table 
writetable(result_table, results_file)
contrast_file = strrep(results_file, '.csv', '_contrast.csv');
writetable(contrast_table, contrast_file)

%% Check assumptions

function [AIC, BIC] = PlotDiagnostics(model, titlestr, covar)

fitted_vals = model.Fitted;
residuals_raw = model.Residuals.Raw;
residuals_standardized = model.Residuals.Standardized;
leverage = model.Diagnostics.Leverage;

% residuals are (approximately) normally distributed (Q-Q plot)
figure;
subplot(3,2,1);
histfit(residuals_raw);
title('Histogram - raw residuals');
subplot(3,2,2);
qqplot(residuals_raw);
title('Q-Q - raw residuals');

% Residual vs fitted to assess possible heteroscedascity
% Plot fitted vs residuals
subplot(3,2,3);
h = plot(fitted_vals, residuals_standardized, 'bx');
xlabel("fitted");
ylabel("residual (standardized)");
hold on;

% Add 0 line and a polyline
p = polyfit(fitted_vals, residuals_standardized, 4);
px = linspace(min(fitted_vals), max(fitted_vals));
py = polyval(p, px);
plot(px, py, 'r-', 'LineWidth', 1);
xlim = [min(fitted_vals) max(fitted_vals)];
line(xlim,[0 0],'LineStyle','--', 'LineWidth', 1);
v = axis; % get current values
lo = min( v(1:2:end) ); % lower limit
up = max( v(2:2:end) ); % uppper limit
axis( [lo up lo up] );
hold off;
title('Fitted-Residual');

% Leverage vs residual
subplot(3,2,4);
plot(leverage, residuals_standardized, 'bx');
hold on;
xlabel("leverage");
ylabel("residual (standardized)");
p = polyfit(leverage, residuals_standardized, 4);
px = linspace(min(leverage), max(leverage));
py = polyval(p, px);
plot(px, py, 'r-', 'LineWidth', 1);
xlim = [min(leverage) max(leverage)];
line(xlim,[0 0],'LineStyle','--', 'LineWidth', 1);
hold off;
title('Leverage-Residual');

% Linearity in covariate
if covar ~= ""
    subplot(3,2,5);
    plot(residuals_raw, table2array(model.Variables(:,covar)), 'bx');
    title('Residual-Covariate');
end

sgtitle(titlestr) 

% Save the figure
fig = gcf;
set(gcf,'Position',[100 100 1000 700])
figname  = titlestr + "_" + covar + ".png";
figname = strrep(figname, " ", "_");
%exportgraphics(fig,figname,'Resolution',300);

AIC = model.ModelCriterion.AIC;
BIC = model.ModelCriterion.BIC;

end