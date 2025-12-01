% Subject-level analysis
% Ran on R2024a, forked version of NIRS toolbox
% https://github.com/alkvi/nirs-toolbox-fork/tree/phd_study_4
%% Load

original_dir = pwd();
cd ../../Park-MOVE_fnirs_dataset_v2/fNIRS_data;

% Load the NIRx probe used
probe_folder = "nirx_probe";
nirx_probe = nirs.io.loadNIRxProbe(probe_folder,true);

% Load BIDS
my_data_dir = 'bids_dataset_snirf';
raw_data = nirs.bids.loadBIDS(my_data_dir, true, true, nirx_probe);

% or load from file
%load('../data/mat_files/raw_data.mat');

% Go back to start
cd(original_dir);

demographics = nirs.createDemographicsTable(raw_data);
demographics.SubjectID = erase(demographics.SubjectID, "sub-");
disp(demographics);

%% Trim baseline

% Trim baseline
job = nirs.modules.TrimBaseline();
job.preBaseline  = 1;
job.postBaseline = 1;
raw_data = job.run(raw_data);

%% Preprocess 

job = nirs.modules.OpticalDensity();
od = job.run(raw_data);

%% QT-NIRS version

% Parameters to reflect default QT-NIRS params
fcut = [0.5 2.5];
window = 5;

% Parameters for data
n_channels = 33;
fs = raw_data(1).Fs;

psp_sci_table = table();
for i=1:size(od,1)

    % Get subject info
    subject_id = od(i).demographics.SubjectID;
    session = od(i).demographics.session;

    % Get SCI and PSP
    result_table = nirs.util.scalp_coupling_index(od(i), fcut, false, 5);
    result_table.subject = repmat({subject_id}, n_channels, 1);
    result_table.session = repmat({session}, n_channels, 1);

    % Append to main table
    psp_sci_table = [psp_sci_table; result_table];

end

writetable(psp_sci_table, "../data/sci_psp_nirs_toolbox.csv");

