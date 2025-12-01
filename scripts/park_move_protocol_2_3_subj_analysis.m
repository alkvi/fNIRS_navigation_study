% Subject-level analysis
% Ran on R2024a, forked version of NIRS toolbox
% https://github.com/alkvi/nirs-toolbox-fork/tree/phd_study_4
%% Load

original_dir = pwd();
cd 'C:\Users\alkvis\OneDrive - Karolinska Institutet\Dokument\Project\fNIRS_projects\Park-MOVE_fnirs_dataset_v2\fNIRS_data';

% Load the NIRx probe used
probe_folder = "nirx_probe";
nirx_probe = nirs.io.loadNIRxProbe(probe_folder,true);

% Load BIDS
my_data_dir = 'bids_dataset_snirf';
raw_data = nirs.bids.loadBIDS(my_data_dir, true, true, nirx_probe);

% Go back to start
cd(original_dir);

% Also save the data to file
%save('data/mat_files/raw_data.mat','raw_data');

% ..or load
%load('data/mat_files/raw_data.mat');

demographics = nirs.createDemographicsTable(raw_data);
demographics.SubjectID = erase(demographics.SubjectID, "sub-");
disp(demographics);

%% Add age to demographics

age_data = readtable('../Park-MOVE_fnirs_dataset_v2/basic_demographics.csv'); 
demographics.age = NaN(height(demographics),1);

% Add age data
for idx=1:height(age_data)
    subj_id_seek = string(age_data.subject(idx)); 
    match_idx = strcmp(demographics.SubjectID, subj_id_seek);
    if sum(match_idx) < 1
        continue
    end
    demographics(match_idx,:).age = ones(sum(match_idx),1) * age_data(idx,:).age;
end

job = nirs.modules.AddDemographics;
job.demoTable = demographics;
job.varToMatch = 'UUID';
raw_data = job.run(raw_data);

%% Trim baseline

% Trim baseline
job = nirs.modules.TrimBaseline();
job.preBaseline  = 1;
job.postBaseline = 1;
raw_data = job.run(raw_data);

%% Prune bad channels in data

quality_file = "../data/sci_psp_nirs_toolbox.csv";
quality_table = readtable(quality_file);

% Combined criteria on SCI and PSP
quality_table_filt = quality_table(quality_table.sci < 0.1 | quality_table.power < 0.1, :);

% Grab the long channels
quality_table_filt = quality_table_filt(quality_table_filt.detector < 7, :);

% Loop through bad channel table for each src-det pair.
for i = 1:size(quality_table_filt,1)
    quality_row = quality_table_filt(i,:);
    
    % Find which subjects this row pertains to.
    subject_seek = string(quality_row.subject);
    protocol_seek = string(quality_row.session);
    
    % raw_data(j).data is TxN where N is number of channels
    % replace raw_data(j).data with NaN in col N for bad channel in N
    subj_idx = find(strcmp(demographics.SubjectID, subject_seek) & strcmp(demographics.session, protocol_seek));

    idx_to_null = raw_data(subj_idx).probe.link.source == quality_row.source ...
        & raw_data(subj_idx).probe.link.detector == quality_row.detector;
    
    fprintf("Nulling src-det %d-%d for subj %s\n", quality_row.source, quality_row.detector, subject_seek);
    raw_data(subj_idx).data(:,idx_to_null) = NaN;
end

%% Preprocess 

% Label short-sep
job = nirs.modules.LabelShortSeperation();
job.max_distance = 8;
raw_data = job.run(raw_data);

% Extra short-sep labelling:
% Some short-detectors (virtual detector 8-15) make up links 
% with not just the source they're attached to, but ones further away.
% these are still short-separation (so e.g. 1-14 and 1-15 detect the same
% signal). Mark these also. NOTE: the column is called "ShortSeperation" (sic).
for i = 1:numel(raw_data)
    short_separation_idx =(raw_data(i).probe.link.detector >= 8);
    raw_data(i).probe.link.ShortSeperation(short_separation_idx) = 1;
end

% Set each duration to 20 seconds, make sure amp is 1.
raw_data = nirs.design.change_stimulus_duration(raw_data,[],20);
for subj_idx=1:length(raw_data)
    for key_idx=1:length(raw_data(subj_idx).stimulus.keys)
        key = raw_data(subj_idx).stimulus.keys{key_idx};
        raw_data(subj_idx).stimulus(key).amp(:) = 1;
    end
end

job = nirs.modules.OpticalDensity();
od = job.run(raw_data);

job = nirs.modules.TDDR();
od = job.run(od);

job = nirs.modules.BeerLambertLaw();
hb = job.run(od);

%% Visualize

nirs.viz.TimeSeriesViewer(hb);  

demographics = nirs.createDemographicsTable(hb);
disp(demographics);

figure;
raw_data(1).probe.defaultdrawfcn='3D mesh';
raw_data(1).probe.draw;

%% Get combined measure

% Add CBSI
job = nirs.modules.CalculateCBSI();
hb = job.run(hb);

job=nirs.modules.KeepTypes;
job.types={'cbsi', 'hbo', 'hbr'};
hb=job.run(hb);

%% Get data for protocol 2

p2_idx = strcmp(demographics.session, "protocol2");
hb_setup2 = hb(p2_idx);
demographics_setup2 = nirs.createDemographicsTable(hb_setup2);

% Skip idx 57 who did not follow protocol
hb_setup2([57]) = [];

% Now we need to trim blocks where participants started off wrong 
% or there were issues at the end
demographics_setup2 = nirs.createDemographicsTable(hb_setup2);
stim_table_setup2 = nirs.createStimulusTable(hb_setup2);

% Cut 1044 to last ST_walk (already missing last ST_nav) - 42
stim_table_setup2{42,3}.onset(6) = [];
stim_table_setup2{42,3}.dur(6) = [];
stim_table_setup2{42,3}.amp(6) = [];

% Cut 1063 to ST_nav 1 - 60
stim_table_setup2{60,3}.onset(1) = [];
stim_table_setup2{60,3}.dur(1) = [];
stim_table_setup2{60,3}.amp(1) = [];

% Cut 1064 to ST_walk 2 - 61
stim_table_setup2{61,2}.onset(1) = [];
stim_table_setup2{61,2}.dur(1) = [];
stim_table_setup2{61,2}.amp(1) = [];
stim_table_setup2{61,3}.onset(1) = [];
stim_table_setup2{61,3}.dur(1) = [];
stim_table_setup2{61,3}.amp(1) = [];

% Cut 2008 to ST_nav 1 - 97
stim_table_setup2{97,3}.onset(1) = [];
stim_table_setup2{97,3}.dur(1) = [];
stim_table_setup2{97,3}.amp(1) = [];

% Cut 2010 end to ST_nav 5 - 99
stim_table_setup2{99,2}.onset(6) = [];
stim_table_setup2{99,2}.dur(6) = [];
stim_table_setup2{99,2}.amp(6) = [];
stim_table_setup2{99,3}.onset(6) = [];
stim_table_setup2{99,3}.dur(6) = [];
stim_table_setup2{99,3}.amp(6) = [];

% Cut 2020 to ST_nav 3 - 108
stim_table_setup2{108,2}.onset(1) = [];
stim_table_setup2{108,2}.dur(1) = [];
stim_table_setup2{108,2}.amp(1) = [];
stim_table_setup2{108,3}.onset(1) = [];
stim_table_setup2{108,3}.dur(1) = [];
stim_table_setup2{108,3}.amp(1) = [];
stim_table_setup2{108,2}.onset(1) = [];
stim_table_setup2{108,2}.dur(1) = [];
stim_table_setup2{108,2}.amp(1) = [];
stim_table_setup2{108,3}.onset(1) = [];
stim_table_setup2{108,3}.dur(1) = [];
stim_table_setup2{108,3}.amp(1) = [];

% Cut 2030 to ST_nav 2 - 118
stim_table_setup2{118,2}.onset(1) = [];
stim_table_setup2{118,2}.dur(1) = [];
stim_table_setup2{118,2}.amp(1) = [];
stim_table_setup2{118,3}.onset(1) = [];
stim_table_setup2{118,3}.dur(1) = [];
stim_table_setup2{118,3}.amp(1) = [];
stim_table_setup2{118,3}.onset(1) = [];
stim_table_setup2{118,3}.dur(1) = [];
stim_table_setup2{118,3}.amp(1) = [];

job = nirs.modules.ChangeStimulusInfo();
job.ChangeTable = stim_table_setup2;
hb_setup2 = job.run(hb_setup2);

% Trim baseline again
job = nirs.modules.TrimBaseline();
job.preBaseline  = 1;
job.postBaseline = 1;
hb_setup2 = job.run(hb_setup2);

nirs.viz.TimeSeriesViewer(hb_setup2);  

%% Get data for protocol 3

p3_idx = strcmp(demographics.session, "protocol3");
hb_setup3 = hb(p3_idx);
demographics_setup3 = nirs.createDemographicsTable(hb_setup3);
stim_table_setup3 = nirs.createStimulusTable(hb_setup3);

% Cut 1023 to last DT_nav - 21
stim_table_setup3{21,2}.onset(6) = [];
stim_table_setup3{21,2}.dur(6) = [];
stim_table_setup3{21,2}.amp(6) = [];
stim_table_setup3{21,3}.onset(6) = [];
stim_table_setup3{21,3}.dur(6) = [];
stim_table_setup3{21,3}.amp(6) = [];

% Cut 1024 to last DT_nav - 22
stim_table_setup3{22,2}.onset(6) = [];
stim_table_setup3{22,2}.dur(6) = [];
stim_table_setup3{22,2}.amp(6) = [];
stim_table_setup3{22,3}.onset(6) = [];
stim_table_setup3{22,3}.dur(6) = [];
stim_table_setup3{22,3}.amp(6) = [];

job = nirs.modules.ChangeStimulusInfo();
job.ChangeTable = stim_table_setup3;
hb_setup3 = job.run(hb_setup3);

% Trim baseline again
job = nirs.modules.TrimBaseline();
job.preBaseline  = 1;
job.postBaseline = 1;
hb_setup3 = job.run(hb_setup3);

nirs.viz.TimeSeriesViewer(hb_setup3);  

%% GLM protocol 2

% Only CBSI
hb_setup2_cbsi = hb_setup2;
job=nirs.modules.KeepTypes;
job.types={'cbsi'};
hb_setup2_cbsi=job.run(hb_setup2_cbsi);
job = nirs.modules.GLM();
job.type='AR-IRLS';
job.AddShortSepRegressors = true;
SubjStats_cbsi = job.run(hb_setup2_cbsi); 
save('data/mat_files/SubjStats_setup_2_cbsi_tddr_nonpruned.mat','SubjStats_cbsi');
clear SubjStats_cbsi;
clear hb_setup2_cbsi;

%% GLM protocol 3

% Only CBSI
hb_setup3_cbsi = hb_setup3;
job=nirs.modules.KeepTypes;
job.types={'cbsi'};
hb_setup3_cbsi=job.run(hb_setup3_cbsi);
job = nirs.modules.GLM();
job.type='AR-IRLS';
job.AddShortSepRegressors = true;
SubjStats_cbsi = job.run(hb_setup3_cbsi); 
save('data/mat_files/SubjStats_setup_3_cbsi_tddr_nonpruned.mat','SubjStats_cbsi');
clear SubjStats_cbsi;
clear hb_setup1_cbsi;
