import os 
import pandas as pd 
import numpy as np 

if __name__ == "__main__":

    pd.set_option('display.max_rows', None)

    # Load data
    ml_data_file = "../data/exported_mobility_lab_data/trial_summaries/Walk_trials.csv"
    ml_data = pd.read_csv(ml_data_file, encoding="utf-8", header=0)  
    regex_pattern = "(\d+)(?=\s*valid gait cycles)"
    ml_data['gait_cycles'] = ml_data['Analysis Log'].str.extract(regex_pattern, expand=False)

    ml_data = ml_data.rename(columns={'Subject Public ID': 'subject', 
                                      'Turns - Turn Velocity (degrees/s) [mean]': 'turns_velocity',
                                      'Turns - Turn Velocity (degrees/s) [std]': 'turns_velocity_std',
                                      'Turns - Steps in Turn (#) [mean]': 'turns_steps',
                                      'Turns - Steps in Turn (#) [std]': 'turns_steps_std',
                                      'Turns - Duration (s) [mean]': 'turns_duration',
                                      'Turns - Duration (s) [std]': 'turns_duration_std',
                                      'Turns - Angle (degrees) [mean]': "turns_angle",
                                      'Turns - Angle (degrees) [std]': "turns_angle_std"})
    ml_data = ml_data[['subject', 'turns_velocity', 'turns_velocity_std', 'turns_steps', 'turns_steps_std', 
                       'turns_duration', 'turns_duration_std', 'turns_angle', 'turns_angle_std', 'gait_cycles']]
    ml_data = ml_data.apply(lambda x: x.str.replace(',', '.'))
    ml_data['subject'] = ml_data['subject'].str.replace('FNP2040NEW', 'FNP2040')
    
    # We know the file order for some subjects is different
    missing_calibration = ["FNP1001", "FNP1005", "FNP1006", "FNP1007"]
    calibration_last = ["FNP2023", "FNP2025", "FNP2037"]
    missing_p1 = ["FNP1054", "FNP1056"]
    missing_p3 = ["FNP1044", "FNP1072", "FNP2036", "FNP2040"]

    # Set up a mapping from file to protocol
    file_id_mapping_normal = {1: "calibration_1", 2: "calibration_2", 3: "protocol_1", 4: "protocol_2", 5: "protocol_3"}
    file_id_mapping_missing_calibration = {1: "protocol_1", 2: "protocol_2", 3: "protocol_3"}
    file_id_mapping_last = {1: "calibration_2", 2: "protocol_1", 3: "protocol_2", 4: "protocol_3", 5: "calibration_1"}
    file_id_mapping_missing_p1 = {1: "calibration_1", 2: "calibration_2", 3: "protocol_2", 4: "protocol_3"}
    file_id_mapping_missing_p3 = {1: "calibration_1", 2: "calibration_2", 3: "protocol_1", 4: "protocol_2"}
    file_id_mapping_1033 = {1: "calibration_1", 2: "calibration_2", 3: "protocol_2", 4: "protocol_1", 5: "protocol_3"}
    file_id_mapping_2040 = {1: "calibration_1", 2: "calibration_2", 3: "protocol_1", 4: "calibration_1_2", 5: "calibration_2_2", 6: "protocol_2"}
    file_id_mapping_normal = {v: k for k, v in file_id_mapping_normal.items()}
    file_id_mapping_missing_calibration = {v: k for k, v in file_id_mapping_missing_calibration.items()}
    file_id_mapping_last = {v: k for k, v in file_id_mapping_last.items()}
    file_id_mapping_missing_p1 = {v: k for k, v in file_id_mapping_missing_p1.items()}
    file_id_mapping_missing_p3 = {v: k for k, v in file_id_mapping_missing_p3.items()}
    file_id_mapping_1033 = {v: k for k, v in file_id_mapping_1033.items()}
    file_id_mapping_2040 = {v: k for k, v in file_id_mapping_2040.items()}

    # Go through each subject and get metrics for protocol 2 and 3
    subjects = np.unique(ml_data['subject'])
    all_subj_data = []
    for subject_id in subjects:
        subj_data = ml_data[ml_data['subject'] == subject_id]

        for protocol in ['protocol_2', 'protocol_3']:

            if protocol == 'protocol_3' and subject_id in missing_p3:
                continue

            # Which protocol are we processing?
            if subject_id == "FNP1033":
                file_protocol = file_id_mapping_1033[protocol]
            elif subject_id == "FNP2040":
                file_protocol = file_id_mapping_2040[protocol]
            elif subject_id in missing_p1:
                file_protocol = file_id_mapping_missing_p1[protocol]
            elif subject_id in missing_p3:
                file_protocol = file_id_mapping_missing_p3[protocol]
            elif subject_id in missing_calibration:
                file_protocol = file_id_mapping_missing_calibration[protocol]
            elif subject_id in calibration_last:
                file_protocol = file_id_mapping_last[protocol]
            else:
                file_protocol = file_id_mapping_normal[protocol]

            subj_df = subj_data.iloc[[file_protocol-1]].copy()
            subj_df['protocol'] = protocol
            all_subj_data.append(subj_df)

    p2_p3_data = pd.concat(all_subj_data)

    # Go through each subject
    ya_ids  = pd.read_csv("../../Park-MOVE_fnirs_dataset_v1/identifiers_YA.csv", encoding="utf-8", header=0)['id_nummer'].to_list()
    oa_ids  = pd.read_csv("../../Park-MOVE_fnirs_dataset_v1/identifiers_OA.csv", encoding="utf-8", header=0)['id_nummer'].to_list()
    pd_ids  = pd.read_csv("../../Park-MOVE_fnirs_dataset_v1/identifiers_PD.csv", encoding="utf-8", header=0)['id_nummer'].to_list()

    # Assign groups
    for id in ya_ids:
        p2_p3_data.loc[p2_p3_data["subject"] == id, 'group'] = "YA"
    for id in oa_ids:
        p2_p3_data.loc[p2_p3_data["subject"] == id, 'group'] = "OA"
    for id in pd_ids:
        p2_p3_data.loc[p2_p3_data["subject"] == id, 'group'] = "PD"

    print(p2_p3_data)

    # Save DF 
    p2_p3_data.to_csv("../data/mobility_lab_turn_parameters.csv", index=False)