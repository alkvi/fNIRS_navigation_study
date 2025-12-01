import os
import mne
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_theme()
import pandas as pd
import mne_nirs
from mne_bids import BIDSPath, read_raw_bids

mne.viz.set_3d_backend('pyvistaqt')

def plot_bar(raw_intensity, quality_data, xlim, vline, label, fig_name):
    channel_names = raw_intensity.ch_names
    wl_channel_names = [name for name in channel_names if "760" in name]
    wl_idx = mne.pick_channels(channel_names, include=wl_channel_names)

    # Just take the first wavelength (both are the same)
    channel_names = [channel_names[i] for i in wl_idx]
    quality_data = [quality_data[i] for i in wl_idx]

    sns.set(rc = {'figure.figsize':(15,10)})
    sns.set_color_codes("pastel")
    ax = sns.barplot(x=quality_data, y=channel_names, data=None, color="b")
    ax.axvline(x=vline, color='r', linestyle='dashed')
    ax.set(xlim=(0, xlim), ylabel="",
        xlabel=label)
    sns.despine(left=True, bottom=True)
    plt.savefig(fig_name)
    plt.close()

if __name__ == "__main__":

    # Specify BIDS root folder
    bids_root = "../../Park-MOVE_fnirs_dataset_v2/fNIRS_data/bids_dataset_snirf"
    plot_subject = False
    sci_threshold = 0.7
    psp_threshold = 0.1

    # Calculated quality data
    sci_psp_df = pd.read_csv("../data/sci_psp_nirs_toolbox.csv")

    # We have 4 file types: events, channels, optodes, and nirs.
    datatype = 'nirs'
    bids_path = BIDSPath(root=bids_root, datatype=datatype)
    nirs_files = bids_path.match()

    # For visualization
    subjects_dir = str(mne.datasets.sample.data_path()) + '/subjects'
    mne.datasets.fetch_fsaverage(subjects_dir=subjects_dir)
    mne.set_config('MNE_NIRS_FOLD_PATH', '/Users/alexander.kvist/OneDrive\ -\ Karolinska\ Institutet/Dokument/Software/fOLD-public/Supplementary')

    # Prepare params
    task = 'complexwalk'
    suffix = 'nirs'
    session = 'protocol3'

    # Only use data in clusters
    study_subjects = pd.read_csv("../data/clusters_oa.csv")['subject'].to_list()
    study_subjects.extend(pd.read_csv("../data/clusters_pd.csv")['subject'].to_list())
    print("Subjects")
    print(study_subjects)

    # Load montage info from one subject
    subject = study_subjects[0]
    bids_path = BIDSPath(subject=subject, task=task, session=session,
                        suffix=suffix, datatype=datatype, root=bids_root)
    raw_intensity = read_raw_bids(bids_path=bids_path, verbose=False)

    # Add landmark positions from digpts.txt
    montage = raw_intensity.get_montage()
    new_montage = mne.channels.make_dig_montage(montage.get_positions()['ch_pos'], 
                                                nasion=np.array([0.400, 85.900, -47.600])/1000,
                                                lpa=np.array([-83.800, -18.600, -57.200])/1000, 
                                                rpa=np.array([83.900, -16.600, -56.700])/1000)

    raw_intensity.set_montage(new_montage)
    coreg = mne.coreg.Coregistration(
        raw_intensity.info, "fsaverage", subjects_dir, fiducials="estimated"
    )
    coreg.fit_fiducials(lpa_weight=1.0, nasion_weight=1.0, rpa_weight=1.0)

    # Which channels are short channels? Have 16 (8 for 850, 8 for 760)
    # Some short-detectors (virtual detector 8-15) make up links 
    # with not just the source they're attached to, but ones further away.
    # these are still short-separation (so e.g. 1-14 and 1-15 detect the same
    # signal). Mark these also. 
    long_ch_names = np.array(raw_intensity.ch_names)
    detector_long = ["D8", "D9", "D10", "D11", "D12", "D13", "D14", "D15"]
    long_ch_names = [str(name) for name in long_ch_names if not any(code in name for code in detector_long)]
    long_ch_idx = mne.pick_channels(raw_intensity.ch_names, include=long_ch_names)

    # Which pairs should we get?
    long_ch_names_unique = np.unique([channel.replace(" 760", "").replace(" 850", "") for channel in long_ch_names]).tolist()
    pairs = [(int(label.split('_')[0][1:]), int(label.split('_')[1][1:])) for label in long_ch_names_unique]

    # Pick only the long channels
    raw_intensity = raw_intensity.copy().pick(picks=long_ch_idx)

    # Save data
    sci_array = []
    bad_ch_array_sci = []
    power_array = []
    bad_ch_array_psp = []
    bad_ch_array_combined = []

    # Go through data for each subject
    for subject in study_subjects:
        
        # Get corresponding SCI and PSP from table
        print(f"Handling subject {subject} session {session}")
        filtered_df = sci_psp_df[(sci_psp_df['subject'] == subject) & (sci_psp_df['session'] == session)]
        filtered_df = filtered_df[filtered_df.apply(lambda row: (row['source'], row['detector']) in pairs, axis=1)]
        if filtered_df.empty:
            print(f"Subject {subject} session {session} has no quality data, skipping")
            continue

        # These are in the same order as MNE channels. 
        overall_sci = filtered_df['sci'].to_numpy()
        peak_power = filtered_df['power'].to_numpy()

        # Repeat each value to get one value per wavelength, as MNE expects
        overall_sci = np.array(np.tile(overall_sci, 2).tolist())
        peak_power = np.array(np.tile(peak_power, 2).tolist())

        # Append to lists
        sci_array.append(overall_sci)
        power_array.append(peak_power)
        
        # Mark channels for subject/session
        # Combine both criteria
        combined_mask = (overall_sci < sci_threshold) | (peak_power < psp_threshold)

        # Store logical list of bad channel idx in array
        bad_ch_array_sci.append(overall_sci < sci_threshold)
        bad_ch_array_psp.append(peak_power < psp_threshold)
        bad_ch_array_combined.append(combined_mask)

        if plot_subject:
            # Plot overall SCI
            fig_name = f"figures_quality/SCI_overall_{bids_path.basename}_{session}.png"
            plot_bar(raw_intensity, overall_sci, 1.0, sci_threshold, "SCI", fig_name)

            # Plot overall PSP
            fig_name = f"figures_quality/PSP_overall_{bids_path.basename}_{session}.png"
            plot_bar(raw_intensity, peak_power, 0.3, psp_threshold, "PSP", fig_name)


# How many subjects did we calculate for?
num_subjects = len(sci_array)
print(f"Calculated for {num_subjects}")

# SCI averages
sci_array = np.array(sci_array)
sci_array_avg = np.mean(sci_array, axis=0)
bad_ch_array_sci = np.array(bad_ch_array_sci)
percent_bad_ch_sci = np.sum(bad_ch_array_sci, axis=0) / num_subjects * 100

# PSP averages
power_array = np.array(power_array)
power_array_avg = np.mean(power_array, axis=0)
bad_ch_array_psp = np.array(bad_ch_array_psp)
percent_bad_ch_psp = np.sum(bad_ch_array_psp, axis=0) / num_subjects * 100

# Combined averages
bad_ch_array_combined = np.array(bad_ch_array_combined)
percent_bad_ch_combined = np.sum(bad_ch_array_combined, axis=0) / num_subjects * 100

# SCI, average
fig = mne_nirs.visualisation.plot_nirs_source_detector(
    sci_array_avg,
    raw_intensity.info,
    surfaces="brain",
    subject="fsaverage",
    subjects_dir=subjects_dir,
    trans="fsaverage",
)
# Bad channels SCI, number of
fig = mne_nirs.visualisation.plot_nirs_source_detector(
    percent_bad_ch_sci,
    raw_intensity.info,
    surfaces="brain",
    subject="fsaverage",
    subjects_dir=subjects_dir,
    trans="fsaverage",
)
# Peak power, average
fig = mne_nirs.visualisation.plot_nirs_source_detector(
    power_array_avg,
    raw_intensity.info,
    surfaces="brain",
    subject="fsaverage",
    subjects_dir=subjects_dir,
    trans="fsaverage",
)
# Bad channels PSP, number of
fig = mne_nirs.visualisation.plot_nirs_source_detector(
    percent_bad_ch_psp,
    raw_intensity.info,
    surfaces="brain",
    subject="fsaverage",
    subjects_dir=subjects_dir,
    trans="fsaverage",
)
# Bad channels combined, number of
fig = mne_nirs.visualisation.plot_nirs_source_detector(
    percent_bad_ch_combined,
    raw_intensity.info,
    surfaces="brain",
    subject="fsaverage",
    subjects_dir=subjects_dir,
    trans="fsaverage",
)
plt.show()
input("ENTER WHEN DONE")