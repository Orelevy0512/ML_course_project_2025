% Problematic participants:
% Resting state - 538 (not valid XDF for Before)
% Oddball - 553 (Problem with segmentation)
% Not exist participants - 512, 517, 525

clc;clear all; close all;

%% trial-level analysis
relevant_subjects = [501:511,513:516,518:524,526:537, 539:552,554:563];
% FOOOF PARAMETERS
alpha_down = 6; % Lower bound of alpha band (Hz)
alpha_up = 15; % Upper bound of alpha band (Hz)
alpha_elec = {'P1','P2','P4','PO3','PO4','POz','Pz'}; % To create the data_matrix_64 we changed the electrodes labels to include all 64 channels. 
beta_down = 15; % Lower bound of beta band (Hz)
beta_up = 30; % Upper bound of beta band (Hz)
beta_elec = {'C2', 'Cz', 'F1', 'F2', 'FC1', 'FC2', 'FCz', 'Fz'}; % To create the data_matrix_64 we changed the electrodes labels to include all 64 channels. 
alpha_and_beta_elec = union(alpha_elec,beta_elec);
data_matrix = {};

for part_id = relevant_subjects

    fooof_path = ['C:\data\ZeusData\Orel\Phd\EEG_Hard_times\Results\Resting state Alpha\' num2str(part_id)]
    load([fooof_path, '\data_spectrum_BA_OC_FOOOF_2_40.mat'])

    %% ================ RESTING STATE ===============
    % For RS - we already have trial level analysis, because for each condition
    % (before and after, open and closed) there is only one trial. 

    open_paradigms_fooof = {'data_spectrum_Before_Open', 'data_spectrum_After_Open'};
    paradigm_id = 1; % RS_OPEN

    for trial_id = 1:length(open_paradigms_fooof)
        feature_matrix_alpha = []; feature_matrix_beta = [];
        fooof_data = eval(open_paradigms_fooof{trial_id});
        % FoooF Vector Extraction 
        [~, electrode_indices] = ismember(alpha_and_beta_elec, fooof_data.label);
        FOOOF_Vector = mean(fooof_data.powspctrm(electrode_indices, :));
        FOOOF_Vector = FOOOF_Vector(~isnan(FOOOF_Vector));

        % Alpha 
        [~, electrode_indices] = ismember(alpha_elec, fooof_data.label); % Find the indices of the specified electrodes in the data
        relevant_freq_indices_fooof = find(fooof_data.freq > alpha_down & fooof_data.freq < alpha_up);
        power_values_in_alpha_fooof = fooof_data.powspctrm(electrode_indices, relevant_freq_indices_fooof);
        mean_power_alpha_fooof = mean(power_values_in_alpha_fooof,2);
        [max_power_fooof, max_index_fooof] = max(power_values_in_alpha_fooof,[],2);
        frequency_at_max_power_fooof = fooof_data.freq(relevant_freq_indices_fooof(max_index_fooof));

        feature_matrix_alpha(:,1) = frequency_at_max_power_fooof'; % Alpha Freq for each electrode
        feature_matrix_alpha(:,2) = max_power_fooof; % Alpha max for each electrode
        feature_matrix_alpha(:,3) = mean_power_alpha_fooof; % Alpha mean for each electrode

        % Beta 
        [~, electrode_indices] = ismember(beta_elec, fooof_data.label); % Find the indices of the specified electrodes in the data
        relevant_freq_indices_fooof = find(fooof_data.freq > beta_down & fooof_data.freq < beta_up);
        power_values_in_beta_fooof = fooof_data.powspctrm(electrode_indices, relevant_freq_indices_fooof);
        mean_power_beta_fooof = mean(power_values_in_beta_fooof,2);
        [max_power_fooof, max_index_fooof] = max(power_values_in_beta_fooof,[],2);
        frequency_at_max_power_fooof = fooof_data.freq(relevant_freq_indices_fooof(max_index_fooof));

        feature_matrix_beta(:,1) = frequency_at_max_power_fooof'; % Alpha Freq for each electrode
        feature_matrix_beta(:,2) = max_power_fooof; % Alpha max for each electrode
        feature_matrix_beta(:,3) = mean_power_beta_fooof; % Alpha mean for each electrode  

        data_matrix(end+1, :) = {
            part_id, paradigm_id, trial_id, feature_matrix_alpha, feature_matrix_beta, FOOOF_Vector
        };
    end


    closed_paradigms_fooof = {'data_spectrum_Before_Closed', 'data_spectrum_After_Closed'};
    paradigm_id = 2; % RS_Closed

    for trial_id = 1:length(closed_paradigms_fooof)
        feature_matrix_alpha = []; feature_matrix_beta = [];
        fooof_data = eval(closed_paradigms_fooof{trial_id});
        % FoooF Vector Extraction 
        [~, electrode_indices] = ismember(alpha_and_beta_elec, fooof_data.label);
        FOOOF_Vector = mean(fooof_data.powspctrm(electrode_indices, :));
        FOOOF_Vector = FOOOF_Vector(~isnan(FOOOF_Vector));

        % Alpha 
        [~, electrode_indices] = ismember(alpha_elec, fooof_data.label); % Find the indices of the specified electrodes in the data
        relevant_freq_indices_fooof = find(fooof_data.freq > alpha_down & fooof_data.freq < alpha_up);
        power_values_in_alpha_fooof = fooof_data.powspctrm(electrode_indices, relevant_freq_indices_fooof);
        mean_power_alpha_fooof = mean(power_values_in_alpha_fooof,2);
        [max_power_fooof, max_index_fooof] = max(power_values_in_alpha_fooof,[],2);
        frequency_at_max_power_fooof = fooof_data.freq(relevant_freq_indices_fooof(max_index_fooof));

        feature_matrix_alpha(:,1) = frequency_at_max_power_fooof'; % Alpha Freq for each electrode
        feature_matrix_alpha(:,2) = max_power_fooof; % Alpha max for each electrode
        feature_matrix_alpha(:,3) = mean_power_alpha_fooof; % Alpha mean for each electrode

        % Beta 
        [~, electrode_indices] = ismember(beta_elec, fooof_data.label); % Find the indices of the specified electrodes in the data
        relevant_freq_indices_fooof = find(fooof_data.freq > beta_down & fooof_data.freq < beta_up);
        power_values_in_beta_fooof = fooof_data.powspctrm(electrode_indices, relevant_freq_indices_fooof);
        mean_power_beta_fooof = mean(power_values_in_beta_fooof,2);
        [max_power_fooof, max_index_fooof] = max(power_values_in_beta_fooof,[],2);
        frequency_at_max_power_fooof = fooof_data.freq(relevant_freq_indices_fooof(max_index_fooof));

        feature_matrix_beta(:,1) = frequency_at_max_power_fooof'; % Alpha Freq for each electrode
        feature_matrix_beta(:,2) = max_power_fooof; % Alpha max for each electrode
        feature_matrix_beta(:,3) = mean_power_beta_fooof; % Alpha mean for each electrode  

        data_matrix(end+1, :) = {
            part_id, paradigm_id, trial_id, feature_matrix_alpha, feature_matrix_beta, FOOOF_Vector
        };
    end

    %% ================ Oddball ===============

    variables_folder = (['\\AFRODITA\LabBackup\Rotem\Results\Oddball\', num2str(part_id)]); 
    %main_output_folder = 'C:\data\ZeusData\Orel\Phd\EEG_Hard_times\Results\Oddball Alpha'
    %cd(main_output_folder)
    %mkdir(num2str(part_id))
    %mainData_folder = ([main_output_folder,'\', num2str(part_id)]); 

    load([variables_folder, '\data_afterICA.mat']) 
    load([variables_folder, '\trl_new.mat']) 
    number_of_trials = size(trl_new,1)/80;

    new_matrix=[];
    for trial = 1:number_of_trials
        start_idx = (trial - 1) * 80 + 1; % Index for the first event of the trial
        end_idx = trial * 80; % Index for the last event of the trial

        % Get the start time of the first event in the trial
        trial_start = trl_new(start_idx, 1);

        % Get the end time of the last event in the trial
        trial_end = trl_new(end_idx, 2);

        % Store the start and end times in the new matrix
        new_matrix(trial, 1) = trial_start;
        new_matrix(trial, 2) = trial_end;
        new_matrix(trial, 3) = -102.4;
        new_matrix(trial, 4) = 0;
    end
    %% segment to trials
    cfg = [];
    cfg.trl = new_matrix;
    dataAfterICA_epochs = ft_redefinetrial(cfg, data_afterICA);

    dataaftercleaning=dataAfterICA_epochs;
    %% Prepare data to analysis - align all trials to be in the same length 
    smallest_trial_in_samples = min(dataaftercleaning.sampleinfo(:,2)-dataaftercleaning.sampleinfo(:,1));
    newTRL = dataaftercleaning.sampleinfo;
    newTRL(:,3) = zeros(length(newTRL),1);
    newTRL(:,4) = dataaftercleaning.trialinfo

    for i=1:length(dataaftercleaning.sampleinfo)
        diff = dataaftercleaning.sampleinfo(i,2)-dataaftercleaning.sampleinfo(i,1)
        if diff-smallest_trial_in_samples>0
            newTRL(i,2) = newTRL(i,2) - (diff-smallest_trial_in_samples)
        end
    end
    cfg = [];
    cfg.trl = newTRL;
    data_afterICA_fixed = ft_redefinetrial(cfg,dataaftercleaning);
    paradigm_id = 3 % Oddball
    %% FOOOF 2-40
    for trial_id=1:length(data_afterICA_fixed.trial)
        fooof_data = {}; feature_matrix_alpha = []; feature_matrix_beta = [];
        cfg= [];
        cfg.method = 'mtmfft';
        cfg.output     = 'fooof_peaks';
        cfg.foilim = [2 40];
        cfg.taper = 'hanning';
        cfg.trials = trial_id
        fooof_data = ft_freqanalysis(cfg,data_afterICA_fixed);

        % FoooF Vector Extraction 
        [~, electrode_indices] = ismember(alpha_and_beta_elec, fooof_data.label);
        FOOOF_Vector = mean(fooof_data.powspctrm(electrode_indices, :));
        FOOOF_Vector = FOOOF_Vector(~isnan(FOOOF_Vector));

        % Alpha 
        [~, electrode_indices] = ismember(alpha_elec, fooof_data.label); % Find the indices of the specified electrodes in the data
        relevant_freq_indices_fooof = find(fooof_data.freq > alpha_down & fooof_data.freq < alpha_up);
        power_values_in_alpha_fooof = fooof_data.powspctrm(electrode_indices, relevant_freq_indices_fooof);
        mean_power_alpha_fooof = mean(power_values_in_alpha_fooof,2);
        [max_power_fooof, max_index_fooof] = max(power_values_in_alpha_fooof,[],2);
        frequency_at_max_power_fooof = fooof_data.freq(relevant_freq_indices_fooof(max_index_fooof));

        feature_matrix_alpha(:,1) = frequency_at_max_power_fooof'; % Alpha Freq for each electrode
        feature_matrix_alpha(:,2) = max_power_fooof; % Alpha max for each electrode
        feature_matrix_alpha(:,3) = mean_power_alpha_fooof; % Alpha mean for each electrode

        % Beta 
        [~, electrode_indices] = ismember(beta_elec, fooof_data.label); % Find the indices of the specified electrodes in the data
        relevant_freq_indices_fooof = find(fooof_data.freq > beta_down & fooof_data.freq < beta_up);
        power_values_in_beta_fooof = fooof_data.powspctrm(electrode_indices, relevant_freq_indices_fooof);
        mean_power_beta_fooof = mean(power_values_in_beta_fooof,2);
        [max_power_fooof, max_index_fooof] = max(power_values_in_beta_fooof,[],2);
        frequency_at_max_power_fooof = fooof_data.freq(relevant_freq_indices_fooof(max_index_fooof));

        feature_matrix_beta(:,1) = frequency_at_max_power_fooof'; % Alpha Freq for each electrode
        feature_matrix_beta(:,2) = max_power_fooof; % Alpha max for each electrode
        feature_matrix_beta(:,3) = mean_power_beta_fooof; % Alpha mean for each electrode  

        data_matrix(end+1, :) = {
            part_id, paradigm_id, trial_id, feature_matrix_alpha, feature_matrix_beta, FOOOF_Vector
        };
    end

    %% ================ Lecture ===================

    load(['C:\data\ZeusData\Orel\Phd\EEG_Hard_times\Results\Lecture\',num2str(part_id),'\full_speech_trf_Noclean\newData']);
    dataaftercleaning=newData;
    %% Prepare data to analysis - align all trials to be in the same length 
    smallest_trial_in_samples = min(dataaftercleaning.sampleinfo(:,2)-dataaftercleaning.sampleinfo(:,1));
    newTRL = dataaftercleaning.sampleinfo;
    newTRL(:,3) = zeros(length(newTRL),1);
    newTRL(:,4) = dataaftercleaning.trialinfo

    for i=1:length(dataaftercleaning.sampleinfo)
        diff = dataaftercleaning.sampleinfo(i,2)-dataaftercleaning.sampleinfo(i,1)
        if diff-smallest_trial_in_samples>0
            newTRL(i,2) = newTRL(i,2) - (diff-smallest_trial_in_samples)
        end
    end
    cfg = [];
    cfg.trl = newTRL;
    data_afterICA_fixed = ft_redefinetrial(cfg,dataaftercleaning);
    paradigm_id = 4 % Lecture
    %% FOOOF 2-40
    for trial_id=1:length(data_afterICA_fixed.trial)
        fooof_data = {}; feature_matrix_alpha = []; feature_matrix_beta = [];
        cfg= [];
        cfg.method = 'mtmfft';
        cfg.output     = 'fooof_peaks';
        cfg.foilim = [2 40];
        cfg.taper = 'hanning';
        cfg.trials = trial_id
        fooof_data = ft_freqanalysis(cfg,data_afterICA_fixed);

        % FoooF Vector Extraction 
        [~, electrode_indices] = ismember(alpha_and_beta_elec, fooof_data.label);
        FOOOF_Vector = mean(fooof_data.powspctrm(electrode_indices, :));
        FOOOF_Vector = FOOOF_Vector(~isnan(FOOOF_Vector));

        % Alpha 
        [~, electrode_indices] = ismember(alpha_elec, fooof_data.label); % Find the indices of the specified electrodes in the data
        relevant_freq_indices_fooof = find(fooof_data.freq > alpha_down & fooof_data.freq < alpha_up);
        power_values_in_alpha_fooof = fooof_data.powspctrm(electrode_indices, relevant_freq_indices_fooof);
        mean_power_alpha_fooof = mean(power_values_in_alpha_fooof,2);
        [max_power_fooof, max_index_fooof] = max(power_values_in_alpha_fooof,[],2);
        frequency_at_max_power_fooof = fooof_data.freq(relevant_freq_indices_fooof(max_index_fooof));

        feature_matrix_alpha(:,1) = frequency_at_max_power_fooof'; % Alpha Freq for each electrode
        feature_matrix_alpha(:,2) = max_power_fooof; % Alpha max for each electrode
        feature_matrix_alpha(:,3) = mean_power_alpha_fooof; % Alpha mean for each electrode

        % Beta 
        [~, electrode_indices] = ismember(beta_elec, fooof_data.label); % Find the indices of the specified electrodes in the data
        relevant_freq_indices_fooof = find(fooof_data.freq > beta_down & fooof_data.freq < beta_up);
        power_values_in_beta_fooof = fooof_data.powspctrm(electrode_indices, relevant_freq_indices_fooof);
        mean_power_beta_fooof = mean(power_values_in_beta_fooof,2);
        [max_power_fooof, max_index_fooof] = max(power_values_in_beta_fooof,[],2);
        frequency_at_max_power_fooof = fooof_data.freq(relevant_freq_indices_fooof(max_index_fooof));

        feature_matrix_beta(:,1) = frequency_at_max_power_fooof'; % Alpha Freq for each electrode
        feature_matrix_beta(:,2) = max_power_fooof; % Alpha max for each electrode
        feature_matrix_beta(:,3) = mean_power_beta_fooof; % Alpha mean for each electrode  

        data_matrix(end+1, :) = {
            part_id, paradigm_id, trial_id, feature_matrix_alpha, feature_matrix_beta, FOOOF_Vector
        };
    end
end