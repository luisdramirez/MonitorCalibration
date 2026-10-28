%%% MonitorCalibration
% ! Make sure there is no ambient light when using analog photometer!

% Measure luminance with photometer and CalibrateMonitorPhotometer().
% The output of CalibrateMonitorPhotometer() will be saved as a .mat file ___.

% To use, load the .mat file ___ and call
% Screen('LoadNormalizedGammaTable', window, correctedGammaTable);

clear all; close all; clc
KbName('UnifyKeyNames');

%% Define display/setup

% which_setup = input('Enter a setup/display name: ');
which_setup = '3329C_ASUS';

%% Measure luminance with photometer

num_lum_measures = 9;

[gamma_table1, ~, display_baseline, display_range, display_gamma, max_level] = CalibrateMonitorPhotometer(num_lum_measures);

%% Verify corrected gamma table from the display's gamma

num_levels = 256;

corrected_gamma_table = linspace(0,1,num_levels)' .^ (1 / display_gamma);
corrected_gamma_table = repmat(corrected_gamma_table, 1, 3);

%% Visualize gamma table

figure('Color',[1 1 1])

plot(gamma_table1)
hold on
plot(corrected_gamma_table(:,1))
box off;
set(gca,'TickDir','out')

legend({'orig gamma', 'new gamma'});

%% Store corrected gamma

if sum(corrected_gamma_table(:,1) == gamma_table1) == length(gamma_table1)
    disp('Gamma corrected.')
    corrected_gamma.table = corrected_gamma_table;
else
    error('Redo Calibration.')
end

corrected_gamma.display_baseline = display_baseline;
corrected_gamma.display_range = display_range;
corrected_gamma.display_gamma = display_gamma;
corrected_gamma.max_level = max_level;

%% Save gamma table

script_dir = pwd;

% Define save dir 
save_dir = '/home/serenceslabexp/Desktop/MonitorCalibration/GammaTables';

% Check if save dir exists
if ~exist(save_dir,'dir')
    error(['Create ' save_dir]);
else
    disp('/GammaTable exists.');
end

% Define filename
corrected_gamma_filename = ['corrected_gamma_table_' which_setup '.mat'];

% Save corrected gamma table
cd(save_dir)
save(corrected_gamma_filename,'corrected_gamma')
cd(script_dir)