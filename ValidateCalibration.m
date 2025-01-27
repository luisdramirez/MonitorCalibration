%%% ConfirmCalibration

clear all; close all; clc;
KbName('UnifyKeyNames');

%% Define display/setup 

which_setup = '3329C_ASUS';

switch which_setup
    case '3329C_ASUS'
        
        w.refresh_rate = 59.95; % display refresh rate in Hertz, Hz
        w.view_distance = 40;  % cm
        w.screen_width = 58.5; %  cm
        w.screen_width_px = 2560; % px
        w.screen_height_px = 1440; % px

    case '3329C_HP'

        w.refresh_rate = 85; % display refresh rate in Hertz, Hz
        w.view_distance = 40;  % cm
        w.screen_width = 40.7; %  cm
        w.screen_width_px = 1600; % px
        w.screen_height_px = 1200; % px

end

%% Calculate visual angle of the display, pixels per degree of visual angle, size of a pixel, and Nyquist frequency

screen_length = w.screen_width;
screen_length_px = w.screen_width_px;

w.visual_angle = 2 * atan2d(screen_length/2,  w.view_distance); % Visual angle of the whole screen in degrees
w.ppd = floor(screen_length_px/w.visual_angle); % Pixels per degree of visual angle
w.px_size = screen_length/screen_length_px; % size of pixel in cm
w.f_Nyquist = 1/(2*w.px_size);

%% Open window

[window, ~] = Screen('OpenWindow', 0, [0, 0, 0]); % Black background

%% Load gamma corrected table

corrected_gamma_dir = '/home/serenceslabexp/Desktop/MonitorCalibration/GammaTables';
corrected_gamma_filename = ['corrected_gamma_table_' which_setup '.mat'];

cd(corrected_gamma_dir)
load(corrected_gamma_filename)
cd(script_dir)

Screen('LoadNormalizedGammaTable', window, corrected_gamma.table)

%% Display a gradient pattern

% Create a horizontal gradient (grayscale values from 0 to 255)
gradient = repmat(linspace(0, 1, 256), 256, 1); % 256 x 256 gradient
gradient_texture = Screen('MakeTexture', window, gradient * 255); % Scale to 0–255 (grayscale)

% Display the gradient
Screen('DrawTexture', window, gradient_texture);
Screen('Flip', window);
