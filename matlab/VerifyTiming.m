%%% VerifyTiming

clear all; close all; clc;
KbName('UnifyKeyNames');

%%

script_dir = pwd;

%% Define display/setup 

which_setup = '3329C_ASUS';

half_screen = 0;

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

if half_screen
    w.screen_width_px = w.screen_width_px / 2; % px
    w.screen_height_px = w.screen_height_px / 2; % px
end

black = [0 0 0];
white = [255 255 255];
gray = floor(white/2);

%% Calculate visual angle of the display, pixels per degree of visual angle, size of a pixel, and Nyquist frequency

screen_length = w.screen_width;
screen_length_px = w.screen_width_px;

w.visual_angle = 2 * atan2d(screen_length/2,  w.view_distance); % Visual angle of the whole screen in degrees
w.ppd = floor(screen_length_px/w.visual_angle); % Pixels per degree of visual angle
w.px_size = screen_length/screen_length_px; % size of pixel in cm
w.f_Nyquist = 1/(2*w.px_size);

%% Open window

[window, screen_size_px] = Screen('OpenWindow', 0, 127*[1 1 1], [0 0 w.screen_width_px w.screen_height_px]);
default_gamma = Screen('LoadNormalizedGammaTable',window);
HideCursor;
commandwindow;

% Turn on alpha blending
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

% Get center coords
centerX = screen_size_px(3) / 2;
centerY = screen_size_px(4) / 2;

%% Load CLUT

corrected_gamma_dir = '/home/serenceslabexp/Desktop/MonitorCalibration/GammaTables';
corrected_gamma_filename = ['corrected_gamma_table_' which_setup '.mat'];

cd(corrected_gamma_dir);
load(corrected_gamma_filename);
cd(script_dir);

Screen('LoadNormalizedGammaTable', window, corrected_gamma.table);

%% Create Fixation
% Draw fixation for t seconds, frame by frame.
% The fixation dot changes color every N frames

% Define fixation
fixation_dot_deg = 0.15;
fixation_dot_px = round(w.ppd * fixation_dot_deg);
if ~mod(fixation_dot_px,2), fixation_dot_px = fixation_dot_px + 1; end
fixation_dot_color = [black; white];
fixation_dot_patch = CenterRectOnPoint([0 0 fixation_dot_px fixation_dot_px], centerX, centerY);

% Set timing
frame_rate = 20; % Hz
frame_dur = 1/frame_rate; % s
fixation_dur = 5; % s
fixation_color_dur = fixation_dur / 10;

% Define frame onsets
frame_onsets_baseline = 0:frame_dur:fixation_dur-frame_dur;
num_frames = length(frame_onsets_baseline);

% Define fixation color per frame
num_frames_per_color = round(fixation_color_dur/frame_dur);

fixation_color_update = zeros(1,num_frames);
fixation_color_update(num_frames_per_color+1:num_frames_per_color:num_frames) = 1;

%% Draw fixation

% Pre-allocate vectors for timing measures
VBL_Timestamp = nan(1,num_frames);
StimulusOnsetTime = nan(1,num_frames);
FlipTimeStamp = nan(1,num_frames);
Missed = nan(1,num_frames);

% Draw fixation
curr_color = 1;
for n_frame = 1:num_frames

    if n_frame == 1
        frame_onsets = frame_onsets_baseline + GetSecs;
    end
    
    % Change fixation color when it's time
    if fixation_color_update(n_frame)
        if curr_color == 1
            curr_color = 2;
        elseif curr_color == 2
            curr_color = 1;
        end
    end

    % Fixation dot
    Screen('FillOval', window, fixation_dot_color(curr_color,:), fixation_dot_patch);

    % Flip
    [VBL_Timestamp(n_frame), StimulusOnsetTime(n_frame), FlipTimeStamp(n_frame), Missed(n_frame)] = Screen('Flip', window, frame_onsets(n_frame));

end

%% Evaluate fixation timing
% VBLTimestamp': a high-precision estimate of the system time (in seconds) when the actual flip has happened
% 'StimulusOnsetTime': An estimate of Stimulus-onset time
% 'FlipTimestamp' is a timestamp taken at the end of Flip's execution. 
% Use the difference between FlipTimestamp and VBLTimestamp to get an estimate of how long Flips execution takes.
 
% 'Missed' indicates if the requested presentation deadline for your stimulus has been missed. 
% A negative value means that deadlines have been satisfied. 
% Positive values indicate a deadline-miss.

% Set figure path
figure_path = 'TimingLogs';
if ~exist(figure_path,'dir')
    sca;
    error([figure_path ' does not exist!']);
end

% Summarize frames missed
frames_missed = sum(Missed > 0);

% Plot frames missed
figure_name = 'Fixation missed frames';
fg = figure('Color', [1 1 1], 'Name', figure_name);
plot(Missed * 10^3);
hold on
line([0 max(xlim)], [0 0], 'LineStyle','--', 'Color', black)

% Format figure
title(['Frames missed: ' num2str(frames_missed) ' (' num2str(100*frames_missed/num_frames) '%)']);
xlabel('Frame #');
ylabel('Deadline offset (ms)')
box off;
set(gca, 'TickDir', 'out');

saveas(gcf, [figure_path '/' figure_name '.pdf']);
close(fg);

% Calculate timing error
timing_error = FlipTimeStamp - VBL_Timestamp;

% Plot flip execution duration
figure_name = 'Fixation flip execution duration';
fg = figure('Color', [1 1 1], 'Name', figure_name);
plot(timing_error*(10^3)); % converting from s to ms

% Format figure
xlabel('Frame #');
ylabel('Flip exe duration (ms)')
box off;
set(gca, 'TickDir', 'out');

saveas(gcf, [figure_path '/' figure_name '.pdf']);
close(fg);

%% Create Grating

% Define grating
grating_size_deg = 5; % ie diameter 
grating_size_px = round(grating_size_deg * w.ppd);
grating_sf = 2 * grating_size_px/w.ppd;
grating_orientation = 90;
grating_phase = 0;
grating_contrast = 1;

% Create grating
[x,y] = meshgrid(-grating_size_px/2:grating_size_px/2-1, -grating_size_px/2:grating_size_px/2-1);
grating_tex = sin(grating_sf*2*pi / grating_size_px * (x.*sin(grating_orientation * (pi/180)) + y.*cos(grating_orientation * (pi/180))) - grating_phase);

% Convert grating to grayscale pixel values
grating_tex = grating_tex * grating_contrast * 127 + 127;
% figure, imshow(grating_tex(:,:,1),[0 255]);

grating = Screen('MakeTexture', window, grating_tex);
grating_patch = CenterRectOnPoint([0 0 grating_size_px grating_size_px], centerX, centerY);

%% Create circular aperture

aperture_size_deg = grating_size_deg * 1.5;
aperture_size_px = round(aperture_size_deg * w.ppd);

aperture_radius_px = round((grating_size_deg-1)/2 * w.ppd);

[x,y] = meshgrid(-aperture_size_px/2:aperture_size_px/2-1, -aperture_size_px/2:aperture_size_px/2-1);
r = sqrt(x.^2 + y.^2);

aperture_mask = ones(aperture_size_px);
aperture_mask(r <= aperture_radius_px) = 0;

edge_smooth_sd = 0.1 * w.ppd;
aperture_mask = imgaussfilt(aperture_mask, edge_smooth_sd);

aperture_tex(:,:,1) = ones(size(aperture_mask)) * 127;
aperture_tex(:,:,2) = aperture_mask * 255;
% figure, imshow(aperture_tex(:,:,2),[0 255]);

aperture = Screen('MakeTexture', window, aperture_tex);
aperture_patch = CenterRectOnPoint([0 0 aperture_size_px aperture_size_px], centerX, centerY); 

%% Make grating

% Define orientation change
num_orientations = 180;
orientations = linspace(0,359,num_orientations);

% Set timing
frame_rate = 20; % Hz
frame_dur = 1/frame_rate; % s
grating_dur = 5; % s
orientation_dur = grating_dur / num_orientations;

% Define frame onsets
frame_onsets_baseline = 0:frame_dur:fixation_dur-frame_dur;
num_frames = length(frame_onsets_baseline);

% Define fixation color per frame
num_frames_per_orientation = round(orientation_dur/frame_dur);

orientation_update = zeros(1,num_frames);
orientation_update(num_frames_per_orientation+1:num_frames_per_orientation:num_frames) = 1;

%% Draw Grating

curr_orientation = 1;
for n_frame = 1:num_frames

    if n_frame == 1
        frame_onsets = frame_onsets_baseline + GetSecs;
    end
    
    % Change fixation color when it's time
    if orientation_update(n_frame)
        curr_orientation = curr_orientation + 1;
    end

    % Draw grating
    Screen('DrawTexture', window, grating, [], grating_patch, orientations(curr_orientation));

    % Draw aperture
    Screen('DrawTexture', window, aperture, [], aperture_patch);

    % Fixation dot
    Screen('FillOval', window, black, fixation_dot_patch);

    % Flip
    [VBL_Timestamp(n_frame), StimulusOnsetTime(n_frame), FlipTimeStamp(n_frame), Missed(n_frame)] = Screen('Flip', window, frame_onsets(n_frame));

end


%% Evaluate grating timing
% VBLTimestamp': a high-precision estimate of the system time (in seconds) when the actual flip has happened
% 'StimulusOnsetTime': An estimate of Stimulus-onset time
% 'FlipTimestamp' is a timestamp taken at the end of Flip's execution. 
% Use the difference between FlipTimestamp and VBLTimestamp to get an estimate of how long Flips execution takes.
 
% 'Missed' indicates if the requested presentation deadline for your stimulus has been missed. 
% A negative value means that deadlines have been satisfied. 
% Positive values indicate a deadline-miss.

% Set figure path
figure_path = 'TimingLogs';
if ~exist(figure_path,'dir')
    sca;
    error([figure_path ' does not exist!']);
end

% Summarize frames missed
frames_missed = sum(Missed > 0);

% Plot frames missed
figure_name = 'Grating missed frames';
fg = figure('Color', [1 1 1], 'Name', figure_name);
plot(Missed * 10^3);
hold on
line([0 max(xlim)], [0 0], 'LineStyle','--', 'Color', black)

% Format figure
title(['Frames missed: ' num2str(frames_missed) ' (' num2str(100*frames_missed/num_frames) '%)']);
xlabel('Frame #');
ylabel('Deadline offset (ms)')
box off;
set(gca, 'TickDir', 'out');

saveas(gcf, [figure_path '/' figure_name '.pdf']);
close(fg);

% Calculate timing error
timing_error = FlipTimeStamp - VBL_Timestamp;

% Plot flip execution duration
figure_name = 'Grating flip execution duration';
fg = figure('Color', [1 1 1], 'Name', figure_name);
plot(timing_error*(10^3)); % converting from s to ms

% Format figure
xlabel('Frame #');
ylabel('Flip exe duration (ms)')
box off;
set(gca, 'TickDir', 'out');

saveas(gcf, [figure_path '/' figure_name '.pdf']);
close(fg);

%% Close window

Screen('LoadNormalizedGammaTable', window, default_gamma);
Screen('CloseAll');
ShowCursor;