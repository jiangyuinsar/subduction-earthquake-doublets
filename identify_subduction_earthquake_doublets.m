clc, clear all, close all

% Add all third-party scripts and helper functions to the MATLAB path.
addpath(genpath('ThirdPartyScript'))

% Create an output folder for saving the detected doublet lists.
mkdir('txt_Doublets')

% Input matrix: output_matrix = [lon lat z Mw YYYY MM DD Type]
% lon  = longitude
% lat  = latitude
% z    = depth
% Mw   = moment magnitude
% YYYY = event year
% MM   = event month
% DD   = event day
% Type = event type or mechanism classification (1: Thrust, 2: Strike-Slip, 3: Normal)
load('gCMT_catalog_19760101_20260101_M6+.mat','earthquake_catalog');

% Avoid exactly zero latitude values, which can occasionally cause numerical
% issues in coordinate-conversion functions.
earthquake_catalog(earthquake_catalog(:,2)==0,2) = earthquake_catalog(earthquake_catalog(:,2)==0,2) + 1e-10;

% Convert YYYY, MM, DD to MATLAB datetime format for time-difference calculation.
dates = datetime(earthquake_catalog(:, 5), earthquake_catalog(:, 6), earthquake_catalog(:, 7));

% Convert datetime values to serial date numbers.
% This makes it easy to calculate event-to-event time differences in days.
serial_dates = datenum(dates);

% Number of earthquakes in the catalog.
[n, ~] = size(earthquake_catalog);

% -------------------------------------------------------------------------
% Step 1: Pre-compute pairwise magnitude, time, and distance differences
% -------------------------------------------------------------------------
% For every pair of events, calculate:
%   mag_diff  = absolute magnitude difference
%   time_diff = time difference in days
%   dist2_km  = 2-D epicentral distance in km
%   dist3_km  = 3-D hypocentral distance in km, including depth difference
%
% Note:
% time_diff(i,j) = serial_dates(i) - serial_dates(j)
% Therefore, time_diff(i,j) > 0 means event i occurred after event j.
% -------------------------------------------------------------------------

for i = 1:n
    for j = 1:n

        % Magnitude difference between event i and event j.
        mag_diff(i,j) = abs(earthquake_catalog(i, 4) - earthquake_catalog(j, 4));

        % Time difference between event i and event j, in days.
        time_diff(i,j) = serial_dates(i) - serial_dates(j);

        % Horizontal distance between event i and event j.
        % llh2local converts longitude/latitude coordinates into local
        % Cartesian coordinates relative to the second event.
        xy = llh2local([earthquake_catalog(i, 1),earthquake_catalog(i, 2)]',[earthquake_catalog(j, 1),earthquake_catalog(j, 2)]);

        % Depth difference between event i and event j.
        dz = earthquake_catalog(i, 3) - earthquake_catalog(j, 3);

        % 2-D distance: horizontal epicentral distance only.
        dist2_km(i,j) = sqrt(sum(xy.^2));

        % 3-D distance: horizontal distance plus depth difference.
        dist3_km(i,j) = sqrt(sum(xy.^2) + dz^2);

    end
end

% -------------------------------------------------------------------------
% Step 2: Define doublet-search parameter combinations
% -------------------------------------------------------------------------
% The code tests multiple doublet definitions using different:
%   - maximum magnitude differences
%   - maximum spatial distances
%   - maximum time separations
%   - 2-D or 3-D distance calculations
% -------------------------------------------------------------------------

% Maximum allowed magnitude difference between the two events.
max_mag_diff_list = [0.25 0.50];

% Maximum allowed distance between the two events, in km.
max_dist_km_list = [50 100 150];

% Maximum allowed time difference between the two events, in days.
% Here the tested windows are 0.5 years, 1 year, and 10 years.
max_time_days_list = 365*[0.5 1 10];

% Distance definition:
%   2 = use 2-D epicentral distance
%   3 = use 3-D hypocentral distance
dimension_list = [2 3];

% -------------------------------------------------------------------------
% Step 3: Loop over all parameter combinations
% -------------------------------------------------------------------------
% For each combination of magnitude, distance, time, and dimension criteria,
% identify candidate doublet pairs and save the results.
% -------------------------------------------------------------------------

for i_mag = 1:length(max_mag_diff_list)
    for i_dist = 1:length(max_dist_km_list)
        for i_time = 1:length(max_time_days_list)
            for i_dimension = 1:length(dimension_list)
                tic

                % Select the current doublet-search thresholds.
                max_mag_diff = max_mag_diff_list(i_mag);
                max_dist_km = max_dist_km_list(i_dist);
                max_time_days = max_time_days_list(i_time);
                dimension = dimension_list(i_dimension);

                % Choose either 2-D or 3-D distance for the current run.
                if dimension == 2
                    dist_km = dist2_km;
                elseif dimension == 3
                    dist_km = dist3_km;
                end

                % Store detected doublet pairs.
                % Each row will contain:
                % [index1 index2 mag_diff time_diff dist_km mean_Mw Type1 Type2]
                doublets = [];

                % -----------------------------------------------------------------
                % Step 4: Search all event pairs for candidate doublets
                % -----------------------------------------------------------------
                % The loop only checks j = i+1:n to avoid duplicate pairs and
                % self-pairs.
                % -----------------------------------------------------------------
                for i = 1:n
                    for j = i+1:n

                        % First-level doublet criteria:
                        %   1. magnitude difference is small enough
                        %   2. spatial distance is small enough
                        %   3. absolute time separation is short enough
                        if mag_diff(i,j) <= max_mag_diff && dist_km(i,j) <= max_dist_km && abs(time_diff(i,j)) <= max_time_days

                            % -----------------------------------------------------
                            % Step 5a: Check whether event i has a larger nearby
                            % preceding earthquake.
                            % -----------------------------------------------------
                            % Criterion 1:
                            % A possible preceding event occurred before event i
                            % and within 180 days.
                            %
                            % Because time_diff(i,k) = time_i - time_k,
                            % time_diff(i,k) > 0 means event k occurred before i.
                            cond1 = time_diff(i, :) > 0 & time_diff(i, :) <= 180;

                            % Criterion 2:
                            % The possible preceding event is within 50 km of event i.
                            cond2 = dist_km(i, :) < 50;

                            % Criterion 3:
                            % The possible preceding event is larger than event i.
                            cond3 = earthquake_catalog(:, 4)' > earthquake_catalog(i, 4);

                            % Events satisfying all three criteria are larger,
                            % nearby, preceding events relative to event i.
                            all_criteria = cond1 & cond2 & cond3;

                            % Exclude event j from the list, because j is the paired
                            % candidate event itself.
                            matching_indices_for_event_i = setdiff(find(all_criteria),j);

                            % -----------------------------------------------------
                            % Step 5b: Check whether event j has a larger nearby
                            % preceding earthquake.
                            % -----------------------------------------------------
                            % Criterion 1:
                            % A possible preceding event occurred before event j
                            % and within 180 days.
                            cond1 = time_diff(j, :) > 0 & time_diff(j, :) <= 180;

                            % Criterion 2:
                            % The possible preceding event is within 50 km of event j.
                            cond2 = dist_km(j, :) < 50;

                            % Criterion 3:
                            % The possible preceding event is larger than event j.
                            cond3 = earthquake_catalog(:, 4)' > earthquake_catalog(j, 4);

                            % Events satisfying all three criteria are larger,
                            % nearby, preceding events relative to event j.
                            all_criteria = cond1 & cond2 & cond3;

                            % Exclude event i from the list, because i is the paired
                            % candidate event itself.
                            matching_indices_for_event_j = setdiff(find(all_criteria),i);

                            % -----------------------------------------------------
                            % Step 5c: Retain or reject the candidate pair
                            % -----------------------------------------------------
                            % The pair is retained if at least one of the two events
                            % does not have a larger nearby preceding event within
                            % 180 days and 50 km, excluding the paired event itself.
                            %
                            % In other words, the candidate pair is rejected only when
                            % both events have larger nearby preceding events.
                            if isempty(matching_indices_for_event_i)==1 | isempty(matching_indices_for_event_j)==1

                                % Save the doublet information:
                                %   column 1: index of event i
                                %   column 2: index of event j
                                %   column 3: magnitude difference
                                %   column 4: time difference in days
                                %   column 5: distance in km
                                %   column 6: mean magnitude of the two events
                                %   column 7: event type of event i
                                %   column 8: event type of event j
                                doublets = [doublets; i j mag_diff(i,j) time_diff(i,j) dist_km(i,j) mean([earthquake_catalog(i, 4), earthquake_catalog(j, 4)]) earthquake_catalog(i, 8) earthquake_catalog(j, 8)];
                            end
                        end

                    end
                end

                % -----------------------------------------------------------------
                % Step 6: Output results
                % -----------------------------------------------------------------

                % Sort doublets by the 6th column, mean magnitude, in descending order.
                doublets = sortrows(doublets, -6);

                % Save the detected doublet list for the current parameter set.
                writematrix(doublets,['txt_Doublets/Doublets_diffM' num2str(max_mag_diff) '_' num2str(max_dist_km) 'km_' num2str(max_time_days/365) 'year_' num2str(dimension) 'D_allmechanism.txt'],'Delimiter','tab')

                % Display the current parameter combination.
                disp([num2str(max_mag_diff) ' M, ' num2str(max_dist_km) ' km, ' num2str(max_time_days/365) ' years, ' num2str(dimension) '-D'])

                % Keep only the variables needed for the next parameter combination.
                clearvars -except max_mag_diff_list max_dist_km_list max_time_days_list dimension_list i_mag i_dist i_time i_dimension mag_diff time_diff dist2_km dist3_km earthquake_catalog n

                toc

            end
        end
    end
end