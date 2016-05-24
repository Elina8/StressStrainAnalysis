# StressStrainAnalysis

% Created by Evangeline Tzatzalos
% Most recently updated on 05/24/2016
% Copyright (C) Evangeline Tzatzalos
% For more information contact etzatzalos@gmail.com
%
%%The goal of this script is to 
%   (1) Load the data of Force Generation for each EHM
%   (2) Plot the overal data for each ring (from one file). Data is in
%       seconds and mN.
%   (3) Plot Twitch Tensile (TT), Resting Tensile (RT), and Active Forces
%       (AF) as a function of time
%
% INPUT FILES
%   (1) Input data file is a *.csv file.
%   (2) The first worksheet is labeled "DATA". This contains the time and
%       force vectors for each sample.  There is no limit on the number of
%       samples that can be included in this worksheet.
%   (3) The second worksheet is labeled "STEPS". This worksheet contains
%       the location of each step (not time).  This information must be
%       manually calculated since the steps are manually performed.
%
% INPUT PARAMETERS
% n                                % INPUT; total # samples
% sampleNum                        % INPUT; selected sample for analysis
% lengthEHM                        % INPUT; mm of length of EHM
% diameterEHM                      % INPUT; mm of diameter of EHM
% DiffMin                          % default = 25, AF calculation: for pairing maxes (twitches) and mins (passives) for substration 
%                                   (AF)
% DiffMax                          % default = 30,
% adjust                           % refer below FIND PEAKS
% MPD                              % refer below FIND PEAKS
% linearRegion                     % refer below ELASTIC MODULUS, (first location of linear region)/(last location of linear region)
% max_y_forstressgraph             % the maximum y-value on the graph for Force vs Time, Stress vs Strain (RT-TT-AF)
% strain_xlimit                    % the maximum x-value on graph for strain
%  
% OUTPUT
%   (1) Elastic modulus is a numeric output.
%   (2) Five graphs are provided in the output
%       (a) Figure 1. Distance (micrometer) vs. Time (seconds)
%       (b) Figure 2. Force (mN) vs. Time (seconds)
%       (c) Figure 3. Stress (mN/mm^2) vs. Strain 
%       (d) Figure 4. Tension (Active, Resting, Twitch (mPa)) vs. Time (seconds)
%       (e) Figure 5. Active Stress (mPa) vs. Time (seconds)
%
% Notes: 
% No stimulation, 1Hz, 2Hz (20sec intervals)
% 4/1/2014 automated the generation of the strain plot
