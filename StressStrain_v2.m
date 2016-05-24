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
%% INPUT PARAMETERS
% n                         % INPUT; total # samples
% sampleNum                 % INPUT; selected sample for analysis
% lengthEHM                 % INPUT; mm of length of EHM
% diameterEHM               % INPUT; mm of diameter of EHM
% DiffMin                   % default=25, before 4, Below in AF calculation: for pairing maxes (twitches) and mins (passives) for substration (AF)
% DiffMax                   % default=30, before 15
% adjust                    % default=10, before 10; refer below FIND PEAKS
% MPD                       % refer below FIND PEAKS, how often you expect to find peaks.  we expect to find 120 peaks per second.
% lowerForce                % default=0, Fig 2
% upperForce                % default=6, Fig 2
% lowerStress               % default=0
% upperStress               % default=3, the maximum y-value on the graph for Force vs Time, Stress vs Strain (RT-TT-AF)
% strain_xlimit             % defaut = 0.21; the maximum x-value on graph for strain
% maxY_distanceVsTime       % default=2
% lowerTime                 % default=0, for plotting
% upperTime                 % default=250, for plotting
% y_AF                      % default=0.3,
% smoothingFactor           % default=10;
%  
%% OUTPUT
%   (1) Elastic modulus is a numeric output.
%   (2) Five graphs are provided in the output
%       (a) Figure 1. Distance (micrometer) vs. Time (seconds)
%       (b) Figure 2. Force (mN) vs. Time (seconds)
%       (c) Figure 3. Stress (mN/mm^2) vs. Strain 
%       (d) Figure 4. Tension (Active, Resting, Twitch (mPa)) vs. Time (seconds)
%       (e) Figure 5. Active Stress (mPa) vs. Time (seconds)
%%
% Additional Notes
% No stimulation, 1Hz, 2Hz (20sec intervals)
% 4/1/2014 automated the generation of the strain plot

close all; clear all; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LOAD FORCE DATA, DEFINE SAMPLE, AND DEFINE INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%INPUT
n=6;                                 % INPUT; total # samples
sampleNum=1;                         % INPUT; selected sample for analysis
force=2*sampleNum;                   % the column on 'compiled_data.xlsx' with raw force data
time=2*sampleNum-1;                  % the column on 'compiled_data.xlsx' with time (in seconds) values corresponding to raw force data
lengthEHM=8;                        % INPUT; mm of length of EHM
diameterEHM=1;                       % INPUT; mm of diameter of EHM
DiffMin=15; %default = 25, before 5   % Below in AF calculation: for pairing maxes (twitches) and mins (passives) for substration (AF)
DiffMax=30; %default = 30, before 15
adjust=10; %before 10               % refer below FIND PEAKS
MPD=0.8;%1/(60*2);                             % refer below FIND PEAKS, how often you expect to find peaks.  we expect to find 120 peaks per second.
lowerForce=0;                       % Fig 2
upperForce=6;                       % Fig 2
lowerStress=0;
upperStress=3;                     % the maximum y-value on the graph for Force vs Time, Stress vs Strain (RT-TT-AF)
strain_xlimit=0.21;                  % the maximum x-value on graph for strain
maxY_distanceVsTime=2;
lowerTime=0;                        % for plotting
upperTime=250;                      % for plotting
y_AF=0.3;
smoothingFactor=10;

%READ-IN DATA
[EHMall,header]=xlsread('2015_0812_rEHM003.xlsx','DATA');
EHM=EHMall(:,time:force);  %Redefine EHM to only look at time and force columns of the sample you are running
EHM=EHM(~any(isnan(EHM),2),:);
EHM(:,2)=smooth(EHM(:,2),smoothingFactor);
[stepIndicesAll,header2]=xlsread('2015_0812_rEHM003.xlsx','STEPS'); % the file with the compiled raw data and stepIndices.  The stepIndices were manually determined from the max peak of force generated at each step of stretch
stepIndices=stepIndicesAll(:,sampleNum+1); % create a variable to pull out stepsIndices related to the sample
stepIndices=stepIndices(~isnan(stepIndices)); %remove NaN from list
stepIndices=stepIndices(1:(end-2)); % stepIndices not including lower and upper limits of the linear region
steps=length(stepIndices);          % the total number of steps of stretch that are being measured
linearRegion=stepIndicesAll(end-1,sampleNum+1)/stepIndicesAll(end,sampleNum+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT STRETCH REGIMEN (Distance (mm) vs. Time (sec))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                              %written 4/1/2014, do not block out, need it for further sections
distance=zeros(stepIndices(end),1);                           % create a vector of zeros which will hold the new values of distance (stepIndices will get converted to distance of stretch)
for i=1:steps-1
    distance(stepIndices(i):stepIndices(i+1),1)=i*0.125;
end
                                                              % distance that the EHM was stretched in mm
figure; % Figure 1
plot(EHM(stepIndices(1):stepIndices(steps),1),distance(stepIndices(1):stepIndices(steps)), 'Color',[0 0 0]);
     title(header(1,time),'fontsize',16,'fontweight','bold');
     xlabel(header(2,time),'fontsize',14,'fontweight','bold'); 
     set(gca,'FontSize',14);
     ylabel('Distance of Stretch (mm)','fontsize',14,'fontweight','bold');
     ylim([0 maxY_distanceVsTime]);

     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%PLOT RAW FORCE DATA (Force (mN) vs. Time (sec))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; % Figure 2
plot(EHM(stepIndices(1):stepIndices(steps),1),EHM(stepIndices(1):stepIndices(steps),2),'color',[0 0 0])   
     title(header(1,time),'fontsize',16,'fontweight','bold');
     xlabel(header(2,time),'fontsize',14,'fontweight','bold'); 
     set(gca,'FontSize',14);
     ylabel(header(2,force),'fontsize',14,'fontweight','bold');
     ylim([lowerForce upperForce]);
     xlim([lowerTime upperTime]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT STRESS (mN/mm^2) vs. STRAIN(deltaL/L0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stress (mPa) and Strain 

ehmStrain=distance/lengthEHM;
ehmStress=EHM(1:stepIndices(steps),2)/(2*pi*(diameterEHM/2)^2); 

figure; % Figure 3
plot(ehmStrain,ehmStress,'k.')   
     title('Force-Length Dependence: Stress vs Strain','fontsize',16,'fontweight','bold')
     xlabel('% Strain ((L-Lo)/Lo)','fontsize',14,'fontweight','bold')
     xlim([0 strain_xlimit]);
     ylim([lowerStress upperStress]);
     set(gca,'FontSize',14)
     ylabel('Stress (mN/mm^2)','fontsize',14,'fontweight','bold')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIND PEAKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
partialEHM=[ehmStrain ehmStress]; % Call Stress Data
% figure; % Figure 4
% plot(partialEHM(:,1),partialEHM(:,2),'r');

[MAXpeaks,MAXlocs]=findpeaks(partialEHM(:,2),'minPeakDistance', MPD); %MPD=17;  % INPUT (set above in 1st section) sets the minPeakDistance     <--- CHANGE to fit MAX and MIN peaks

temp_EHM= -1.*partialEHM(:,2);% Flip force data to find minimums
[MINpeaks,MINlocs]=findpeaks(temp_EHM(:,1),'minPeakDistance', MPD);

figure; hold on; % Figure 2b
plot(EHM(1:stepIndices(steps),1),ehmStress(1:stepIndices(steps),1),'color',[0 0 0])   
     title(header(1,time),'fontsize',16,'fontweight','bold');
     xlabel(header(2,time),'fontsize',14,'fontweight','bold'); 
     set(gca,'FontSize',14);
     ylabel('Stress kPa','fontsize',14,'fontweight','bold');
     ylim([lowerStress upperStress]);
     xlim([lowerTime upperTime]);
plot(EHM(MAXlocs,1),MAXpeaks,'r*');   % Maxpeaks
plot(EHM(MINlocs,1),-MINpeaks,'b*');   % Minpeaks
hold off;
     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ACTIVE FORCE = TWITCH(MAX) - RESTING(MIN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AF=zeros(length(MAXlocs),1);
AF_time=zeros(length(MAXlocs),1);

% for i=1:length(MAXlocs)
%     for j=1:length(MINlocs)
%         if (MINlocs(j)>MAXlocs(i)) && (MINlocs(j)-MAXlocs(i)>DiffMin) && (MINlocs(j)-MAXlocs(i)<=DiffMax);
%             AF(i)=partialEHM(MAXlocs(i),2)-partialEHM(MINlocs(j),2);
%             AF_time(i)=partialEHM(MAXlocs(i),1);
%             break
%         end
%     end
% end

for i=1:length(MAXlocs)
    for j=1:length(MINlocs)
        if (MINlocs(j)>MAXlocs(i)) && (MINlocs(j)-MAXlocs(i)>DiffMin) && (MINlocs(j)-MAXlocs(i)<=DiffMax);
            AF(i)=MAXpeaks(i)+MINpeaks(j);
            AF_time(i)=partialEHM(MAXlocs(i),1);
            break
        end
    end
end


figure; hold on; % Figure 5
plot(AF_time,AF,'r.');
plot(partialEHM(MAXlocs,1),MAXpeaks,'g.');   % Total tension
plot(partialEHM(MINlocs,1),-MINpeaks,'b.');   % Passive tension
   title('Total, Passive, and Active Stresses', 'fontSize',16,'fontWeight','bold');
   xlabel('Strain','fontSize',14,'fontWeight','bold');
   ylabel('Stress (kPa)','fontsize',14,'fontweight','bold');
   legend('Passive Stress');
   legend('Active Stress', 'Total Stress', 'Passive Stress');
   xlim([0 strain_xlimit]);
   ylim([lowerStress upperStress]);
   set(gca,'fontSize',14,'fontWeight','bold');
hold off;

% plot only active stress, Figure 6
figure; hold on;
plot(AF_time,AF,'r.')
   title('Active Stress', 'fontSize',16,'fontWeight','bold');
   xlabel('Strain','fontSize',14,'fontWeight','bold');
   ylabel('Stress (kPa)','fontsize',14,'fontweight','bold');
   xlim([0 strain_xlimit]);
   ylim([0 y_AF]);
   set(gca,'fontSize',14,'fontWeight','bold');

sortAF=sort(AF,'descend');
maxAF=max(sortAF);
maxAF_top10=sortAF(1:10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Elastic Modulus
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Tensile tension / Strains at those stresses

tensileStrain=partialEHM(MAXlocs(1:length(MAXlocs)*linearRegion),1);
tensileStress=MAXpeaks(1:length(MAXpeaks)*linearRegion);
Y=polyfit(tensileStrain,tensileStress,1)   % Tensile tension
