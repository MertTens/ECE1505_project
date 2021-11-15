% ========================================================================================
%
% NAME: 	readTRICKS_DCE.m
% PURPOSE: 	main program to import TRICKS DCE-MRI data 
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     July 15, 2009
% MODIFIED: July 15, 2009
%
% ========================================================================================

warning off;

% Matlab m-file directory
addpath_dir = 'C:\D''ocuments and Settings\''Ha''i-Ling Cheng\''M''y Documents\''q''MRI tools\''DCE\';
eval(['addpath ' addpath_dir]);
eval(['addpath ' addpath_dir '\Tofts']);
eval(['addpath ' addpath_dir '\Patlak']);
eval(['addpath ' addpath_dir '\extractAIF']);

% -------------------------------------------------------------------------
% User input: directory to save output data; DCE sequence parameters.
% -------------------------------------------------------------------------
numSlices = input('How many slices? ');
del_t = input('Time resolution (sec)? ');
numTimePts = input('How many time-points? ');
numBaseline = input('How many pre-contrast time-points? ');
OutputDir = uigetdir('C:\D''ocuments and Settings\''Ha''i-Ling Cheng\''M''y Documents\','Choose final directory to save data');
DCEDir = uigetdir('C:\D''ocuments and Settings\''Ha''i-Ling Cheng\''M''y Documents\','Choose DCE raw data directory');

% -------------------------------------------------------------------------
% Define constants.
% -------------------------------------------------------------------------
MODEL = 'AUC';  % analysis = {'Patlak','Tofts','Tofts-biexp','AUC'}
Dose = 0.1;             % Gd dose (mmol/kg) 
r1 = 4.1;               % T1 relaxivity (liter/mmol/sec) -- Magnevist
t = [0:del_t:(numTimePts-1)*del_t] / 60; % time vector of DCEMRI (min)

% -------------------------------------------------------------------------
% Import AIF (into variable Cp). Can be either measured or assumed. 
% NOTE: temporal resolution of AIF must be the same as that of DCE data.
% -------------------------------------------------------------------------
if (input('Is there an AIF matlab data file? [y/n] ','s')=='y')     
    [filename,PathName] = uigetfile('*.mat','Select the AIF matlab data file');
    cd (PathName); load (filename);
end

switch MODEL     
     case {'Tofts-biexp'} % Tofts with bolus delay (fit biexp to AIF)        
         A=[3.99 4.78];        % initial guess for amplitudes (kg/liter)
         m=[0.144 0.0111];    % initial guess for decay constants (1/min)
         x0=[A m]; lb=[0.1 0.1 0.001 0.001]; ub=[50 50 20 20];
         [x,resnorm] = fit_biexpCp(t,Cp,x0,lb,ub,options,Dose);          
end
        
% -------------------------------------------------------------------------
% Import T1maps (into variable T1).
% -------------------------------------------------------------------------
[filename,PathName] = uigetfile('*.mat','Select the T1map matlab data file');
cd (PathName); load (filename);

% -------------------------------------------------------------------------
% Read in TRICKS DCE data.
% -------------------------------------------------------------------------
for m = 1:numTimePts
    if m < 10
        newDir=strrep(DCEDir,' - Tricks', ['0' num2str(m) ' - Tricks']);
    else
        newDir=strrep(DCEDir,' - Tricks', [num2str(m) ' - Tricks']);
    end
    
    for k = 1:numSlices
    filename = [newDir '\' num2str(k) '.dcm'];
    info = dicominfo(char(filename));
    A(:,:,k,m) = double(dicomread(info));    
    
    if m == 1
        SliceLoc(k) = info.SliceLocation;
        if k == 1
            TR = info.RepetitionTime;       % msec
            FA = info.FlipAngle*pi/180;     % radians
        end
    end
    end
end
clear info;

figure(1); [x,y]=getpts;
while (x<256 & x>0 & y>0 & y<256)
figure(2); plot(t,squeeze(A(round(y(1)),round(x(1)),4,:)));
figure(1); [x,y]=getpts;
end

% -------------------------------------------------------------------------
% Select ROI on DCE images
% -------------------------------------------------------------------------
if (input('Is there a T1w ROI file to import? [y/n] ','s')=='y')     
    [filename,PathName] = uigetfile('*.mat','Select the ROI matlab data file');
    cd (PathName); load (filename);
else
    for k = 1:numSlices
        figure(1); imagesc(A(:,:,k,1)); title(['Slice Position: ' num2str(SliceLoc(k))]);
        mask(:,:,k) = roipoly;
    end
end
mask = double(mask);

% -------------------------------------------------------------------------
% DCE-MRI analysis
% -------------------------------------------------------------------------
Ktrans_map = zeros(size(A,1),size(A,2),numSlices); 
vp_map = zeros(size(A,1),size(A,2),numSlices);
ve_map = zeros(size(A,1),size(A,2),numSlices);
Ktrans = 0; vp = 0; ve = 0;

for k = 1:numSlices
disp(['processing slice#' num2str(k)]);
for nx = 1:size(A,2)
for ny = 1:size(A,1)
    if mask(ny,nx,k)==1  
        
        % convert tissue signal into tissue concentration curve
        SI = squeeze(A(ny,nx,k,:)).';
        SIo = mean(SI(1:numBaseline));
        E1 = exp(-TR/T1(ny,nx,k));
        gain = SIo*(1-cos(FA)*E1)/(1-E1);
        R1 = log((gain-SI*cos(FA))./(gain-SI)) / TR;
        Ct = 1e3*(R1 - 1/T1(ny,nx,k)) / r1;            % (mmol/liter)
        
        switch MODEL
            case {'Patlak'}     % (AIF measured): extract Ktrans, vp  
                [Ktrans,vp,Ktrans_max,vp_max] = PTLK_DCEMRI(Cp,Ct,t);
                
            case {'Tofts'}      % Tofts with bolus delay (AIF measured): extract Ktrans, vp, ve, tau  
                [Ktrans,vp,ve,tau,chi_err] = mTK_bolus_delay_DCEMRI(Cp,Ct,t);
                
            case {'Tofts-biexp'} % Tofts with bolus delay (fitted biexp AIF): extract Ktrans, vp, ve, tau 
                [Ktrans,vp,ve,tau,chi_err] = mTK_biexpCp_DCEMRI([x,Dose],Ct,t);
            
            case {'AUC'}         % AUC analysis
                [AUC] = AUC_DCEMRI(Ct,t,interval);    
        end
      
        Ktrans_map(ny,nx,k) = Ktrans;
        vp_map(ny,nx,k) = vp;
        ve_map(ny,nx,k) = ve;
        AUC_map(ny,nx,k) = AUC;
    end
end
end
disp(['finished slice ' num2str(k)]);
end


% -------------------------------------------------------------------------
% Save output
% -------------------------------------------------------------------------
cd(DirName)
save DCEMaps Ktrans_map ve_map vp_map AUC_map


