% This test script loads invivo B1+ datasets of the human body at 7T 
% and computes frequency robust universal kT-point pulses.
%
% Christoph S. Aigner, Sebastian Dietrich, Felix Krüger, Max Lutz and Sebastian
% Schmitter, Towards frequency robust tailored and universal pulses in the 
% human heart at 7T, submitted to ISMRM 2022
%
% The 36 channel-wise invivo B1+ datasets of the human body at 7T are 
% available at: https://doi.org/10.6084/m9.figshare.17074724.v1
%
% The optimization of the kT-points is performed using code by Will Grissom
% and Zhipeng Cao (https://bitbucket.org/wgrissom/acptx/) who have given 
% permission for inclusion within this package. Please cite appropriately.
% 
% Created by Christoph S. Aigner, PTB, November 2021.
% Email: christoph.aigner@ptb.de
%
% This code is free under the terms of the GPL-3.0 license.

addpath("ktutil")
load("cmap.mat");
load("kTrandphases.mat");

rangeofkTpoints = 1:8;
frequ = [ -1129, -1010, -772, -576, -115, 0, 178];             %in Hz

frequ=0;
for c_kTnumb=rangeofkTpoints
    prbp.Npulse = c_kTnumb                  % number of kT-points subpulses
    prbp.pathDat        = 'B1R'; % set the folder that contains the in vivo B1+ datasets
    prbp.Nc = 8;
    c_diffrand= 165;
    c_dat = 1; % set dataset counter to 1
    prbp.allIndices = 1:31;
    disp(['load ', num2str(length(prbp.allIndices)), ' invivo B1+ maps']);
    count_maps=0;
    count_maps2=0;
    allmaps = cell(1,2);
    mapsUP=[];

    for c_subj=prbp.allIndices
        %load the Matlab container for dataset #countsubj
        if ~isfile([prbp.pathDat '/lightB1R_' num2str(prbp.allIndices(c_subj)) '.mat'])
            disp('NO B1 MAPS FOUND!');
            return;
        end
        load([prbp.pathDat '/lightB1R_' num2str(prbp.allIndices(c_subj)) '.mat']);
        B1R = non_respiration_resolved_B1R;

        %rearrange the dimensions 
        B1ptemp.cxmap = permute(squeeze(B1R.B1Rp),[3 2 1 4]);
        B1ptemp.cxmap = B1ptemp.cxmap(end:-1:1,:,end:-1:1,:);
                
        %normalization wrt to the mean in the heart ROI
        sumabsB1 = sum(abs(B1ptemp.cxmap),4);
        meansumabsB1mask = mean(sumabsB1(B1R.kTpoints.maps.mask));
        maps.b1 = B1ptemp.cxmap/meansumabsB1mask;
    
        if c_subj == 1
            B1Dim = size(maps.b1);
            Nc = B1Dim(4); %number of transmit channels
            fov = B1R.kTpoints.maps.fov; %fov in cm
        else
            if B1Dim ~= size(maps.b1)
                disp(['loaded B1 map #',num2str(c_subj), ...
                      ' has a different size as B1 map #1 !!!']);
            end
        end

        %save the other parameters. So far, B0 maps are set to 0
        maps.b0      = maps.b1*0;
        maps.mask    = B1R.kTpoints.maps.mask;
        maps.fov     = B1R.kTpoints.maps.fov;
        maps.DatNum  = prbp.allIndices(c_subj);
        B1p.cxmap = maps.b1;
        ROI.masks = maps.mask;

        %% create the variables to simulate the optimized pulse with different freq
        %general Parameters
        maps.numberofmaps = length(frequ);  % two different B1 maps with different frequencies 
        maps.fov = [312.5000  312.5000  250.0000]/10; %fits to the 3D RPEB1R (sag) 
        Nc = size(B1p.cxmap,4); % number of Tx coils (Ncoils)
        Nz = size(B1p.cxmap,3); % number of Tx coils (Ncoils)
        phsinitmode = 'randphase';

        %get the B1 maps out of B1p.cxmap (TBD: normalize the B1 maps to the median?) 
        maps.b1 = double(B1p.cxmap); %be sure to have doubles!
        maps.b1(isnan(maps.b1)) = 0;     %set potential NANs to 0

        %get the ROI out of ROI.mask
        maps.mask = logical(ROI.masks); %be sure to have logicals!

        %get the B0 map (We dont use B0 here)
        maps.b0=maps.mask*0; %so far use the mask for B0

        %flip the B1 map, B0 map and the ROI
        maps.b1 = maps.b1(end:-1:1,end:-1:1,:,:);
        maps.b0 = maps.b0(end:-1:1,end:-1:1,:);
        maps.mask = maps.mask(end:-1:1,end:-1:1,:);

        %initialize the maps (store the original maps)
        for count_freq=1:maps.numberofmaps
            count_maps=count_maps+1;
            allmaps{count_maps}.b1 = maps.b1;
            allmaps{count_maps}.mask = maps.mask;
            allmaps{count_maps}.b0 = maps.b0;
        end

        %build one set of maps and modify the B0 offset
        for count_freq=1:maps.numberofmaps
            count_maps2=count_maps2+1;
            mapsUP.b1(:,:,(1:Nz)+(count_maps2-1)*Nz,:)=allmaps{count_maps2}.b1;
            mapsUP.mask(:,:,(1:Nz)+(count_maps2-1)*Nz)=allmaps{count_maps2}.mask;
            mapsUP.b0(:,:,(1:Nz)+(count_maps2-1)*Nz)=allmaps{count_maps2}.b0+frequ(count_freq);
        end

    end
    maps=mapsUP;
    maps.numberofmaps = count_maps2;  % two different B1 maps with different frequencies 
    maps.fov = [312.5000  312.5000  250.0000]/10; %fits to the 3D RPEB1R (sag) 
    prbp.Ncred = inf;
    
    % set initial target phase                 
    switch phsinitmode
        case 'zerophase'
            disp('zero phase initial')
            maps.phsinit = zeros(size(maps.mask)); 
        case 'defaultphase'
            disp('default phase initial')
            maps.phsinit = angle(sum(maps.b1,4));%default phase 
        case 'quadmode' %quad mode does not perform in the body
            bcb1 = 0;
            for ii = 1:prbp.Nc
               bcb1 = bcb1 + maps.b1(:,:,:,ii)*...
                             exp(1i*(ii-1)*2*pi/prbp.Nc).*...
                             exp(-1i*angle(maps.b1(:,:,:,1)));
            end
            maps.phsinit = angle(bcb1);
        case 'randphase'
            bcb1 = 0;
            for ii = 1:prbp.Nc
                bcb1 = bcb1 + maps.b1(:,:,:,ii)*...
                             exp(1i*randphases(c_diffrand,ii)).*...
                             exp(-1i*angle(maps.b1(:,:,:,1)));
            end
            maps.phsinit = angle(bcb1);
    end

    % Algorithm and problem parameters
    prbp.delta_tip = 10;              % flip angle, degrees
    prbp.ndims = ndims(maps.mask);    % # spatial dimensions 
    prbp.kmaxdistance = [Inf Inf Inf];% maximum kT-point location
    prbp.beta = 500;                % initial RF regularization parameter
    prbp.betaadjust = 0;              % automatically adjust RF regularization parameter
    prbp.dt = 10e-6;                  % dwell time, sec
    prbp.Nsubpts = floor(0.1/prbp.dt/1000);     % # of time points in the subpulses (0.1 ms subpulses)
    prbp.nblippts = 20;               % number of time points between hard pulses to accomodate gradient blips
    prbp.dimxyz = size(maps.b0);
    prbp.filtertype = 'Lin phase';    % alternately add kT-points on either size of the DC point
    prbp.trajres = 2;                 % maximum spatial frequency of OMP search grid (was 2 for most August 2015 results)
    %prbp.Npulse = 4;                  % number of kT-points subpulses
    algp.nthreads  = 2;               % number of compute threads (for mex only)
    algp.computemethod = 'mex';
    algp.ncgiters = 3; 
    algp.cgtol = 0.9999;

    % Run the kT point design 
    %   m ... STA solution 
    %   wvfrms ... optimized RF and gradient blips)        
    [m, wvfrms] = dzktpts(algp,prbp,maps); 
    rfw = wvfrms.rf;

    %evaluate the optimized results
    farmse = sqrt(mean((abs(m.images(maps.mask))/pi*180 - prbp.delta_tip).^2));
    rfrms = norm(rfw);
    fprintf('Flip angle RMSE: %.4f, RMS RF power: %.4f.\n\n',farmse,rfrms);
    
    save(['UP_Frob_', num2str(length(frequ)),'freq_', num2str(c_kTnumb),'kT_4.mat'],'wvfrms');
end
