% This is code by Zhipeng Cao and Will Grissom who have given permission 
% for inclusion within this package. This code can be found also in their
% repo https://bitbucket.org/wgrissom/acptx/ Please cite appropriately.

function [all_images, waveforms, nSVT] = dzktpts(algp,prbp,maps,kinit)

fov = maps.fov;  % Field of View in each dim, cm
dimxyz = prbp.dimxyz; % pixels, dim of design grid
Ns = prod(dimxyz); % total # pixels
delta_tip = prbp.delta_tip; % flip angle, degrees

Np = prbp.Npulse; % Number of (hard) subpulses
Nsubpts = prbp.Nsubpts;
nblippts = prbp.nblippts;
trajres = prbp.trajres;

ndim = ndims(maps.mask);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load and interpolate b1 and b0 maps onto design grid
Nc = size(maps.b1,4); % # tx channels

% load B1 maps
B1in = repmat(maps.mask,[1 1 1 Nc]).*maps.b1;

% in case its not already done, subtract off first coil's phase
sens = B1in.*exp(-1i*repmat(angle(B1in(:,:,:,1)),[1 1 1 Nc]));

mask = maps.mask;
fmap = maps.b0 .* mask;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% design a slice-selective pulse
% on this gradient to init basis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b1d = [ones(Nsubpts,1);zeros(nblippts,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get desired pattern for design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ndim == 3
    [xx,yy,zz]=ndgrid(-fov(1)/2:fov(1)/dimxyz(1):fov(1)/2-fov(1)/dimxyz(1), ...
        -fov(2)/2:fov(2)/dimxyz(2):fov(2)/2-fov(2)/dimxyz(2), ...
        -fov(3)/2:fov(3)/dimxyz(3):fov(3)/2-fov(3)/dimxyz(3));
elseif ndim == 2
    [xx,yy]=ndgrid(-fov(1)/2:fov(1)/dimxyz(1):fov(1)/2-fov(1)/dimxyz(1), ...
        -fov(2)/2:fov(2)/dimxyz(2):fov(2)/2-fov(2)/dimxyz(2));
else
    disp('Error');
    return;
end;

% desired excitation pattern
d = ones(Ns,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mask b1 maps, target pattern 
% for design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fmapd = fmap(logical(mask(:)));
sensd = reshape(sens,[Ns Nc]);
sensd = sensd(logical(mask(:)),:);
dd = d(logical(mask(:)),:);
xxd = xx(logical(mask(:)));
yyd = yy(logical(mask(:)));
if ndim == 3
    zzd = zz(logical(mask(:)));
    xyzd = [xxd yyd zzd];
elseif ndim == 2
    xyzd = [xxd yyd];
else
    disp('Error')
    return    
end
if isfield(maps,'phsinit')
  prbp.phsinit = maps.phsinit(logical(mask(:)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   design small-tip pulses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% OMP pre-definitions
if ndim == 3
    dimsearch = 24; kfovsearch = 1/trajres;
    [kxs,kys,kzs] = ndgrid(-kfovsearch/2:kfovsearch/dimsearch:kfovsearch/2-kfovsearch/dimsearch);
    kxs = kxs(:);kys = kys(:);kzs = kzs(:);
    [~,zi] = find(kxs == 0 & kys == 0 & kzs == 0);
    ks = [kxs([1:zi-1 zi+1:end]) kys([1:zi-1 zi+1:end]) kzs([1:zi-1 zi+1:end])];
else
    dimsearch = 32; kfovsearch = 1/trajres;
    [kxs,kys] = ndgrid(-kfovsearch/2:kfovsearch/dimsearch:kfovsearch/2-kfovsearch/dimsearch);
    kxs = kxs(:);kys = kys(:);
    zi = dimsearch*(dimsearch/2) + dimsearch/2 + 1;
    kxs = [kxs(1:zi-1);kxs(zi+1:end)]; % remove DC as OMP candidate
    kys = [kys(1:zi-1);kys(zi+1:end)]; % remove DC as OMP candidate
    ks = [kxs kys];
end;

% run the design
[rf,kpe,mout,compWts,nSVT] = localKtPts(dd*delta_tip*pi/180, ...
    sensd,xyzd,ks,fmapd,b1d,algp,prbp);
waveforms.compWts = compWts;

% embed resulting m into full array
m = zeros(size(mask));
m(logical(mask(:))) = mout;
rf = reshape(rf,[Np Nc]);

all_images.images = m;
waveforms.rf = rf;
waveforms.k = kpe;


