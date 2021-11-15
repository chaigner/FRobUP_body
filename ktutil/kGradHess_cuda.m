% This is code by Zhipeng Cao and Will Grissom who have given permission 
% for inclusion within this package. This code can be found also in their
% repo https://bitbucket.org/wgrissom/acptx/ Please cite appropriately.

function [vx,Sx,vy,Sy,vz,Sz] = kGradHess_cuda(A,d,xx,Nrmax)

persistent gpukGradHess
persistent vxg Sxg vyg Syg vzg Szg

[Ns,Nd] = size(xx);
Nr = size(A,2);

if isempty(gpukGradHess)

  disp 'Initializing GPU trajectory derivative calculations';
  % compile with: nvcc -ptx kGradHess.cu
  gpukGradHess = parallel.gpu.CUDAKernel('kGradHess.ptx','kGradHess.cu','kGradHess');
  tmp = gpuDevice;
  gpukGradHess.GridSize = [min(1024,tmp.MaxGridSize(1)) 1]; % min(1024, maxblocks)
  gpukGradHess.ThreadBlockSize = [min(128,tmp.MaxThreadsPerBlock) 1 1]; % min(128, maxthreadsperblock)

  vxg = gpuArray(zeros(Ns,Nrmax));
  Sxg = gpuArray(zeros(Ns,Nrmax*(Nrmax+1)/2));
  vyg = gpuArray(zeros(Ns,Nrmax));
  Syg = gpuArray(zeros(Ns,Nrmax*(Nrmax+1)/2));
  vzg = gpuArray(zeros(Ns,Nrmax));
  Szg = gpuArray(zeros(Ns,Nrmax*(Nrmax+1)/2));
  
end

if size(vxg,1) ~= Ns | size(vxg,2) ~= Nrmax
  disp('The problem size has changed. Rebuilding GPU arrays');
  vxg = gpuArray(zeros(Ns,Nrmax));
  Sxg = gpuArray(zeros(Ns,Nrmax*(Nrmax+1)/2));
  vyg = gpuArray(zeros(Ns,Nrmax));
  Syg = gpuArray(zeros(Ns,Nrmax*(Nrmax+1)/2));
  vzg = gpuArray(zeros(Ns,Nrmax));
  Szg = gpuArray(zeros(Ns,Nrmax*(Nrmax+1)/2));  
end

if Nd < 3
  xx = [xx zeros(Ns,3-Nd)];
end

m2v = -ones(Nr);
ind = 0;
for jj = 1:Nr
  for ii = jj:Nr
    m2v(ii,jj) = ind;
    ind = ind + 1;
  end
end
m2v = max(m2v,m2v');

[vxg,Sxg,vyg,Syg,vzg,Szg] = feval(gpukGradHess,vxg,Sxg,vyg,Syg,vzg,Szg,real(A),imag(A),real(d),imag(d),...
                                  xx(:,1),xx(:,2),xx(:,3),m2v(:),Ns,Nd,Nr);

vx = gather(sum(vxg,1))';
Sx = gather(sum(Sxg,1))';
vy = gather(sum(vyg,1))';
Sy = gather(sum(Syg,1))';
vz = gather(sum(vzg,1))';
Sz = gather(sum(Szg,1))';

vx = vx(1:Nr);
vy = vy(1:Nr);
vz = vz(1:Nr);
Sx = Sx(1:Nr*(Nr+1)/2);
Sy = Sy(1:Nr*(Nr+1)/2);
Sz = Sz(1:Nr*(Nr+1)/2);
tmp = Sx;
Sx = reshape(tmp(m2v+1),[Nr Nr]);
tmp = Sy;
Sy = reshape(tmp(m2v+1),[Nr Nr]);
tmp = Sz;
Sz = reshape(tmp(m2v+1),[Nr Nr]);



