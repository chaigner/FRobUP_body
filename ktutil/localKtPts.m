% This is code by Zhipeng Cao and Will Grissom who have given permission 
% for inclusion within this package. This code can be found also in their
% repo https://bitbucket.org/wgrissom/acptx/ Please cite appropriately.

function [rf,k,m,compWts,nSVT,beta] = localKtPts(d,sens,xx,ks,f0,rfss,algp,prbp,kinit)

% prbp contains:
%                 kmaxdistance, Nrungs, dt, beta, betaadjust, filtertype
% algp contains:
%                 computemethod,nthreads,cgtol

kmaxdistance = prbp.kmaxdistance;
Npulse = prbp.Npulse;
dt = prbp.dt;
beta = prbp.beta;
filtertype = prbp.filtertype;

computemethod = algp.computemethod;
nthreads = algp.nthreads;
cgtol = algp.cgtol;

% Interleaved greedy and local kT-points pulse design. 2D or 3D supported.

% input:
%   d [Nx x Nb] - target magnitude patterns
%   sens [Nx x Nc] - transmit sensitivities
%   xx [Nx x Ngc] - spatial design grid
%   ks [] - candidate phase encode locations
%   f0 [Nx x 1] - spatial frequency offset/off-resonance map
%   fb [Nb x 1] - frequency offset of each band
%   rfss [Nt x 1] - slice-selective rf pulse for one (non-DC) spoke or point
%   dt [scalar] - RF sampling period (seconds)
%   beta [scalar] - RF Tikhonov penalty parameter (dynamically adjusted during iterations)
%   filtertype [char] - Switch to set spectral characteristic

% output:
%   rf [Npulse x Nc]  - designed rung weights
%   k  [Npulse x Ngc] - designed phase encodes

Ngc = size(xx,2);
[Nx,Nc] = size(sens); % Nx: # spatial locs, Nc: # tx coils
Nrp = length(rfss); % number of samples in one rung

compWts = eye(Nc); % initial coil compression matrix for no compression

gambar = 4257;             % gamma/2pi in Hz/T
gam = gambar*2*pi;         % gamma in radians/g

tr = 0:dt:(length(rfss)-1)*dt; % time vector for one rung
A = 1i*gam*dt*exp(1i*2*pi*f0*tr);
m1rung(:,1) = A*rfss;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% design initial pulse
disp('Designing small-tip pulses');

% initialize desired patterns (independent)
if ~isfield(prbp,'phsinit')
    phs = zeros(Nx,1);
else
    phs = prbp.phsinit;
end
dwphs = d.*exp(1i*phs);

% precalculate OMP sens+one rung excitation system matrix
Aomp = pinv(bsxfun(@times,sens,abs(m1rung(:))));

% spoke addition loop
% start with center rung only, then progressively add 2 at a time
% on either side
if ~exist('kinit','var')
    k = zeros(1,Ngc); % get dc spoke to start with
else
    k = kinit;
end
Npulseinit = size(k,1);
rf = zeros(Npulseinit*Nc,1);

nSVT = 0;

Arungsep = zeros(Nx,Npulse);
for Npulset = Npulseinit:Npulse
        
    cost = [Inf Inf];    
    mincost = Inf;
    A = zeros(Nx,Nc*Npulset);
    
    % construct design matrix by modifying one rung excitation
    % patterns with appropriate k-space locations, sensitivities, and
    % off-resonance time offset
    for ii = 1:Npulset
        % blip-induced phase shift
        kphs = xx*k(ii,:)';
        % off res-induced phase shift - account for phase accrual to
        % end of pulse
        totphs = exp(1i*2*pi*(f0*((ii-1)*Nrp - Npulset*Nrp)*dt+kphs));
        tmp = m1rung.*totphs;
        for kk = 1:Nc
            % apply sens, stick it in the design matrix
            A(:,(kk-1)*Npulset+ii) = sens(:,kk).*tmp;
        end
    end

    while ((length(cost) == 2) || (cost(end) < cgtol*cost(end-1)))
        
        % design pulses using pseudoinverse
        if Npulset <= prbp.Ncred % if not enough subpulses to do compression
            
            if Npulset > 1
                rf = (A'*A+beta*speye(Npulset*Nc))\(A'*dwphs);
            else
                rf = A\dwphs;
            end     
            
        else
            
            % take a couple CG iterations
            [xS,~] = qpwls_pcg(rf,A,1,dwphs,0,sqrt(beta),1,...
                algp.ncgiters,ones(numel(rf),1));
            rf = xS(:,end);
      
            rfw = reshape(rf,[Npulset Nc]);
            
            % SVT Compression
            if ~isfield(prbp,'coilMapping')
                [u,s,v] = svd(rfw.','econ');
                s((prbp.Ncred+1):end,(prbp.Ncred+1):end) = 0; % Hard thresholding
                rfw = (u*(s*v')).';rf = rfw(:);
                compWts = u(:,1:prbp.Ncred); % compression weights
            else
                rfwt = rfw;
                compWts = [];
                for ii = 1:size(prbp.coilMapping,2)
                    [u,s,v] = svd(rfw(:,prbp.coilMapping(:,ii)).','econ');
                    s((prbp.Ncred+1):end,(prbp.Ncred+1):end) = 0; % Hard thresholding
                    rfwt(:,prbp.coilMapping(:,ii)) = (u*(s*v')).';
                    compWts(:,ii) = u(:,1:prbp.Ncred);
                end
                rf = rfwt(:);
            end
            nSVT = nSVT + 1;
            
        end
        
        % get excitation patterns and update target phase
        m = reshape(A * rf,[Nx 1]);
        phs = angle(m);
        dwphs = d.*exp(1i*phs);
        
        % optimize phase encodes
        if Npulset > 1
            % get A*b, separated for each rung
            for ii = 1:Npulset
                Arungsep(:,ii) = A(:,ii:Npulset:end)*rf(ii:Npulset:end);
            end
            if isreal(dwphs)
                dwphs = complex(dwphs);
            end
            % get locally-optimal gradient blip area changes
            switch(computemethod)
                case 'mex'
                    if size(xx,2) == 2
                        [vx,Sx,vy,Sy,vz,Sz] = kGradHess(Arungsep(:,1:Npulset),dwphs,[xx zeros(size(xx(:,1)))],nthreads);
                    else
                        [vx,Sx,vy,Sy,vz,Sz] = kGradHess(Arungsep(:,1:Npulset),dwphs,xx,nthreads);
                    end
                case 'gpu'
                    [vx,Sx,vy,Sy,vz,Sz] = kGradHess_cuda(Arungsep(:,1:Npulset),dwphs,xx,Npulse);
                otherwise
                    error 'Unrecognized compute method'
            end;
            dk = [];
            dk(:,1) = inv(Sx)*vx;
            dk(:,2) = inv(Sy)*vy;
            if Ngc == 3
                dk(:,3) = inv(Sz)*vz;
            end
            k = k + dk;
            % add phase to A and update m
            for ii = 1:Npulset
                % A(:,(ii-1)*Nc+1:ii*Nc) = bsxfun(@times,A(:,(ii-1)*Nc+1:ii*Nc),exp(1i*2*pi*(xx*dk(ii,:).'))); % Wrong
                A(:,ii:Npulset:end) = bsxfun(@times,A(:,ii:Npulset:end),exp(1i*2*pi*(xx*dk(ii,:).')));
            end
            m = reshape(A*rf,[Nx 1]);
        end
        
        % evaluate error
        cost(end+1) = 1/2*norm(m(:) - dwphs(:))^2;
        if Npulset > 1
            RFreg = 1/2*beta*real(rf'*rf);
        else
            RFreg = 0;
        end
        if prbp.betaadjust & (rem(length(cost)-2,50) == 0) & ...
                ((RFreg > 1.25*cost(end)) | (RFreg < 0.75*cost(end))) ...
               & (cost(end) < 1.25*mincost) & (Npulset > 1) 
            disp('Adjusting RF penalty');
            beta = cost(end)/RFreg*beta;
            mincost = min([cost mincost]);
            cost = [inf inf];
        else
            cost(end) = cost(end) + RFreg;
        end
        
        if rem(length(cost)-2,10) == 0
            disp(sprintf('Iter %d, %d points. Mean flip angle: %0.4f degrees. Flip angle RMSE: %0.4f degrees. Peak RF: %0.2f.',...
                length(cost)-2, Npulset, mean(abs(m(abs(d) > 0)))*180/pi, sqrt(mean(abs(m - dwphs).^2))/pi*180, max(abs(rf))));
        end
        
    end % while error decreasing loop
    
    
    % add a spoke
    if Npulset < Npulse
        odds = rem(Npulset,2);
        
        % get excitation basis matrix for this spoke location (without
        % k-space loc phase)
        m1rungphs = [];
        
        if (odds & strcmp(filtertype,'Lin phase')) | strcmp(filtertype,'Min phase') % add phase for first rung in pulse
            offresphs = exp(1i*2*pi*f0*(-(Npulset+1)*Nrp)*dt);
        else % add phase for last rung in pulse
            offresphs = exp(1i*2*pi*f0*(-Nrp)*dt);
        end;
        m1rungphs(:,1) = m1rung(:,1).*offresphs;
        
        m1rungphs = m1rungphs(:); % stack frequency bands
        
        if (~odds & strcmp(filtertype,'Lin phase')) | strcmp(filtertype,'Max phase') % if adding a rung on the end of the pulse
            % advance m phase to account for later observation time
            offresphs = exp(1i*2*pi*f0*(-Nrp)*dt);
            m(:,1) = m(:,1).*offresphs;
            
            % advance desired phase to account for later observation time
            f0phsadd = -2*pi*f0*Nrp*dt; % additional phase due to main
            % field offsets
            
            % case 'independent'
            phs(:,1) = phs(:,1) + f0phsadd;
            % completely independent phase relaxation across bands
            dwphs(:,1) = d(:,1).*exp(1i*phs(:,1));
        end
        
        disp('Adding a kt-point using orthogonal matching pursuit');
        
        % calculate magnitude of pulses for each k-loc
        % try basing selection on actual error after adding this spoke?
        % can i calculate terms for solo for the rest of the spokes
        % before looking at new candidates, to reduce cost?
        res = dwphs(:) - m(:); % initialize residual to current error
        res = exp(-1i*angle(m1rungphs)).*res;
        rfnorm = [];
        if (odds & strcmp(filtertype,'Lin phase')) | strcmp(filtertype,'Min phase')
            % spoke will be added to beginning of pulse
            kclosest = k(1,:);
        else
            % add spoke to end of pulse
            kclosest = k(end,:);
        end;
        for jj = 1:length(ks)
            if ~any(abs(ks(jj,:) - kclosest) > kmaxdistance) % if this point isnt too far
                rfnorm(jj) = norm(Aomp*(exp(-1i*2*pi*(xx*ks(jj,:).')).*res));
            else
                rfnorm(jj) = 0;
            end
        end
        % add the k-loc with max norm to the list
        [~,maxi] = max(rfnorm);
        
        if (odds & strcmp(filtertype,'Lin phase')) | strcmp(filtertype,'Min phase')
            % add spoke to beginning of pulse
            k = [ks(maxi,:); k];
            rf = [zeros(1,Nc); reshape(rf,[Npulset Nc])];
            rf = rf(:);
        else
            % add spoke to end of pulse
            k = [k; ks(maxi,:)];
            rf = [reshape(rf,[Npulset Nc]); zeros(1,Nc)];
            rf = rf(:);
        end
    end
end



