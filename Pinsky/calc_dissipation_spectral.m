% DEPRECATED: use Pinksy/calc_dissipation.m
% Calculate cloud top TKE dissipation from Pinsky-adjusted Doppler velocity spectrum
% Simon de Szoeke 2013.04.26

%% set environment and constants
%cd /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband
warning('off','MATLAB:interp1:NaNinY')
run('../../read_parameters')
% way_Pinsky_retrieval='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky/retrieval/';
way_Pinsky_retrieval='/Volumes/Pegasus2/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky/retrieval/';

yyyy=sprintf('%04d',year);
% stday=318; % why missing time adjusted motion till 318? % 281;
% enday=336;
stday=318;
enday=335;
fs=3.5; %Hz
crosscomp=3/4;
kolmogorov=0.54;
factr=crosscomp/kolmogorov;

% universal constant for 2nd order structure function
% longitudinal structure function constant C2 from Pope 6.2 (p.193) after Saddoughi and Veeravalli (1994)
C2ll=2.0;
factrz=1/C2ll;
factrx=3/4/C2ll;

% Kolmogorov (viscous) scale for atmospheric boundary layer is
% eta=(nu^3/epsilon); epsilon~3e-4 m^2/s^2/s, kinematic viscosity nu=2e-5 m^2/s
% --> eta= 2.3 mm
%     25 m = 1000 eta;  1000 m = 440 000 eta

%% w spectra window and FFT parameters used for dissipation
%nwin=64;
%keep=2:8; % indices of spectra to keep - not good if detrended
%cut=10; % index at which to start regarding as noise
nwin=128;
lok=4;
hik=16;
% keep=4:hik; % indices of spectra to keep, start at 4 for detrended & windowed spectra
% spectrum seems to roll off after 16, 32 is in the noise floor.
% Use Hildebrand and Sekhon white noise discrimination below
% to prevent noise contamination.
cut=20; % index at which to start regarding as noise
detredit=1; % 1 detrends, 0 just subtracts mean before calculating spectra

F=(1:nwin/2)'/nwin*fs; % frequencies, Hz; 64: first is nonzero, last is Nyquist
dF=F(2)-F(1); % scalar
F53=F.^(5/3);
dt=10*60; % seconds per window (10 min})
di=floor(fs*dt); % samples per window
S=NaN(size(F));

%% structure function parameters
dxapprox=100; % m [1e3 1e2 1e1]

%% load 10-min ship-relative wind at cloud top interpolated from soundings
%stem='~/Data/cruises/VOCALS_2008/RHB'
stem='/Volumes/Pegasus2/cruises/VOCALS_2008/RHB';
load([stem '/Scientific_analysis/programs/VOCALS2008_programs/wband/uv_cloudtop10.mat']);
load([stem '/Scientific_analysis/programs/VOCALS2008_programs/wband/uv10radar.mat']);
% ship maneuvers removed
% horizontal ship-relative speed
U=sqrt(utoprel.^2+vtoprel.^2); % m/s interpolated relative U from soundings and ship velocity
Uradar=sqrt(uradarrel.^2+vradarrel.^2);

retrfile=dir([way_Pinsky_retrieval 'Pinsky_wRetrieval2008_318_23.mat']); % any file
load([way_Pinsky_retrieval retrfile.name]) % need for h

% DELETE if possible...
% %% After best recalibration, Ken Moran says to subtract 1.72 dB from all recorded values.
% % receiver gain is 1.72 higher than previously estimated, which reduces dBZ
% % and dBmW signals. This is an increase in sensitivity. Note recalibration noise
% % source error of +-0.6 dB.
% dB_offset=-1.72; % dB
% % define a noise level for the reflectivity based on Ken Moran's manual
% radar_const=19.613; % dB_
% analog_noise_signal=-120.74; % dBmW noise peak diagnosed by Simon
% digital_noise_signal=-115.4;  % dBmW
% noise_margin=+3.5; % dB
% adhocthreshold=-43+dB_offset;
% dBZnoiseavg=20*log10(h)+radar_const+analog_noise_signal;
% dBZnoisefloor=20*log10(h)+radar_const+digital_noise_signal+noise_margin;
% noisefloor=max(adhocthreshold,dBZnoisefloor);

%% loop length numbers
nhr=24*(enday-stday+1);
n10min=nhr*6;
%% allocate output data variables
kmax=100;
ndk=5;
count=zeros(n10min,1);
countk=zeros(n10min,kmax);
epsilon=NaN+zeros(n10min,1);
epsilonk=NaN+zeros(n10min,kmax);
epsx=NaN(n10min,kmax);
epsz=NaN(n10min,ndk);
ydayw10=NaN+zeros(n10min,1);
sigv2=NaN+zeros(n10min,1);
sigv2k=NaN+zeros(n10min,kmax);
spec=NaN+zeros(n10min,floor(nwin/2));
speck=NaN+zeros(n10min,floor(nwin/2),kmax);
knoisetop=zeros(n10min,1,'uint8');
Snoisetop=zeros(n10min,1);
Sthreshtop=zeros(n10min,1);
knoisek=zeros(n10min,kmax,'uint8');
Snoisek=zeros(n10min,kmax);
Sthreshk=zeros(n10min,kmax);

%% loop over files by day, hour
for iday=stday:enday
    fprintf(1,'\n%3d ',iday)
    for idhr=0:23
        fprintf(1,' %02d ',idhr)
        retrfile=dir([way_Pinsky_retrieval sprintf('Pinsky_wRetrieval%04i_%03i_%02i.mat',2008,iday,idhr)]);
        if ~isempty(retrfile)
            load([way_Pinsky_retrieval retrfile.name]);
            nt=size(W,1);
            % calculate cloud-top dissipation from inertial range of Pinsky-adjusted vel spectra                        
            % loop through, interp and calc spectra for 10-min chunks
            for ind=1:round(nt/di) % index of the 10 minute interval in the hour
                bigind=24*6*(iday-stday)+6*idhr+ind;
                i10=(ind-1)*di+1; % first index in the 10 minute interval
                ydayw10(bigind,1)=iday+(idhr+(ind-1)/6)/24; % current time of this 10-minute window
                Utop=double(interp1(yday10,U,ydayw10(bigind)+5/(24*60))); % interp ship-rel wind speed to the middle of the current time window
                % find the max coverage of nonnoise returns --> nmax
                ii10=i10:min(i10+di-1,nt);
                % if length(ii10)<di/4; continue; end % break from loop if too few obs.
                if length(ii10)>=di/4 % skip to next iteration if too few obs.
                    nw=sum(isfinite(W(ii10,:)));
                    [nmax,knmax]=max(nw);
                    % choose the highest vertical level that has >0.5*nmax of this coverage --> ik
                    ik=find(nw>0.5*nmax,1,'last');  % highest
                    hk=find(nw>0.5*nmax,1,'first'); % lowest
                    
                    %% compute dissipation at cloud top
                    if ~isempty(ik) && ik<100
                        w10=W(ii10,ik-1:ik+1);
                        iwt=isfinite(w10);
                        % take highest finite w within a level of this 50%ile cloud top
                        wtop=w10(:,1); % w is positive-down
                        wtop(iwt(:,2))=w10(iwt(:,2),2); % overwrites
                        wtop(iwt(:,3))=w10(iwt(:,3),3); % overwrites if there's a higher one
                        % filter out extreme values
                        wtop(abs(wtop-nanmean(wtop))>3*nanstd(wtop))=NaN;
                        % interpolate missing values
                        isn=isnan(wtop);
                        if sum(~isn)>di/4
                            wtop(isn)=interp1(time_offset(ii10(~isn)),wtop(~isn),time_offset(isn));
                            count(bigind,1)=sum(isfinite(wtop));
                            %[S,F]=pwelch(wtop(isfinite(wtop)),nwin,floor(nwin/2),nwin,fs); % uses hamming window 0.54-0.46cos % 0.036 s elapsed
                            %[Sm,Fm]=pwelch(wtop(isfinite(wtop)),hann(nwin),floor(nwin/2),nwin,fs); % hanning cos window (makes no difference)
                            S(:)=gappy_psd(wtop,nwin,fs,70,detrendit); % 2x faster than pwelch, Hanning window and detrends before FFT, handles gappy data
                            % Ignore Nyquist-frequency power estimate!
                            spec(bigind,:)=S';
                            % Snoise=min(S(cut:end-1));
                            % Snoise=median(S(1:end-1));
                            % Sthresh=Snoise+Snoise-min(S(1:end-1)); % symmetric high noise threshold
                            [Sthresh, Snoise, knoise]=HildebrandSekhon(S);
                            hfk=min(hik,knoise); % high-frequency cutoff index (inclusive)
                            keep=lok:hfk;
                            vls=factr.*(2*pi/Utop)^(2/3).*mean(F53(keep).*(S(keep)-Snoise)); % dissipation ^2/3
                            knoisetop(bigind)=knoise;
                            Sthreshtop(bigind)=Sthresh;
                            Snoisetop(bigind)=Snoise;
                            epsilon(bigind,1)=vls.^1.5; % dissipation
                            sigv2(bigind,1)=dF*sum(S(1:knoise)-Snoise); % calculate variance from the whole resolved nonnoise spectrum
                            %                         % Emily Shroyer's method of fitting
                            %                         % f^-5/3 and a noise floor
                            % does not seem to work here
                            %                         Form=[ones(size(F(4:end-1))),F(4:end-1).^(-5/3)];
                            %                         FIT2=lsqlin(Form,S(4:end-1),[],[],[],[],[0 0],[  ],[epsilon(bigind,1) Snoise]);
                            %                         eeps=2*pi/Utop*(factr*FIT2(2))^(3/2);
                            %                         % test plot spectrum and fits
                            %                         loglog(F,S,F,1./F53*4/3*0.54*(Utop/2/pi*epsilon(bigind,1))^(2/3)+Snoise,...
                            %                             F,1./F53*4/3*0.54*(Utop/2/pi*eeps)^(2/3)+FIT2(1))
                            %                         hold on
                            %                         loglog(F,FIT2(2)./F53+FIT2(1),'m')
                        end % there are spectra
                    end % there is 10 min data for cloud top
                    
                    %% compute dissipation elsewhere in+below cloud
                    kk=(max(hk-1,4):min(ik+1,100))';
                    nk=length(kk);
                    Ulev=double(interp1(yday10,Uradar(:,min(kk,100)),ydayw10(bigind)+5/(24*60))); % interp ship-rel wind speed to the middle of the current time window
                    w10=W(ii10,kk);
                    iwt=isfinite(w10);
                    % filter out extreme values
                    w10(bsxfun(@gt,abs(bsxfun(@plus,w10,-nanmean(w10))),4*nanstd(w10)))=NaN;
                    for k=1:nk % loop through levels
                        if countk(bigind,kk(k))>di/4
                            % gappy_psd interpolates missing values
                            S(:)=gappy_psd(w10(:,k),nwin,fs,70,detrendit); % Ignore Nyquist-frequency power estimate!
                            speck(bigind,:,kk(k))=S;
                            %[Sthresh, Snoise, knoise]
                            [sthr,meannoise,kthr]=HildebrandSekhon(S);
                            hfk=min(hik,kthr); % high-frequency cutoff index (inclusive)
                            keep=lok:hfk;
                            vls=factr.*(2*pi/Utop)^(2/3).*mean(F53(keep).*(S(keep)-Snoise)); % dissipation ^2/3 !!mistake to use Utop!!
                            knoisek(bigind,kk(k))=kthr;
                            Snoisek(bigind,kk(k))=meannoise;
                            Sthreshk(bigind,kk(k))=sthr;
                            epsilonk(bigind,kk(k))=vls.^1.5; % dissipation
                            sigv2k(bigind,kk(k))=dF*sum(S(1:knoise)-Snoise); % calculate variance from the whole resolved nonnoise spectrum
                        end % there is data for spectra
                    end % k level where there is data
                    
                    %% transverse 2nd order structure function
                    % make velocity dx increments based on mean U
                    % ?want to do for full hour?
                    if any(isfinite(Ulev))
                        dt=dxapprox./Ulev;
                        di=round(dt*fs); % index increment
                        dx=di.*Ulev/fs;
                        wincrx=w10(1+di:end,:)-w10(1:end-di,:);
                        epsx(bigind,kk)=((factrx*nanmean(wincrx.*wincrx)).^(3/2))./dx;
                        epsx(bigind,dt<4/fs & dt>600)=NaN;
                        epsx(bigind,kk(isnan(Ulev)))=NaN;
                    else
                        epsx(bigind,:)=NaN;
                    end
                    
                    %% longitudinal 2nd order structure function
                    % vertical homogeneity could be a problem
                    for dk=1:ndk
                        dz=h(100)-h(100-dk);
                        wincrz=w10(:,1+dk:end)-w10(:,1:end-dk);
                        epsz(bigind,dk)=nanmean((factrz*nanmean(wincrz.*wincrz,2)).^(3/2)/dz);
                    end
                end % skip if too few obs
            end % 10 min
        end % there is hourly data
    end % hr
end % day

% truncate data
ydayw10=ydayw10(1:bigind);
sigv2=sigv2(1:bigind);
epsilon=epsilon(1:bigind);
count=count(1:bigind);
spec=spec(1:bigind,:);

epsilonk=epsilonk(1:bigind,:);
sigv2k=sigv2k(1:bigind,:);
countk=countk(1:bigind,:);
speck=speck(1:bigind,:,:);
%% QC for enough points in the spectrum
isenuf=count>1800;
isenufk=countk>1800;

% Already QCd Pinsky retrieval for missing motion data and motion
% correction failure, (85% good) and restored some motion values using
% the Crossbow accelerometer/gyro.
%  and, a posteriori, variance <0.4 m^2/s^2
isgdtmp=isenuf;
% nanstd(sigv2(isgdtmp))*3 ~ 0.4
isgd=isgdtmp & sigv2<0.4;
sp=spec(isgd,:);
ep=epsilon(isgd);
yd=ydayw10(isgd);

spk=speck(isgd,:);
epk=epsilonk(isgd,:);

%% save data
save cloudtopturb.mat ydayw10 F spec epsilon sigv2 count knoisetop Snoisetop Sthreshtop
% vertically-resolved data
save radarturb.mat ydayw10 F speck epsilonk sigv2k countk knoisek Snoisek Sthreshk

%% recover noise for each spectrum
%snoise=min(spec(:,cut:end-1),[],2);

%%% plot

%% test QC procedure
plot(ydayw10(isgd),sigv2(isgd),'r.')
hold on
plot(ydayw10,sigv2,'b')
plot(ydayw10,epsilon/1e-2,'k')
plot(ydayw10(isgd),epsilon(isgd)/1e-2,'m.')
set(gca,'ylim',[0 1])

%% scatter plot epsilon and power
plot(epsilon(isgd),sigv2(isgd),'.')
hold on
set(gca,'yscale','log','xscale','log')
axis([0 0.01 0 1])
plot([1.7e-3 1.7e-2].^1.5,[2e-2 2e-1],'r') % epsilon=power^1.5 dependence
xlabel('dissipation')
ylabel('<w''w''>')

%% time series
subplot(2,1,1)
x=ydayw10(isfinite(ydayw10)); y=epsilon(isfinite(ydayw10));
isx=isgdtmp(isfinite(ydayw10)) & sigv2(isfinite(ydayw10))<0.4;
% y(~isx)=NaN;
% plot(x-5/24,y)
y(~isx)=-.01e-3;
area(x-5/24,y,'facecolor',[0.5 0.6 0.7],'edgecolor','none')
hold on
plot(ydayw10(isgd)-5/24,epsilon(isgd),'.','markersize',4)
plot(ydayw10-5/24,5e-4*max(0,-cos(2*pi*(ydayw10-5/24))))
x6=x(1:6:end);
y6=all(~isx([1:6:end; 2:6:end; 3:6:end; 4:6:end; 5:6:end; 6:6:end]))';
plot(x6-5/24,.7e-3*y6-.1e-3,'r.')
set(gca,'tickdir','out','ylim',[0 3e-3],'fontsize',14,'xminortick','on')
xlabel('2008 yearday')
ylabel('TKE dissipation m^2 s^{-3}')
print -depsc CTdissipation/dissip_timeseries.eps

%% QC outcome displayed diurnally
plot(mod(ydayw10*24-5,24),epsilon,'rx') % all points
hold on
plot(mod(ydayw10(isgd)*24-5,24),epsilon(isgd),'.') % only good points
% histogram
subplot(2,1,1)
hist(log10(epsilon),-5.5:0.05:-0.5)
title('all')
set(gca,'xlim',[-5.5 -0.5])
subplot(2,1,2)
hist(log10(ep),-5.5:0.05:-0.5)
set(gca,'xlim',[-5.5 -0.5])
title('QCd')
xlabel('log_{10}(\epsilon)')
% dissipation is lognormally distributed
% QCing based on count and motion correction has little effect on its distribution

%% shade spectra and plot 10x<w'w'>
ii=isfinite(sigv2);
pcolor(ydayw10(ii),F,log10(spec(ii,:)')); shading flat
hold on
plot(ydayw10,10*sigv2,'k')

%% sort the spectra by power in low, high freq tails
clf
subplot(2,1,1)
[ssp,iord]=sort(sp(:,2));
pcolor(log10(sp(iord,:))'); shading flat
hold on
plot(1e4*ep(iord),'k')
set(gca,'xlim',[1 length(iord)],'ylim',[0 33])

subplot(2,1,2)
[ssp,iord]=sort(sp(:,end-1));
pcolor(log10(sp(iord,:))'); shading flat
hold on
plot(1e4*ep(iord),'k')
set(gca,'xlim',[1 length(iord)],'ylim',[0 33])

%% divide spectra into modes with an EOF analysis
ii=isfinite(sigv2);
X=log10(spec(ii,2:end-1))-repmat(mean(log10(spec(ii,2:end-1))),[sum(ii) 1]); % de-mean log10 power series
[u,s,v]=svd(X);

subplot(2,2,1)
pcolor(X'); shading flat
title('full matrix')
subplot(2,2,2)
p=zeros(size(s)); p(1,1)=s(1,1);
pcolor((u*p*v')'); shading flat
title('mode 1')
hold on
plot(v(:,1)/max(v(:,1))*1e3,1:31,'color','k','linewidth',1.4)
subplot(2,2,3)
p=zeros(size(s)); p(2,2)=s(2,2);
pcolor((u*p*v')'); shading flat
title('mode 2')
hold on
plot(v(:,2)/max(abs(v(:,2)))*5e2+5e2,1:31,'color','k','linewidth',1.4)
subplot(2,2,4)
p=zeros(size(s)); p(3,3)=s(3,3);
pcolor((u*p*v')'); shading flat
title('mode 3')
hold on
plot(v(:,3)/max(abs(v(:,3)))*5e2+5e2,1:31,'color','k','linewidth',1.4)

%% composite diurnal cycle of dissipation on local hour: mean, median histogram
indx=bin2(0:24,mod(yd*24-5,24)); % local hour index 1:24; 1 corresp to 0-1 local
ledges=-5.5:0.1:-0.5;
edges=0:1e-4:1e-2;
for ih=1:24
    ii=indx==ih;
    epse(ih)=nanmean(ep(ii));
    epsm(ih,:)=quantile(ep(ii),[.25 .50 .75]);
    % histogram
    hlep(ih,:)=histc(log10(ep(ii)),ledges);
    hep(ih,:)=histc(ep(ii),edges);
    % histogram with median diurnal cycle removed
    y=ep(ii)-epsm(ih,2);
    hdlep(ih,:)=histc(log10(y),ledges);
    hdep(ih,:)=histc(y,edges-2e-3);
end

clf
plot(mod(yd*24-5,24),ep*1e3,'.','color',[0.7 0.7 0.7])
hold on
plot(0.5:23.5,epse*1e3,'ko') % local hour
lh=plot(-0.5:24.5,epsm([end 1:end 1],:)*1e3,'k'); % local hour
set(lh(2),'linewidth',2)
set(gca,'fontsize',16,'xlim',[0 24],'xtick',0:6:24,'xticklabel',[0:6:18 0])
ylabel('TKE dissipation (10^{-3} m^2 s^{-3})')
xlabel('local hour')
set(gca,'ylim',[0 1.5])
%print -depsc CTdissipation/dissipation_diel.eps

colormap(1-gray(16))
% log scale diurnal histogram
subplot(2,1,1)
pcolor(0:36,ledges,hlep([1:end 1:end/2+1],:)'); shading flat
colorbar
set(gca,'ylim',[-5 -2],'xtick',0:6:36,'xticklabel',[0:6:18 0:6:12],'tickdir','out')
subplot(2,1,2)
pcolor(0:36,ledges,hdlep([1:end 1:end/2+1],:)'); shading flat
colorbar
set(gca,'ylim',[-5 -2],'xtick',0:6:36,'xticklabel',[0:6:18 0:6:12],'tickdir','out')

% linear scale diurnal histogram
subplot(2,1,1)
pcolor(0:36,edges,log2(hep([1:end 1:end/2+1],:))'); shading flat
colorbar
hold on
plot(-0.5:36.5,epsm([end 1:end 1:end/2+1],2),'g')
set(gca,'ylim',[0 2e-3],'xtick',0:6:36,'xticklabel',[0:6:18 0:6:12],'tickdir','out')
subplot(2,1,2)
pcolor(0:36,edges-2e-3,hdep([1:end 1:end/2+1],:)'); shading flat
colorbar
set(gca,'ylim',[-0.5e-3 2e-3],'xtick',0:6:36,'xticklabel',[0:6:18 0:6:12],'tickdir','out')

%% log2-frequency diurnal histogram
clf
pcolor(0:36,edges*1e3,log2(hep([1:end 1:end/2+1],:))'); shading flat
a=b2rcolormap(20);
caxis([-0.5 5])
cm=colormap([[1 1 1]; a([11 11:19],:)]);
hc=colorbar;
set(hc,'ylim',[-0.1 4.6],'ytick',[0 1:0.5:4.5],'yticklabel',[1 2 3 4 6 8 12 16 23] ,'fontsize',14)
hold on
plot(mod(yd*24-5,24),ep*1e3,'k.','markersize',3)
plot(24+mod(yd*24-5,24),ep*1e3,'k.','markersize',3)
hml=plot(-0.5:36.5,epsm([end 1:end 1:end/2+1],:)*1e3,'color',[0 0 0]);
set(hml(2),'linewidth',2)
plot(0.5:35.5,epse([1:end 1:end/2])*1e3,'ko') % local hour
set(gca,'ylim',[0 2],'xlim',[0 36],'xtick',0:6:36,'xticklabel',[0:6:18 0:6:12],...
    'tickdir','out','fontsize',16)
xlabel('local hour')
ylabel('cloud top dissipation (10^{-3} m^2 s^{-3})')
%print -depsc CTdissipation/dissipation_diel_shade.eps

% de-meaning by shifting the peak of the distribution takes out recognizable
% cycle of the peak but broadens the distribution of the tails of the distribution of epsilon.
% In the raw series it may be bounded below by noise.

%% distribution shifts over the diurnal cycle
% start at 17 local, 24
b2rcolormap(25);
subplot(2,1,1)
area(ledges,hlep([18:end 1:17],:)')
set(gca,'xlim',[-5 -2],'fontsize',16)
xlabel('log_{10}\epsilon')
subplot(2,1,2)
area(edges,hep([18:end 1:17],:)')
set(gca,'xlim',[0 5e-3],'fontsize',16)
xlabel('\epsilon')

% shifed to remove cycle of diurnal median
subplot(2,1,1)
area(ledges,hdlep([18:end 1:17],:)')
set(gca,'xlim',[-6 -2],'fontsize',16)
xlabel('log_{10}\epsilon')
title('diurnal median removed')
subplot(2,1,2)
area(edges,hdep([18:end 1:17],:)')
set(gca,'xlim',[1e-3 5e-3],'fontsize',16)
xlabel('\epsilon')

