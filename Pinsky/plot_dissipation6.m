% Calculate cloud top TKE dissipation from Pinsky-adjusted Doppler velocity spectrum
% Simon de Szoeke 2013.04.26

%% set environment and constants
%cd /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky
%cd /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/wband/Pinsky
warning('off','MATLAB:interp1:NaNinY')
run('../../read_parameters')
machine=char(java.net.InetAddress.getLocalHost.getHostName);
if regexp(machine,'squall')==1
    stem='~/Data/cruises/VOCALS_2008/RHB';
    way_Pinsky_retrieval='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/wband/Pinsky/retrieval/';
elseif regexp(machine,'fog')==1
    stem='/Volumes/Pegasus2/cruises/VOCALS_2008/RHB';
    way_Pinsky_retrieval='/Volumes/Pegasus2/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky/retrieval/';
end

yyyy=sprintf('%04d',year);
% stday=318; % why missing time adjusted motion till 318? % 281;
% enday=336;
stday=318;
enday=335;
fs=3.5; %Hz
crosscomp=3/4;
kolmogorov=0.54;
factr=crosscomp/kolmogorov;

%% structure function parameters
dxapprox=[1e3 5e2 1e2 5e1 1e1]; % m 
ndx=length(dxapprox);
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
% nwin=128;
% keep=4:16; % indices of spectra to keep, start at 4 for detrended & windowed spectra, ?32 eats into noise floor?
% cut=20; % index at which to start regarding as noise
nwin=256;
keep=8:32; % indices of spectra to keep, start at 4 for detrended & windowed spectra, ?32 eats into noise floor?
cut=40; % index at which to start regarding as noise


F=(1:nwin/2)'/nwin*fs; % frequencies, Hz; 64: first is nonzero, last is Nyquist
dF=F(2)-F(1); % scalar
F53=F.^(5/3);
dt=10*60; % seconds per window (10 min})
di=floor(fs*dt); % samples per window
S=NaN(size(F));

%% load 10-min ship-relative wind at cloud top interpolated from soundings
if regexp(machine,'squall')==1
    stem='~/Data/cruises/VOCALS_2008/RHB';
    load([stem '/Scientific_analysis/programs/wband/uv_cloudtop10.mat']);
    load([stem '/Scientific_analysis/programs/wband/uv10radar.mat']);
elseif regexp(machine,'fog')==1
    stem='/Volumes/Pegasus2/cruises/VOCALS_2008/RHB';
    load([stem '/Scientific_analysis/programs/VOCALS2008_programs/wband/uv_cloudtop10.mat']);
    load([stem '/Scientific_analysis/programs/VOCALS2008_programs/wband/uv10radar.mat']);
end
% ship maneuvers removed
% horizontal ship-relative speed
U=sqrt(utoprel.^2+vtoprel.^2); % m/s interpolated relative U from soundings and ship velocity
Uradar=sqrt(uradarrel.^2+vradarrel.^2);

%% load data
load('cloudtopturb.mat')
% old version cloudtopturb_nfft128_2013-05-07.mat used nfft=128 points

% vertically-resolved data
load('radarturb.mat'); % for backward compatibility of some plots
v1=load('radarturb.mat');  % constant thresholds for noise, too big Ulev
v2=load('radarturb2.mat'); % HS74-like R_2=2 noise threshold, too conservative
v3=load('radarturb3.mat'); % probabilistic noise threshold p1side=0.98
v4=load('radarturb4.mat'); % arbitratily cut some high frequencies out of dissipation calc.
v4p5=load('radarturb4p5.mat'); % nwin=192 (do not cut frequencies)
v5=load('radarturb5.mat'); % nwin=256
v6=load('radarturb6.mat'); % nwin=192, ninterpgap=20 results identical to v4p5 with ninterpgap=70
% longer windows give larger dissipation. This could be a technical issue,
% but Jim Moum couldn't think of any reason for it to be so. Another
% possibility is that dissipation is higher in the middle of clouds, and
% shorter windows sample weaker dissipation closer to the edges of the
% cloud.
DZ=(h(100)-h(99))*(1:5);
% old version radarturb_100m_2013-05-07.mat used only 100m structure function

%% recover noise for each spectrum
snoise=min(spec(:,cut:end-1),[],2);
snoisek=min(speck(:,cut:end-1),[],2);

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
% used for histogram

%%% plot
%% vertical structure of epsilon in cloud
clf
subplot(3,1,1)
hsol=area(ydayw10(isfinite(ydayw10))-304,max(0,-cos(2*pi*(ydayw10(isfinite(ydayw10))-5/24))),0,'edgecolor','w','facecolor','w');
set(gca,'color',0.8*[1 1 1])
hold on
% pcolor(ydayw10-304,h(1:100)/1e3,log10(epsx(:,:,3)'));shading flat
% pcolor(ydayw10-304,h(1:100)/1e3,log10(epsilonk'));shading flat; caxis([-5 -2.5])
pcolor(ydayw10-304,h(1:100)/1e3,log10(v2.epsilonk'));shading flat; caxis([-5 -2.5])
hc=colorbar; set(hc,'fontsize',14)
set(hsol,'edgecolor','w','facecolor','w','visible','on')
set(get(hsol,'children'),'edgecolor','w','facecolor','w')
set(gca,'xminortick','on')
ylim([0 2])
set(gca,'fontsize',14,'tickdir','out','color',0.8*[1 1 1])
ylabel('height (km)')
xlabel('November 2008')
title('v2 spectral dissipation log_{10}(\epsilon/[m^2 s^{-3}])')
orient landscape
box on
set(gca,'xlim',[15 32],'xtick',15:32)
set(gcf,'inverthardcopy','off','color','w')

subplot(3,1,2)
hsol=area(ydayw10(isfinite(ydayw10))-304,max(0,-cos(2*pi*(ydayw10(isfinite(ydayw10))-5/24))),0,'edgecolor','w','facecolor','w');
set(gca,'color',0.8*[1 1 1])
hold on
pcolor(ydayw10-304,h(1:100)/1e3,log10(v3.epsilonk'));shading flat; caxis([-5 -2.5])
hc=colorbar; set(hc,'fontsize',14)
set(hsol,'edgecolor','w','facecolor','w','visible','on')
set(get(hsol,'children'),'edgecolor','w','facecolor','w')
set(gca,'xminortick','on')
ylim([0 2])
set(gca,'fontsize',14,'tickdir','out','color',0.8*[1 1 1])
ylabel('height (km)')
xlabel('November 2008')
title('v3 spectral dissipation log_{10}(\epsilon/[m^2 s^{-3}])')
orient landscape
box on
set(gca,'xlim',[15 32],'xtick',15:32)

subplot(3,1,3)
hsol=area(ydayw10(isfinite(ydayw10))-304,max(0,-cos(2*pi*(ydayw10(isfinite(ydayw10))-5/24))),0,'edgecolor','w','facecolor','w');
set(gca,'color',0.8*[1 1 1])
hold on
pcolor(ydayw10-304,h(1:100)/1e3,(v3.epsilonk-v2.epsilonk)');shading flat; caxis(1e-3*[-1 1])
% pcolor(ydayw10-304,h(1:100)/1e3,(v2.epsilonk-epsilonk)');shading flat; caxis(1e-3*[-1 1])
hc=colorbar; set(hc,'fontsize',14)
set(hsol,'edgecolor','w','facecolor','w','visible','on')
set(get(hsol,'children'),'edgecolor','w','facecolor','w')
set(gca,'xminortick','on')
ylim([0 2])
set(gca,'fontsize',14,'tickdir','out','color',0.8*[1 1 1])
ylabel('height (km)')
xlabel('November 2008')
title('v3-v2 spectral dissipation (\epsilon/[m^2 s^{-3}])')
orient landscape
box on
set(gca,'xlim',[15 32],'xtick',15:32)
set(gcf,'inverthardcopy','off','color','w')

%{
subplot(3,1,3)
plot(yday10,Uradar)
set(gca,'xlim',[320, 337])
colorbar
%}

%% compare v5 epsilon with v3
clf
% loglog(v3.epsilonk(:),v5.epsilonk(:),'.','markersize',0.3)
loglog(v4p5.epsilonk(:),v6.epsilonk(:),'.','markersize',0.3)
axis square
axis([1e-6 1e-1 1e-6 1e-1])
hold on
plot([1e-6 1e-1],[1e-6 1e-1],'k-')
ylabel('v5'); xlabel('v3')
title('dissipation [m^2 s^{-3}]')
saveas(gcf,'compare_v3-v5_eps.png','png')

ledges=-8:.1:2;
% h1=histc(log10(v1.epsilonk(:)),ledges); 
% h2=histc(log10(v2.epsilonk(:)),ledges);
% h3=histc(log10(v3.epsilonk(:)),ledges); % mean=2.4684e-04, skewness(log10(v3.epsilonk(isfinite(v3.epsilonk(:))))) = -0.30
h4p5=histc(log10(v4p5.epsilonk(:)),ledges); % mean=2.8351e-04
% h5=histc(log10(v5.epsilonk(:)),ledges);     % mean=3.1235e-04
h6=histc(log10(v6.epsilonk(:)),ledges);     % mean=3.1235e-04
clf
plot(ledges,[h4p5,h6]); xlim([-6 -2])
set(gca,'fontsize',14)
% legend('v3: 128','v4p5: 192','v5: 256')
xlabel('log10 dissipation')
ylabel('counts')

% compare spectra
% st2=permute(v2.speck,[2,1,3]);
st3=permute(v3.speck,[2,1,3]);
st4p5=permute(v4p5.speck,[2,1,3]);
st5=permute(v5.speck,[2,1,3]);
st6=permute(v6.speck,[2,1,3]);
ii=isfinite(squeeze(st3(1,:)))' & isfinite(squeeze(st5(1,:)))' & v3.knoisek(:)>4 & v5.knoisek(:)>4;
fii=find(ii);
s3=st3(:,fii); s4p5=st4p5(:,fii); s5=st5(:,fii);
knoi3=v3.knoisek(fii)';
knoi4p5=v4p5.knoisek(fii)';
knoi5=v5.knoisek(fii)';
% clf
subplot(3,1,1); set(gca,'fontsize',14)
for i = 1:500:length(s3)
    loglog(v5.F(1:knoi5(i)),s5(1:knoi5(i),i),'b');
    hold on
    loglog(v5.F(knoi5(i):end),s5(knoi5(i):end,i),'r');
end
xlim([1e-2,2])
title('nwin=256')
set(gca,'fontsize',14)
subplot(3,1,2)
for i = 1:500:length(s3)
    loglog(v4p5.F(1:knoi4p5(i)),s4p5(1:knoi4p5(i),i),'b');
    hold on
    loglog(v4p5.F(knoi4p5(i):end),s4p5(knoi4p5(i):end,i),'r');
end
xlim([1e-2,2])
title('nwin=192')
ylabel('S_{ww}(f)')
set(gca,'fontsize',14)
subplot(3,1,3)
for i = 1:500:length(s3)
    loglog(v3.F(1:knoi3(i)),s3(1:knoi3(i),i),'b');
    hold on
    loglog(v3.F(knoi3(i):end),s3(knoi3(i):end,i),'r');
end
xlim([1e-2,2])
title('nwin=128')
set(gca,'fontsize',14)
xlabel('frequency (s^{-1})')
text(7e-1,1.4e-1,'white noise','fontsize',14)
saveas(gcf,'S_f_noise.eps','epsc')

%% vertical top of dissipation meas.
isf=isfinite(epsx(:,:,3));
jmax=zeros(size(isf,1),1);
for j=1:size(isf,2)
    jmax(isf(:,j))=j;
end
hmax(jmax>0)=h(jmax(jmax>0));
hmax(jmax==0)=NaN;

% composite cloud-top relative diurnal epsilon
[epstmp,eps2tmp]=deal(zeros(144,100));
cnt=zeros(144,100);
tbin=(0:10:60*24)/60/24; % diurnal 10-min bins
isf=isfinite(ydayw10) & any(isfinite(epsx(:,:,3)),2);
fisf=find(isf);
iit=floor(interp1(tbin,1:length(tbin),mod(ydayw10(isf),1),'nearest'));
for it=1:length(iit)
    tmp=epsx(fisf(it),max(1,jmax(fisf(it))+(-99:0)),3);
    ii=isfinite(tmp);
    tmp(~ii)=0;
    epstmp(iit(it),:)=epstmp(iit(it),:)+tmp;
    eps2tmp(iit(it),:)=eps2tmp(iit(it),:)+tmp.*tmp;
    cnt(iit(it),:)=cnt(iit(it),:)+ii;
end
epscmp=epstmp./cnt;
epscmp(cnt<1)=NaN;
epsstd=sqrt(eps2tmp./cnt-epscmp.*epscmp);
epsstd(cnt<4)=NaN;

clf
subplot(2,1,1)
hsol=area(tbin*24,-2+0.4*max(0,-cos(2*pi*(tbin-5/24))),-2,'edgecolor','w','facecolor','w');
hold on
pcolor(tbin(1:end)*24,(h(19:100)-h(100))/1e3,log10(epscmp([1:end 1],19:end))'); shading flat
colorbar
axis([0 24 -2 0])
caxis([-5 -2.5])
set(gca,'color',0.8*[1 1 1],'fontsize',14,'xtick',0:4:24,'tickdir','out')
set(hsol,'edgecolor','none')
ylabel('height below cloud top (km)')
title('log10 dissipation [m^2 s^{-3}]')

subplot(2,1,2)
hsol=area(tbin*24,-2+0.4*max(0,-cos(2*pi*(tbin-5/24))),-2,'edgecolor','w','facecolor','w');
hold on
pcolor(tbin(1:end)*24,(h(19:100)-h(100))/1e3,(epsstd([1:end 1],19:end)./epscmp([1:end 1],19:end))'); shading flat
colorbar
axis([0 24 -2 0])
caxis([0 2])
set(gca,'color',0.8*[1 1 1],'fontsize',14,'xtick',0:4:24,'tickdir','out')
set(hsol,'edgecolor','none')
title('<dissipation''^2>^{1/2}/dissipation')
xlabel('hour UTC')

set(gcf,'inverthardcopy','off','color','w')
orient tall
print -depsc dissipation_diurnal_cmp.eps

% Standard deviation of epsilon is on the order of 0.5*
% epsilon itself, indicating to me turbulence is not terribly episodic.

% note inversion properties are saved in
% /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/balloon/Processed/sonde_inversion.mat