% Composite W-band reflectivity Z on column max. reflectivity CMR
read_parameters;
load([way_proc_data_wband '1min_stat/Z_1min.mat']);

% Define dimensions
% vertical coordinate (center of ranges, m)
% h=load([way_proc_data_wband '1min_stat/height.txt']);
% % dBZ and w bin edges
% Zbins=[-inf -40:2:30 inf];
% wbins=[-inf -5:0.2:5 inf];

nheight=length(Z.height); % 75
nmin=60;
nhour=24;
nZbin=length(Z.bins);

% load 1-min W-band cloud top height (Simon de Szoeke, Oregon State U.)
% base_time=1225916508;
A=load([way_proc_data_wband 'cloudheight/CloudHeight_1min_2008310-336.txt']);
A(A==-999)=NaN; % NaN out missing value of -999
% already synchronzied to time_yday from 1min_stat/Z_1min.mat
% time=A(:,1); % seconds since base_time
cloud_fraction=A(:,2);
cloud_top_mean=A(:,3);
cloud_top_std=A(:,4);
cloud_top_val15=A(:,5);
cloud_top_median=A(:,6);
cloud_top_val85=A(:,7);
cloud_top_valmean=A(:,8); %mean of tops within 15-85 percentile
num_returns=A(:,9);
num_cloud_returns=A(:,10);
num_clear_returns=A(:,11);
num_nonmet_returns=A(:,12);
flag=A(:,13);
% flag values
% 0 no errors detected, cloud top computed.
% 1 cloud fraction=0, cloud top set to NaN.
% 2 cloud fraction>=0 but found no 3 cloud returns with consistent height,
%   possible glitch, provisional cloud top computed.
% NaN out percentile-based cloud tops=0 using flag
lf=logical(flag);
cloud_top_val15(lf)=NaN;
cloud_top_val85(lf)=NaN;
cloud_top_valmean(lf)=NaN;
cloud_top_median(lf)=NaN;
% time_yday=datenum(0,0,0,0,0,base_time+time_aggreg)+datenum(1970,1,1,0,0,0)-datenum(2008,1,0,0,0,0); %yearday

% define a noise level for the reflectivity based on Ken Moran's manual
radar_const=19.613; % dB_
dB_offset=-1.72; % dB
min_detectable_signal=-113.304+dB_offset; % dBmW
analog_noise=-120.74;  % peak, dBmW
digital_noise=-115.4;   % peak, dBmW
threshold_margin=3.5; % dB
noisefloor=20*log10(Z.height)+radar_const+digital_noise+threshold_margin;
% do not plot colors/contours below the noise threshold
isnoise=Z.mean<repmat(noisefloor',[size(Z.mean,1) 1]);
Z.plotmean=Z.mean;
Z.plotmean(isnoise)=NaN;
Z.plotvar=Z.var;
Z.plotvar(isnoise)=NaN;

% Column maximum reflectivity (CMR)
% vertical max reflectivity of 1-minute means
[zmx,k]=max(Z.mean,[],2);

% composite 1-min Z profiles on cloud_top_mean
topbin=[400 1000 1200 1400 1600 1800 2000]; % 6 bins defined between these 7 endpoints
ntime=length(time_yday);
nheight=length(Z.height);
for indbin=1:length(topbin)-1
    ii=cloud_top_mean>=topbin(indbin) & cloud_top_mean<topbin(indbin+1);
    countinbin(indbin)=sum(ii);
    Z.topbin(indbin,:)=nanmean(Z.mean(ii,:)); % simple level-mean
    Z.topcfad(indbin,:,:)=nansum(Z.cfad(ii,:,:));
    hadj=repmat(Z.height',[sum(ii) 1])-repmat(cloud_top_valmean(ii),[1 nheight]); % 0 at cloud top
end
% composite based on distance from W-band cloud top

h=area([-60; noisefloor],Z.height([1 1:end]),2000,'facecolor',0.7+[0 0 0],'edgecolor','none');
hold on
plot(Z.topbin,Z.height,'linewidth',2)
axis([-50 20 0 2e3])
set(gca,'fontsize',16)
ylabel('altitude (m)')
xlabel('cloud top height composite mean reflectivity (dBZ)')
%print('-dpng',[way_proc_images_wband 'composite_Z_on_cloudtop.png'])

% CFAD of Z
for indbin=1:6
    ax(indbin)=subplot(2,3,indbin,'align');
    
    [cc,hc]=contourf(Z.bins(2:end-1)+(Z.bins(3)-Z.bins(2))/2,Z.height/1e3,squeeze(Z.topcfad(indbin,:,2:end-1)),40);
    set(hc,'edgecolor','none')
    colormap(b2rcolormap(40));
    hb=colorbar('southoutside');
    set(get(hb,'xlabel'),'string','count')
    title(sprintf('z_{top}=[%3.1f %3.1f) km',topbin(indbin+[0 1])/1e3));
    hold on
    plot(noisefloor,Z.height/1e3,'w-')
    plot(Z.topbin(indbin,:),Z.height/1e3,'m-')
    xlabel('dBZ')
    ylabel('height (km)')
end
%print('-dpng',[way_proc_images_wband 'comp_Zcfad_on_cloudtop.png'])

% plot to compare to LES
clf
[cc,hc]=contourf(Z.bins(2:end-1)+(Z.bins(3)-Z.bins(2))/2,Z.height/1e3,squeeze(Z.topcfad(indbin,:,2:end-1)),40);
set(hc,'edgecolor','none')
colormap(b2rcolormap(40));
hb=colorbar('southoutside');
set(get(hb,'xlabel'),'string','count')
title(sprintf('z_{top}=[%3.1f %3.1f) km',topbin(indbin+[0 1])/1e3));
hold on
% contour(Z.bins(2:end-1)+(Z.bins(3)-Z.bins(2))/2,Z.height/1e3,squeeze(Z.topcfad(indbin,:,2:end-1)),6,'k')
plot(noisefloor,Z.height/1e3,'w-')
plot(Z.topbin(indbin,:),Z.height/1e3,'m-')
set(gca,'fontsize',16)
title({'frequency of occurence' 'cloud top height 1.6-1.8 km'})
xlabel('dBZ')
ylabel('height (km)')
print('-dpng',[way_proc_images_wband 'Z_CFAD_cloutop1600-1800.png'])

    
% movie of CFADs in 1 dBZ CMR bins
edges=500:50:1900;
count=histc(cloud_top_mean,edges);

hf=figure; dock(-hf);

ax(1)=subplot(2,1,1); hold on;
ax(2)=subplot(3,1,3);
hbar=bar(edges,count/60,'histc');
set(hbar,'facecolor',0.7*[1 1 1])
axis([500 2000 0 50])
set(gca,'ytick',0:10:50)
set(gca,'tickdir','out','fontsize',14)
ylabel('frequency (hours)')
xlabel('cloud top height (km)') 
set(gca,'tickdir','out')
hold on

indbin=1;
    axes(ax(2));
    ha=plot(edges(indbin)+25,-3,'k^','markersize',8,'markerfacecolor','k','clipping','off');
    
    axes(ax(1));
    ii=cloud_top_mean>=edges(indbin) & cloud_top_mean<edges(indbin+1);
    cfad(indbin,:,:)=nansum(Z.cfad(ii,:,2:end-1),1);
    if any(any(cfad(indbin,:,:)~=0))
        [cc,hc]=contourf(Z.bins(2:end-1)+(Z.bins(3)-Z.bins(2))/2,Z.height/1e3,squeeze(cfad(indbin,:,:)),40);
        set(hc,'edgecolor','none')
        ht=text(1,1,'');
    else
        hc=image(Z.bins(2:end-1)+(Z.bins(3)-Z.bins(2))/2,Z.height/1e3,squeeze(cfad(indbin,:,:))); % blue screen
        axis xy; axis([-39.0000   31.0000    0.1336    1.9818]);
        ht=text(-4,1,'no scenes','fontsize',16,'color','w','horizontalalignment','center');
    end
    set(gca,'fontsize',16)
    colormap(b2rcolormap(40));
    title(sprintf('z_{top}=[%3d %d) m',edges(indbin+[0 1])));
    xlabel('dBZ')
    ylabel('height (km)')
    hb=colorbar('southoutside','fontsize',14);
    set(hb,'position',get(hb,'position')-[0 .04 0 0])
    set(get(hb,'xlabel'),'string','count','fontsize',16)
    hold on
    plot(noisefloor,Z.height/1e3,'w-')
    % make adjustments
    %print('-djpeg',[way_proc_images_wband sprintf('CFAD_movie/ZCFAD%02dcloudtop%+03d.jpeg',indbin,edges(indbin))])
    frame(indbin)=getframe(hf);
    delete(hc,ha,ht,hb);

for indbin=2:length(edges)-1
    axes(ax(2));
    ha=plot(edges(indbin)+25,-3,'k^','markersize',8,'markerfacecolor','k','clipping','off');
    
    axes(ax(1));
    ii=cloud_top_mean>=edges(indbin) & cloud_top_mean<edges(indbin+1);
    cfad(indbin,:,:)=nansum(Z.cfad(ii,:,2:end-1),1);
    if any(any(cfad(indbin,:,:)~=0))
        [cc,hc]=contourf(Z.bins(2:end-1)+(Z.bins(3)-Z.bins(2))/2,Z.height/1e3,squeeze(cfad(indbin,:,:)),40);
        set(hc,'edgecolor','none')
        ht=text(1,1,'');
    else
        hc=image(Z.bins(2:end-1)+(Z.bins(3)-Z.bins(2))/2,Z.height/1e3,squeeze(cfad(indbin,:,:))); % blue screen
        axis xy; axis([-39.0000   31.0000    0.1336    1.9818]);
        ht=text(-4,1,'no scenes','fontsize',16,'color','w','horizontalalignment','center');
    end
    set(gca,'fontsize',16)
    colormap(b2rcolormap(40));
    title(sprintf('z_{top}=[%d %d]',edges(indbin+[0 1])));
%     xlabel('dBZ')
%     ylabel('height (km)')
    hb=colorbar('southoutside','fontsize',14);
    set(hb,'position',get(hb,'position')-[0 .04 0 0])
    set(get(hb,'xlabel'),'string','count','fontsize',16)
    hold on
    plot(noisefloor,Z.height/1e3,'w-')
    %print('-djpeg',[way_proc_images_wband sprintf('cloudtop_movie/ZCFAD%02dcloudtop%+03d.jpeg',indbin,edges(indbin))])
    frame(indbin)=getframe(hf);
    delete(hc,ha,ht,hb);
end

movie2avi(frame,[way_proc_images_wband 'cloudtop_movie/ZCFAD_on_cloudtop.avi'],'fps',5);