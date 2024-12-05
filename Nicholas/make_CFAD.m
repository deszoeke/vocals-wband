% CFADS have been compiled by proc_1min_stat_v6.m --> compile_1min_stat_v6.m.
% among the 1-minute statistics in
% .../RHB/radar/wband/Processed/1min_stat/Z_1min.mat
%                                         w_1min.mat
%                                         Dw_1min.mat
%
% June 2011 Simon de Szoeke

read_parameters
load([way_proc_data_wband '/1min_stat/Z_1min.mat'])
load([way_proc_data_wband '/1min_stat/w_1min.mat'])
load([way_proc_data_wband '/1min_stat/Dw_1min.mat'])

it=7000; % choose a time index

clf
b2rcolormap(21);

subplot(1,3,1)
imagesc(Z.bins(2:end-1),Z.height,squeeze(Z.cfad(it,:,2:end-1)));
set(gca,'ydir','normal')
xlabel('reflectivity (dBZ)')

subplot(1,3,2)
imagesc(w.bins(2:end-1),w.height,squeeze(w.cfad(it,:,2:end-1)));
set(gca,'ydir','normal')
xlabel('Doppler velocity (ms^{-1})')
title(datestr(Z.time_yday(it)+datenum('0-jan-2008')))

subplot(1,3,3)
imagesc(Dw.bins(2:end-1),Dw.height,squeeze(Dw.cfad(it,:,2:end-1)));
set(gca,'ydir','normal')
xlabel('Doppler width (ms^{-1})')

orient landscape
print -depsc CFAD.eps