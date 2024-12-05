
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code interpolates bin values that were set to %
% NaN, then takes the FFT of the data, squares the   %
% Fourier coefficients, takes the average of         %
% frequencies over all profiles in a burst, converts %
% to a one-sided spectral density, and changes units %
% from cycles to radians.                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
dirname='GliderTest';
fname='AquadoppQC1';
load([dirname '/' fname]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up freq scale
Fs = 1/out.cell_size; %sampling rate
[y0,k0] = centeredFFT(0.*out.range,Fs);
N=length(y0);
if mod(N,2)==0
    k1 = k0(N/2:length(y0));
else
    k1 = k0((N+1)/2:length(y0));
end
% k0 is wavenumber for two-sided spectrum
% k1 is wavenumber for one-sided spectrum

% k is one-sided wavenumber in radians/m
k = 2*pi*k1;

Block=18; %Averaging Number for Spectral Estimate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[l,m,n]=size(DATAburst);
FinalSpec_block_cpm = [];
FinalSpec_block = [];
yday_block = [];
DataSetArray=[];    
numprofiles=[];

for burst = 1:l;
    numprofilesi=0; % counter for number of profiles used in spectra for each burst
    %DATAburst was created in Aquadopp_QC
    DATAburstvect=squeeze(DATAburst(burst,:,:));
    TIMEburstvect=squeeze(TIMEburst(burst,:));
    P=nan(size(DATAburstvect));
    
    
    for profile = 1:m;
        if isnan(nanmean(DATAburstvect(profile,:)));
            %P is used for squaring coefficients
            P(profile,:) = nan(size(DATAburstvect(profile,:)));
            %DataSetArray stores the FFT data
            DataSetArray(burst,profile,:) = nan(size(DATAburstvect(profile,:)));
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Prep DATAburstvect for FFT%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % interpolate values for individual bins that were set to NaN
            % (not whole profiles)
            y = DATAburstvect(profile,:);
            yi = fill_data(y);
            
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Detrend data, apply window, remove mean again %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %DataSet = detrend(yi,'linear');
            %DataSetArray stores the prepared data that is used in
            %the FFT
            %DataSet = detrend(yi,'constant');
            %%%%%%%%%%%%%%%%%
            % Line below is to examine phase spectra to investigate "roll up"
            %DataSet=tan(DataSet.*2.*pi.*2.*0.0011./(1500./(2.*10^6)));
            
            DataSet = detrend(yi);
            %Try averaging vel bins (running avg, w/endpoints zeroed)
            %DataSet=ave(DataSet,7,3);DataSet(1)=0;DataSet(end)=0;
            DataSetArray(burst,profile,:) = DataSet;
            
            
            %%%%%%%%%%%%%%%
            % FFT of Data %
            %%%%%%%%%%%%%%%
            
            [YfreqDomain,frequencyRangei] = centeredFFT(DataSet,Fs);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Square FFT Coefficients %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Coeff for 2-sided PSD
            A = out.cell_size/length(DataSet);
            FCoeff = A*real(conj(YfreqDomain).*YfreqDomain);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Create Matrix for Averaging %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Store the squared coefficients
            P(profile,:) = FCoeff;
            
            numprofilesi=numprofilesi+1;
        end
    end
    
    numprofiles(burst)=numprofilesi;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create one-sided spectra %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Average Coefficients over Blocks
    PP=nan(m/Block,n);
    TIME=nan(m/Block,1);
    
    for pp=1:m/Block;
        PP(pp,:)=nanmean(P(Block*(pp-1)+1:Block*pp,:));
        TIME(pp)=nanmean(TIMEburstvect((Block*(pp-1)+1:Block*pp)));
        TIME(pp)=nanmean(TIMEburstvect((Block*(pp-1)+1:Block*pp)));
    end
    
    %1-sided spectra includes the zero frequency
    %Factor of two is because variance is preserved even though we consider
    %only the positive wavenumber band.
    if mod(N,2)==0
        FreqSpec_1sided = 2*PP(:,N/2:N);
    else
        FreqSpec_1sided = 2*PP(:,(N+1)/2:N);
    end

    FinalSpec_block_cpm = [FinalSpec_block_cpm FreqSpec_1sided'];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Convert to Radians (see Veron and Melville paper) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    FinalSpec_block = [FinalSpec_block FreqSpec_1sided'./(2*pi)];
    
    yday_block = [yday_block TIME'];
    pressure_block = [pressure_block P'];
end

spectra.yday=yday_block;
spectra.pressure=pressure_block;
spectra.block=Block;
spectra.k=k1;
spectra.krad=k;
spectra.PSD=FinalSpec_block_cpm;
spectra.PSDrad=FinalSpec_block;

save([dirname '/AquadoppSpectra'],'spectra');

% %%%%%%%%%%%%%%%%%%
% % Variance Check %
% %%%%%%%%%%%%%%%%%%
% %by Parseval's theorem integral should equal variance
% %int gives variance of the spectral density data
% int = trapz(spectra.k,spectra.PSD);
% 
% %var_in gives the variance of the data used in the FFT
% %square velocities, take average over all profiles, then average together
% var_in = nanmean(nanmean(squeeze(DataSetArray(1000,:,:).^2)));




