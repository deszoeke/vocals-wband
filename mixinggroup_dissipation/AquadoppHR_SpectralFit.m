%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Least Squares Fit to Spectra %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all;
dirname='HealyTest/Setup3';
load([dirname '/AquadoppQC1']);
load([dirname '/AquadoppSpectra']);

plotflagspec=1;

[m,n]=size(spectra.PSD);

eps1=nan(1,n);
eps2=nan(1,n);
eps3=nan(1,n);


for burst = 1:n
    
    k = spectra.krad;
    PSD = squeeze(spectra.PSDrad(:,burst))';
    
    if isnan(PSD)==1
        
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % One parameter fit (y-inter only) %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        start = 4;
        stop = length(k); %~15 rad/m to ~97 rad/m
        %         start = 2;
        %         stop = 12; %~15 rad/m to ~97 rad/m
        %       start = 6;
        %        stop = 30; %~15 rad/m to ~97 rad/m
        
        
        x = k(start:stop);
        y = PSD(start:stop);
        
        FIT1 = polyfit(x,y.*x.^(5/3),0);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Two parameter fit                %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Form=[ones(size(x));x.^(-5/3)];
        
        
        [FIT2] = lsqlin(Form',y,[],[],[],[],[0 0],[  ],[0 10^-9]);
        %FIT2 = lsqnonneg(Form',y);

        Noise = FIT2(1);
        FIT2 = FIT2(2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % One parameter fit in Linear Space %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        FIT3 = 10.^nanmean(log10(PSD(start:stop))+(5/3)*log10(k(start:stop)));
        
        %Epsilon=(9*alpha/8)*{(55/18)*PSD*k^5/3}^3/2, alpha=0.4 (Heisenberg's Constant)
        %Eq. 16 in Veron and Melville, 1999, JAOT
        
        eps1(burst) = (9*0.4/8)*((55/18)*FIT1)^(3/2);
        eps2(burst) = (9*0.4/8)*((55/18)*FIT2)^(3/2);
        eps3(burst) = (9*0.4/8)*((55/18)*FIT3)^(3/2);
        
       
        if FIT2<0
            warning('Epsilon cannot be calculated from slope');
            FIT2=nan;
            eps2(burst)=nan;
        end
        
        
        
        if plotflagspec==1
            % Spectral Correction Edson and Fairall, 1998, JAOS (Eq. 25 and Eq. 20) 25
            % Cu_sq=4.*0.55.*eps(burst).^(2./3);
            % L=pi.*0.026;
            % phic=0.25.*Cu_sq.*k.^(-5./3).*(1-0.083.*(k.*L).^2);
            
            eps=nanmean([eps1(burst);eps2(burst);eps3(burst)]);
            
            figure(burst);
            loglog(k,FIT1.*(k).^(-5/3),'b','linewidth',2); hold on;
            loglog(k,FIT2.*(k).^(-5/3)+Noise,'g','linewidth',2);
            loglog(k,FIT3.*(k).^(-5/3),'r','linewidth',2);
            loglog(k(start:stop),PSD(start:stop),'ko','markerfacecolor','k')
            title(['Burst number ' num2str(burst) ', ' datestr(spectra.yday(burst),0) '; \epsilon=' num2str(eps)])
            
            text(.7,.8,['noise = ' num2str(Noise)], 'units','normalized');
             text(.7,.7,['\epsilon = ' num2str(FIT2)], 'units','normalized');
            
            if isnan(FIT2)
                text(.1,.1,'Epsilon cannot be calculated from slope','units','normalized');
            end
            
            [phi11,phi22] = nasmyth(k./2./pi,eps);
            loglog(k,phi11/2/pi,'k--','linewidth',1);
            
            xlabel('Wavenumber (rad/m)')
            ylabel('Spectral density (m^3/s^2)')
            pause;
        end
    end
end






