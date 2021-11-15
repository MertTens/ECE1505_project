% ========================================================================================
%
% NAME: 	generate_Cp.m
% PURPOSE: 	main program to generate Cp (or AIF)
%
% AUTHOR:   Hai-Ling Margaret Cheng
% DATE:     February 26, 2007
% MODIFIED: June 27, 2007
%
% INPUT:	AIF: {'biexponential', 'Yang', 'Leach', 'Cheng', 'Parker', 'Buckley'}
%           Dose: Gd dose (mmol/kg)
%           peak_scale: scaling factor for first-pass peak
%           t: time vector for Cp(t) (minutes).
% OUTPUT:	Cp: plasma concentration curve (mmol/liter)
%
% NOTE:     time resolution must be 5 seconds or less in order to calculate
%           an accurate Cp
%
% ========================================================================================

function [Cp] = generate_Cp (AIF, Dose, peak_scale, t)

del_t = t(2)-t(1);
%LOCAL_AIF = 'y';

switch AIF
    case 'biexponential'
        % ===== from Weinmann's data (Tofts PS. MRM 1991; 17:357) =====
        a = [3.99 4.78];        % kg/liter
        m = [0.144 0.0111];     % min^(-1)
        Cp = Dose*(a(1)*exp(-m(1)*t) + a(2)*exp(-m(2)*t));   % (mmol/liter)
        
    case 'Yang'
        % ===== define AIF commonly seen in aorta (Yang C. MRM 2004; 52:1110) =====
        Cp_aorta = 3.99*exp(-t*0.144) + 4.78*exp(-t*0.0111) + peak_scale*240*(exp(-t/0.125)-exp(-t/0.1)) ...
            + peak_scale*4*((1-exp(-5*(t-0.5))).*unitstep(t-0.5) - (1-exp(-5*(t-0.75))).*unitstep(t-0.75));
        Cp_aorta = Dose * Cp_aorta;                             % (mmol/liter)

        if exist('LOCAL_AIF')
        % ===== define transport function using Gamma variate assumption =====
        to = 0;                 % bolus arrival time (min)
        a = 4;                  % order of Gamma variate function
        b = 0.03;               % dispersion of first-pass (min)
        h = (1/(gamma(a)*(b^a))) * ((t-to).^(a-1) .* exp(-(t-to)./b));
        h(find((t-to)<0)) = 0;
        % ===== convolution of Cp_aorta and h yields Cp in tissue (or AIF) =====
        Cp = conv(Cp_aorta, h) * del_t;
        Cp = Cp(1:length(Cp_aorta));                            % (mmol/liter)

        else
        Cp = Cp_aorta;
        end 
        
    case 'Leach'
        % ===== Martin Leach's 2007 ISMRM abstract (#3498) =====
        ab=2.7e5;     
        ub=40;      % 1/min (Calamante F, et al. MRM 2000; 44:466-73)
        ae=1.5;     % 
        ue=0.144;    % equilibration with whole-body EES (1/min) -- Weinmann
        ar=15;
        ur=8;       % dispersion of recirculation peak relative to bolus peak (1/min)
        tr=0.23;    % recirculation peak delay (min) -- Calamante MRM 2000
        tau=1.8/60;   % dispersion min

        % define initial bolus function 
        cb = ab*ub*ub.*t.*exp(-ub.*t);
        
        % define body transfer functions
        B = ae.*exp(-ue.*t) + ar.*unitstep(t-tr).*(t-tr).*exp(-ur.*unitstep(t-tr).*(t-tr));
     
        % define dispersion of AIF to local tissue
        h=(1/tau).*exp(-t./tau);

        cp_e = conv(cb,B)*del_t; cp_e=cp_e(1:length(t));
        Cp = cp_e + cb*peak_scale;
 %       Cp = conv(Cp,h)*del_t; Cp = Cp(1:length(t));
         
    case 'Cheng'
        % ===== combination of Weinmann, Leach, Calamante =====
        ab=1.08e5;     % ab & ub adjusted to match width & area of bolus in Parker GJ. MRM 2006;56:993
        ub=30;         % 1/min -- equal to (time-to-peak/3)^(-1)
        ae=0.399;      % amplitude for 0.1 mmol/kg
        ue=0.144;    % equilibration with whole-body EES (1/min) -- Weinmann
        ak=0.478;     % amplitude for 0.1 mmol/kg
        uk=0.0111;    % equilibration with kidney elimination (1/min) -- Weinmann
        ar=15;
        ur=8;       % dispersion of recirculation peak relative to bolus peak (1/min)
        tr=0.23;    % recirculation peak delay (min) -- Calamante MRM 2000

        % define initial bolus function (using defn from Calamante MRM 2000).
        cb = ab.*(t.^3).*exp(-ub.*t);
        
        % define body transfer functions for biexponential equilibration
        % with EES and body elimination (Weinmann)
        B1 = ae.*exp(-ue.*t) + ak.*exp(-uk.*t);
           
        % define body transfer function for recirculation (Leach)
        B2 = ar.*unitstep(t-tr).*(t-tr).*exp(-ur.*unitstep(t-tr).*(t-tr));
        
        cp_e = conv(cb,B1)*del_t; cp_e=cp_e(1:length(t));
        cp_r = conv(cb,B2)*del_t; cp_r=cp_r(1:length(t));
        Cp = cb*peak_scale + cp_r + cp_e;
        
    case 'Parker'
        % ===== from Geoff Parker's paper (MRM 2006; 56:993-1000)  =====
        A1=0.809;     % mmol*min
        sigma1=0.0563;  %min 
        T1=0.17046;   % min
        A2=0.330;     % mmol*min
        sigma2=0.132;   % min
        T2=0.365;     % min
        alpha=1.05;     % mmol 
        beta=0.1685;    % 1/min
        s=38.078;       % 1/min
        tau=0.483;      % min
        
        Cp = (A1/sigma1/sqrt(2*pi)).*exp(-(t-T1).^2/(2*sigma1^2));
        Cp = Cp + (A2/sigma2/sqrt(2*pi)).*exp(-(t-T2).^2/(2*sigma2^2));
        Cp = Cp + alpha.*exp(-beta.*t)./(1+exp(-s.*(t-tau)));
        
    case 'Buckley'
        % ===== from Buckley's simulated data =====
        load Glioma-vasc.txt
        t = Glioma_vasc(:,1).';    % time in minutes (time resolution=1 sec)
        Ca = Glioma_vasc(:,2).';   % arterial input in large vessel
        Cp = Glioma_vasc(:,3).';   % plasma curve (dispersion and delay)
        
end

