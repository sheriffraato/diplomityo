function sootload = regen_model(dataStruct, dpfDiam, dpfLen, sootLoadInit, timeStep, p)
                                %     (dataStruct, dpfDiam, dpfLen, sootLoadInit, timeStep, p)
% Inputs:
%    dataStruct: Input dataset in the form:
%        dataStruct.ExhMassFlow: Exhaust massflow [kg/h]
%        dataStruct.ConcC3H6Us: C3H6 concentration DPF upstream [ppm]
%        dataStruct.ConcNoUs: NO concentration DPF upstream [ppm]
%        dataStruct.ConcNo2Us: NO2 concentration DPF upstream [ppm]
%        dataStruct.ConcOxyUs: O2 concentration DPF upstream [ppm]
%        dataStruct.ConcCo2Us: CO2 concentration DPF upstream [ppm]
%        dataStruct.ConcCoUs: CO concentration DPF upstream [ppm]
%        dataStruct.DpfTemp: DPF substrate temperature [°C]
%        dataStruct.ConcSootUs: Engine out soot concentration [ppm]
%        dataStruct.Ashload: Ashload model [g/l]
%        dataStruct.Time: Time vector [s]
%   DpfDiam: diameter of the Dpf [in]
%   DpfLen: length of the Dpf [in]
%   sootLoadInit: Initial soot loading [g/l] (Default: 0)
%   timeStep: Simulation time step [s] (Default: 0.1)
%   paramStruct: DPF kinetic parameters

if nargin < 6
    p = DpfWM_params_HAC();
end

if nargin < 5
    timeStep = 1;
end

if nargin < 4
    sootLoadInit = 0;
end

time = dataStruct.Time;

dpfDiam = dpfDiam * 0.0254;                 % in -> m
dpfLen = dpfLen * 0.0254;                   % in -> m
V_total = dpfLen * pi*(dpfDiam/2)^2;        % Total volume of DPF [m3];
% Number of DPF slices
N = 1;
V_slice = V_total / N;

mf = dataStruct.ExhMassFlow/3600;           % Massflow [kg/h] -> [kg/s]
temp = dataStruct.DpfTemp + 273.15;         % degC -> K
ConcSootUs = 1e-6* dataStruct.ConcSootUs;   % [kg/kg] 
ConcOxyUs = 1e-6* dataStruct.ConcOxyUs;
ConcNo2Us = 1e-6* dataStruct.ConcNo2Us;

msoot = sootLoadInit;

sootload = ones(size(time))*msoot;
%%
k1 = p.DpfWM_GainSootCakeOxyAcvnERgn_P;
k2 = p.DpfWM_GainSootCakeNo2AcvnERgn_P;
Ea1 = p.DpfWM_ExpnSootCakeOxyAcvnERgn_P;
Ea2 = p.DpfWM_ExpnSootCakeNo2AcvnERgn_P;


NM_Coeff = NM_CoeffClcnAddl1(p);
% R = p.DpfWM_Rgas_P;
rho_soot = p.DpfWM_RhoSoot_P;
M_soot = 0.012;

for k = 1:length(time)
    
    [~,CExhGaz] = DpfWall_constants(temp(k),...
            p,NM_Coeff,mf(k),V_slice);

    Ractive = k1 * exp(- Ea1 ./ ( temp(k)) ) * CExhGaz * temp(k);
    Rpassive = k2 * exp(- Ea2 ./ ( temp(k)) ) * CExhGaz;
    
    
    dsoot = mf(k)*ConcSootUs(k) / V_slice -  msoot/rho_soot * M_soot * (Ractive  * ConcOxyUs(k) + ...
                                                                        Rpassive * ConcNo2Us(k));
    msoot = msoot + timeStep * dsoot;
    sootload(k) = msoot;
end


end

function NM_Coeff = NM_CoeffClcnAddl1(p)
%NM_COEFF returns a strcut that contains coefficients of the DpfWM-model.

%Inputs:
%p: parameters of the DpfWM-model

%Outputs:
%NM: coefficients of the DpfWM-model
    
    NM_Coeff.CPSM = (max(p.DpfWM_RhoCell_P,0.00000001))/0.0254^2;
    NM_Coeff.s = sqrt(4/NM_Coeff.CPSM);
    NM_Coeff.d2 = (-2*p.DpfWM_ThcknsWall_P + NM_Coeff.s)/(p.DpfWM_cdr_P + 1);
    NM_Coeff.d1 = p.DpfWM_cdr_P*NM_Coeff.d2;
    NM_Coeff.OFAic = (2*NM_Coeff.d1^2)/NM_Coeff.s^2;
    NM_Coeff.GSAic = (8*NM_Coeff.d1)/NM_Coeff.s^2;
    NM_Coeff.LengthDFL = (p.DpfWM_MaxSootLoadingDFL_P/max(p.DpfWM_RhoSoot_P,0.00000001))/max(NM_Coeff.GSAic,0.00000001);
end


function [Velocity,CExhGaz] = DpfWall_constants(TCell_K,p,NM,MfExh,VolCell)
%DPFWALL_CONSTANTS calculates constants that are needed to calculate soot
%loading.

%Inputs:
%TCell_K: temperature in the Dpf cell [K]
%p: parameters of the DpfWall-model
%NM: coefficients of the DpfWall-model
%MfExh: mass flow per cell [kg/s]

%Outputs:
%Velocity: exhaust gas velocity [m/s]
%CExhGaz: molarity of exhaust gas [mol/m^3]
    
    %The incoming pressure is assumed to be the standard atmosphere.
    PAmbIn = 1013 + 0; %hPa
    PAmbOut = 100*PAmbIn; %Pa
    Pressure = min(max(PAmbOut,1),1e7);
    
    FlowMolExh = (1e3*(MfExh))/p.DpfWM_MmolOfAir_P;
    Velocity = ((FlowMolExh*TCell_K*p.DpfWM_Rgas_P)/Pressure)/(VolCell*NM.GSAic);
    CExhGaz = Pressure/(TCell_K*p.DpfWM_Rgas_P);

end
