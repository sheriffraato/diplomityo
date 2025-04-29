

function sootload = regeneration_model(dataStruct, dpfDiam, dpfLen, sootLoadInit, timeStep)

if nargin < 5
    timeStep = 0.1;
end

if nargin < 4
    sootLoadInit = 0;
end
time = dataStruct.Time;

dpfDiam = dpfDiam * 0.0254;                 % in -> m
dpfLen = dpfLen * 0.0254;                   % in -> m
V_total = dpfLen * pi*(dpfDiam/2)^2;        % Total volume of DPF [m3];
mf = dataStruct.ExhMassFlow/3600;            % Massflow [kg/h] -> [kg/s]
temp = dataStruct.DpfTemp + 273.15;         % degC -> K
ConcSootUs = 1e-6* dataStruct.ConcSootUs;

R = 8.31446261815324; % Molar gar constant [J * 1/K * 1/mol]
Ea1 = 24654.699; % soot cake active regeneration activation energy 
Ea2 = 10792; % soot cake passive regeneration activation energy 
k1 = 2072679936; % Soot cake active regeneration gain [O2]
k2 = 3012830538; % Soot cake passive regeneration gain [NO2]
soot_M = 12.01070; % molar mass of soot [g/mol]
ro_soot = 34.633; % density of the soot layer [kg/m3]

mSoot = sootLoadInit;

sootload = ones(size(time))*mSoot;
for k = 1:length(time)
    Vsoot = mf(k) * ConcSootUs(k) ./ ro_soot;
    Ractive = k1 * exp( Ea1 ./ (R * temp(k)) );
    Rpassive = k2 * exp( Ea2 ./ (R * temp(k)) );
    
    dsoot = ConcSootUs(k) - Vsoot * soot_M * (Ractive  * dataStruct.ConcOxyUs(k) + ...
                                                         Rpassive * dataStruct.ConcNoUs(k));
    mSoot = mSoot + timeStep * dsoot / V_total;
    sootload(k) =mSoot;% max(mSoot, 1e-8);
end


end
% mSootOut = (MfWall*ConcExhGasUs*(1/NM.VolCell) -...
%     0.012*NM.GSAic*(RSootCake.R1Sootcake + RDFL.R1DFL +...
%     RSootCake.R2Sootcake + RDFL.R2DFL + RSootCake.R3Sootcake + RDFL.R3DFL +...
%     RSootCake.R4Sootcake + RDFL.R4DFL))*p.DpfWM_Ts_C + mSootDly;

