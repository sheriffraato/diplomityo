load("tools/ValidationData450230026_39-281_pressures.mat")
data=rmmissing(data);

dataStruct.ExhMassFlow = data.mf_exhaust;
dataStruct.ConcC3H6Us = data.conc_C3H6_ERE;
dataStruct.ConcNoUs = data.conc_NOx_ERE - data.conc_NO2_DOCds_mdl;
dataStruct.ConcNo2Us = data.conc_NO2_DOCds_mdl;
dataStruct.ConcOxyUs = data.conc_NOXus_O2*1e4;
dataStruct.ConcCo2Us = data.conc_CO_ERE;
dataStruct.ConcCoUs = data.conc_CO_ERE;
dataStruct.DpfTemp = data.T_DOC_ds_corr;
dataStruct.ConcSootUs = data.conc_Soot_ERE;
dataStruct.Ashload = data.m_ash_mdl;
dataStruct.Time = data.t_EngineHours_sec;

% DPF model parameters (DpfWM_params used as default if paramStruct is not
% given to function)
%paramStruct = DpfWM_params();
paramStruct = DpfWM_params_HAC();

% Initial soot load
sootLoadInit = data.conc_soot_mdl(1);
timeStep = 1;

tic
sootload_THa = regen_model(dataStruct, DPFDiam, DPFLen, sootLoadInit, timeStep);
toc

%%
tic
N = height(data);
msoot = sootLoadInit;
sootload_new = ones(N,1)*msoot;
[dsoot_in, dsoot_oxi] = soot_jacobian(dataStruct, DPFDiam, DPFLen, paramStruct);

for k = 2:N
    
    dsoot_oxi(k) = dsoot_oxi(k).*msoot;

    dsoot = dsoot_in(k) + dsoot_oxi(k);
    sootload_new(k) = msoot + timeStep*dsoot;
    msoot = sootload_new(k);
end
toc
%%



figure('Name', 'Soot loading')
hold on
% plot(dataStruct.Time, sootLoad, 'DisplayName', 'Original')
plot(dataStruct.Time, sootload_THa, 'DisplayName','THa')
plot(dataStruct.Time, sootload_new, 'DisplayName', 'New')
%plot(dataStruct.Time, data.conc_DpfClogDiag_SootLoadAvg, 'DisplayName', 'SootLoadAvgDV')
% plot(dataStruct.Time, data.conc_soot_mdl, 'DisplayName', 'SootMdlDV')
plot(DPFWeighingEngineHours, DPFWeighingSootLoad, 'LineStyle','none', 'Marker','*', 'DisplayName','Weighings')
legend('Location','northwest')
hold off
ylabel('Soot loading [g/l]')
xlabel('Time [s]')