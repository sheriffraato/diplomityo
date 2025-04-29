%% DPF simulation template
% Template for simulating DPF soot loading with Matlab-function
% DPFWall_timeseries

% Ari-Pekka Kinnunen, ari-pekka.kinnunen@agcocorp.com
% 14.10.2022

%% Define input signal structure and parameters

% DPF dimensions
dpfLen = 5;
dpfDiam = 9.5;

load('data_nrtc_20220422.mat')

% Input data structure
dataStruct.ExhMassFlow = data.EreI_MfExh;
dataStruct.ConcC3H6Us = data.DocM_ConcC3H6Ds;
dataStruct.ConcNoUs = data.DocM_ConcNoDs;
dataStruct.ConcNo2Us = data.DocM_ConcNo2Ds;
dataStruct.ConcOxyUs = data.DocM_ConcOxyDs;
dataStruct.ConcCo2Us = data.DocM_ConcCo2Ds;
dataStruct.ConcCoUs = data.DocM_ConcCoDs;
dataStruct.DpfTemp = data.ZIZU_SENSOR_DOCDSTEMP_value_mp;
dataStruct.ConcSootUs = data.EreSAM_MconcSoot;
dataStruct.Time = seconds(data.Time-data.Time(1));

%% Define input signal structure and parameters: THa; 450.23.0026
load("masters-thesis-main\code\ValidationData450230026.mat")
dpfLen = DPFLen;
dpfDiam = DPFDiam;

dataStruct.ExhMassFlow = data.mf_exhaust;
dataStruct.ConcC3H6Us = data.conc_C3H6_ERE;
dataStruct.ConcNoUs = data.conc_NOx_ERE - data.conc_NO2_DOCds_mdl;
dataStruct.ConcNo2Us = data.conc_NO2_DOCds_mdl;
dataStruct.ConcOxyUs = data.conc_NOXus_O2;
dataStruct.ConcCo2Us = data.conc_CO_ERE;
dataStruct.ConcCoUs = data.conc_CO_ERE;
dataStruct.DpfTemp = data.T_DOC_ds_corr;
dataStruct.ConcSootUs = data.conc_Soot_ERE;
dataStruct.Time = data.t_EngineHours_sec;

% DPF model parameters (DpfWM_params used as default if paramStruct is not
% given to function)
%paramStruct = DpfWM_params();
paramStruct = DpfWM_params_HAC();

% Initial soot load
sootLoadInit = 0.01;

%% Run simulation

timeStep = 0.1;

[sootLoad, R_SootCake, R_DFL] = DPFWall_timeseries(dataStruct, dpfDiam, dpfLen, sootLoadInit, timeStep, paramStruct);

%% Plot results


figure('Name', 'Soot loading')
plot(dataStruct.Time, sootLoad)
ylabel('Soot loading [g/l]')
xlabel('Time [s]')
    
figure('Name', 'Reaction rates')
tiledlayout(1,2)
nexttile
plot(dataStruct.Time, R_SootCake)
xlabel('Time [s]')
ylabel('r_{cake}')
legend('r_1', 'r_2', 'r_3', 'r_4')
nexttile
plot(dataStruct.Time, R_DFL)
xlabel('Time [s]')
ylabel('r_{DFL}')
legend('r_1', 'r_2', 'r_3', 'r_4')

