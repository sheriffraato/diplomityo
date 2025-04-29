

exhtemp = data.T_DOC_ds;
mf = data.mf_exhaust;
pout = data.p_exhaust_BP/10-data.p_dP_DPF_offsetCorr/10; %DPF outlet pressure as absolute pressure [kPa]
sootload = data.conc_soot_mdl;
ashload = data.mf_exhaust;

ws=0.241/25.4*1000; % thickness of filter wall [mil]
cpsi=325; % DPF channel density [cpsi]
R_alpha = 1.33/1.01; % alphaIn/alphaOut [-] 

deltaP_THa = deltaP_model(DPFDiam, DPFLen,mf,exhtemp,pout,sootload,ashload, ws,cpsi, R_alpha);
deltaP_orig = HacDpfPressDropStageV9_300_C640_v20241010(DPFDiam, DPFLen,mf,exhtemp,pout,sootload,ashload, ws,cpsi, R_alpha);
%%
t = data.t_EngineHours_sec/3600;
figure
hold on                                 
plot(t, deltaP_orig, 'LineWidth',1,'DisplayName','Original')
plot(t, deltaP_THa, 'DisplayName','THa')
legend

norm(deltaP_THa-deltaP_orig)