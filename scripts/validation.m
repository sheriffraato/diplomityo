load("AGCO-files\tools\ValidationData450230026_39-281_pressures.mat")
cutData = data(:,:);

%% SOOT VALIDATION

ws=0.241/25.4*1000; % thickness of filter wall [mil]
cpsi=325; % DPF channel density [cpsi]
R_alpha = 1.33/1.01; % alphaIn/alphaOut [-] 

exhtemp = cutData.T_DOC_ds_corr;
mf = cutData.mf_exhaust;
pout = cutData.p_exhaust_BP/10-cutData.p_dP_DPF_offsetCorr/10; %DPF outlet pressure as absolute pressure [kPa]

exhtemp = ones(size(cutData.T_DOC_ds_corr))*mean(exhtemp);
mf = ones(size(cutData.T_DOC_ds_corr))*mean(mf);
pout=ones(size(cutData.T_DOC_ds_corr))*mean(pout);

dataStruct.ExhMassFlow = mf;
dataStruct.DpfTemp = exhtemp;         % degC 
dataStruct.OutletPress = pout;

sootload = cutData.conc_soot_mdl ;
ashload = cutData.m_ash_mdl*1000*0;

deltaP_orig = HacDpfPressDropStageV9_300_C640_v20241010(DPFDiam, DPFLen,mf,exhtemp,pout,sootload,ashload, ws,cpsi, R_alpha);
deltaP_THa = deltaP_model(dataStruct, DPFDiam,DPFLen,sootload,ashload,ws,cpsi,R_alpha);
Jdp_og = deltaP_jacobian(DPFDiam,DPFLen,mf,exhtemp,pout,sootload,ashload,ws,cpsi,R_alpha);
Jdp_tot = get_deltaP_jacobian2(dataStruct, DPFDiam,DPFLen,sootload,ashload,ws,cpsi,R_alpha);

%%
t = cutData.t_EngineHours_sec/3600;
figure
hold on                                 
plot(t, deltaP_orig, 'LineWidth',1,'DisplayName','Original mdl')
% plot(t, deltaP_THa_og, 'DisplayName','THa mdl')
legend
xlabel('engine hours [h]')
ylabel('DeltaP [kPa]')
title('Model comparison')

% norm(deltaP_THa-deltaP_orig)


%%

% DPF-dimensions used to compute sootload --> msoot
rim = 0.003875; % [m]
plugdepth=0.005; % [m]
D=DPFDiam*0.0254-2*rim; %inch->m
L=DPFLen*0.0254-2*plugdepth; %inch->m
V=pi*D^2/4*L; 

for m_e = 0:0.5:3.5
    x = find(abs(sootload - m_e) < 1e-2, 1);   
    if isempty(x)
        disp('no feasible x')
        continue
    end
    
    Jdp = Jdp_og(x,:);
    S = sootload((m_e-1<sootload) & (m_e+1>sootload));

    figure  
    hold on
    plot(sootload, deltaP_THa, 'r.', 'DisplayName','DeltaP model')
    
    plot(S, deltaP_THa(x)+Jdp(1)*V * (S-m_e), 'b.', ...
        'DisplayName',['Linearized at ' num2str(m_e) ' g/l'])
    legend(Location="northwest")
    ylabel('DeltaP [kPa]')
    xlabel('sootload [g/l]')
    hold off
    % a = deltaP_THa(x)+Jdp(1)* random_constant *(sootload-m) - deltaP_THa;
    % a = abs(a);
    % sootload(a==min(a))
    
end

%% ASH VALIDATION

sootload = cutData.conc_soot_mdl*0;
ashload = cutData.m_ash_mdl*1000;

deltaP_orig = HacDpfPressDropStageV9_300_C640_v20241010(DPFDiam, DPFLen,mf,exhtemp,pout,sootload,ashload, ws,cpsi, R_alpha);
deltaP_THa = deltaP_model(dataStruct, DPFDiam,DPFLen,sootload,ashload,ws,cpsi,R_alpha);
Jdp_tot = get_deltaP_jacobian2(dataStruct, DPFDiam,DPFLen,sootload,ashload,ws,cpsi,R_alpha);

%%
t = cutData.t_EngineHours_sec/3600;
figure
hold on                                 
plot(t, deltaP_orig, 'LineWidth',1,'DisplayName','Original mdl')
plot(t, deltaP_THa, 'DisplayName','THa mdl')
legend
xlabel('engine hours [h]')
ylabel('DeltaP [kPa]')
title('Model comparison')

%%

% DPF-dimensions used to compute sootload --> msoot
rim = 0.003875; % [m]
plugdepth=0.005; % [m]
D=DPFDiam*0.0254-2*rim; %inch->m
L=DPFLen*0.0254-2*plugdepth; %inch->m
V=pi*D^2/4*L; 

for m_e = 0:6
    x = find(abs(ashload - m_e) < 1e-2, 1);   
    if isempty(x)
        disp('no feasible x')
        continue
    end
    
    % Jdp = get_deltaP_jacobian(DPFDiam, DPFLen,mf(x),exhtemp(x),pout(x),sootload(x),ashload(x), ws,cpsi, R_alpha);
    Jdp = Jdp_tot(x,:);
    A = ashload((m_e-1<ashload) & (m_e+1>ashload));

    figure  
    hold on
    plot(ashload, deltaP_THa, 'r.', 'DisplayName','DeltaP model')
    
    % plot(A, deltaP_THa(x) + Jdp(2) * (A-m_e), 'b.', ...
    %     'DisplayName',['Ver1: Linearized at ' num2str(m_e) ' g/l'])
    plot(A, deltaP_THa(x) + Jdp(2) * (A-m_e), 'k.', ...
        'DisplayName',['Ver2: Linearized at ' num2str(m_e) ' g/l'])
    legend(Location="northwest")
    ylabel('DeltaP [kPa]')
    xlabel('ashload [g]')
    hold off
    % a = deltaP_THa(x)+Jdp(1)* random_constant *(sootload-m) - deltaP_THa;
    % a = abs(a);
    % sootload(a==min(a))
    
end