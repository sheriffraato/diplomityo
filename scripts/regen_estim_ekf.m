
load("tools/ValidationData450230026_39-281_pressures.mat")
data=rmmissing(data);
 data = data(1:100000,:);
mf = data.mf_exhaust;
exhtemp = data.T_DOC_ds_corr;
pout = data.p_exhaust_BP/10-data.p_dP_DPF_offsetCorr/10; %DPF outlet pressure as absolute pressure [kPa]

dataStruct.ExhMassFlow = mf;
dataStruct.ConcC3H6Us = data.conc_C3H6_ERE;
dataStruct.ConcNoUs = data.conc_NOx_ERE - data.conc_NO2_DOCds_mdl;
dataStruct.ConcNo2Us = data.conc_NO2_DOCds_mdl;
dataStruct.ConcOxyUs = data.conc_NOXus_O2*1e4;
dataStruct.ConcCo2Us = data.conc_CO_ERE;
dataStruct.ConcCoUs = data.conc_CO_ERE;
dataStruct.DpfTemp = exhtemp;
dataStruct.ConcSootUs = data.conc_Soot_ERE;
% dataStruct.Ashload = data.m_ash_mdl;
dataStruct.Time = data.t_EngineHours_sec;
timeStep = 1;

% DPF model parameters (DpfWM_params used as default if paramStruct is not
% given to function)
%paramStruct = DpfWM_params();
paramStruct = DpfWM_params_HAC();

ws=0.241/25.4*1000; % thickness of filter wall [mil]
cpsi=325; % DPF channel density [cpsi]
R_alpha = 1.33/1.01; % alphaIn/alphaOut [-] 

diff_ash = [0; diff(data.m_ash_mdl)]*1000;


x_init = [data.conc_soot_mdl(1); data.m_ash_mdl(1)*1000];
P_init = [1 0; 0 1];

% Process and meas covariances. Make something up here
Q = [1 0; 0 0.0001];
R = 0.01;

N = height(data);

x = zeros(2, N);
P = zeros(4, N);


x(:,1) = x_init;
P(:,1) = reshape(P_init,4,1);

msoot = x(1,1);
mash = x(2,1);


[dsoot_in, dsoot_oxi] = soot_jacobian(dataStruct, DPFDiam, DPFLen, paramStruct);
z_meas = data.p_dP_DPF_offsetCorr/10; % kPa
h = zeros(N,1);
%%
tic
for k = 2:N
    
    Fk = diag([dsoot_oxi(k-1), 1]);

    dsoot = dsoot_in(k-1) + dsoot_oxi(k-1).* msoot;
    dash = diff_ash(k);
    x_priori = x(:,k-1) + timeStep*[dsoot; dash];
    
    Qk_minus_1 = Q;
    Pk_minus_1 = reshape(P(:,k-1), 2,2);
    P_priori = Fk*Pk_minus_1*Fk' + Qk_minus_1;
    
    h(k) = deltaP_model_og(DPFDiam,DPFLen,mf(k),exhtemp(k),pout(k),x_priori(1),x_priori(2),ws,cpsi,R_alpha);
    if imag(h(k)) ~= 0
        disp(k)
        break
    end
    z_tilde = z_meas(k)-h(k);
    Hk = deltaP_jacobian(DPFDiam,DPFLen,mf(k),exhtemp(k),pout(k),x_priori(1),x_priori(2),ws,cpsi,R_alpha);
    Rk = R;
    Sk = Hk*P_priori*Hk' + Rk;
    
    Kk = P_priori * Hk' * inv(Sk);

    x_post = x_priori + Kk*z_tilde;
    x_post = max(x_post, 0);
    P_post = (eye(2) - Kk*Hk)*P_priori;
    
    P_post(1,2)=0;
    P_post(2,1)=0;
    P_post = 1/2*(P_post + P_post'); % Force symmetry


    x(:,k) = x_post;
    P(:, k) = reshape(P_post,4,1);

    msoot = x_post(1,1);
    mash = x_post(2,1);
end
toc


%%
figure
hold on
plot(z_meas)
plot(h)

figure 
plot(z_meas- h)

figure
hold on
plot(x(1,:))
plot(data.conc_soot_mdl)
plot(data.conc_DpfClogDiag_SootLoadAvg)
legend

