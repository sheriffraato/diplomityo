function deltaP_total=deltaP_model(Dtrap,Ltrap,mf,exhtemp,pout,sootload,ashload,ws,cpsi,R_alpha)
%syms Dtrap Ltrap mf exhaustTemp pressureOut sootload ashload ws cpsi R_alpha

 
% DPF-dimensions
inch2meter = 0.0254;  
rim = 0.003875; % [m]
plugdepth=0.005; % [m]
D=Dtrap*0.0254-2*rim; %inch->m
L=Ltrap*0.0254-2*plugdepth; %inch->m
V=pi*D^2/4*L; %filter total volume m3

% filter-properties
ws=ws*2.54e-5; %mil->m
cpsi=cpsi/0.0254^2; %1/inch^2 -> 1/m2
nopen=cpsi*pi*D^2/4/2;
alphaMean=sqrt(1/cpsi)-ws; %m mean side lenght of the filter opening
alphaIn = 2*alphaMean/(1+1/R_alpha); %m side length of the filter inlet opening 
alphaOut = 2*alphaMean/(1+R_alpha); %m side length of the filter inlet opening 

% Fitted parameters
%kwall = 3.189e-13;% permeability of filter wall [m2]
kwall=3.4e-13;% permeability of filter wall [m2]
ksoot=4.399e-14;% permeability of soot layer [m2]
ksoot=9.14130161581979e-14;% permeability of soot layer [m2]

kash=1.3204e-09; % permeability of ash layer [m2]
kash=1.65e-13; % permeability of ash layer [m2]. 

ro_soot=34.633; % density of the soot layer [kg/m3]
ro_ash=5.1832e+05; % density of ash [g/m3]
F = 28.454;


%% Gas properties
mf=mf/(60*60); %mass flow kg/s
mu = viscosityD(exhtemp); 
ro = gasdensity(exhtemp, pout);
Q = mf ./ ro;
channel_factor = (Q .* mu .* F / 6) .* (L.^2 ./ V);
%%
% Particulate layers
msoot = sootload.*V; % kg
wash = 1/2 .* (alphaIn - sqrt(alphaIn.^2 - ashload./(nopen .* L .* ro_ash)));
wsoot = 1/2 .* (alphaIn - 2.*wash - sqrt((alphaIn - 2.*wash).^2 - msoot./(nopen .* L .* ro_soot)));

% Channel pressure drops
deltaP_inlet = channel_factor .* ((alphaIn + alphaOut + 2 .* ws).^2 ./ (alphaIn -2*wash -2*wsoot ).^4);
deltaP_outlet = channel_factor .* ((alphaIn + alphaOut + 2 .* ws).^2 ./ alphaOut.^4);

% Wall pressure drop
deltaP_wall = Q .* mu .* ws ./ (4 .* nopen .* L .* alphaOut .* kwall);

% Particulate cake pressure drops
deltaP_ash = Q.*mu ./ (8 * nopen .* L .* kash) .* log(alphaOut ./ (alphaOut - 2.*wash));
deltaP_soot = Q.*mu ./ (8 * nopen .* L .* ksoot) .* log( (alphaOut - 2.*wash) ./ (alphaOut - 2.*wash - 2.*wsoot));

% Total pressure drop of the DPF
deltaP_total = deltaP_inlet + deltaP_outlet + deltaP_wall + deltaP_soot + deltaP_ash;
deltaP_total = deltaP_total/1000; % Pa -> kPa

end
%

function ro=gasdensity(T,p,M)
    %calculates ideal gas density under given temperature and pressure
    %inputs:
    %M = molecular weight [g/mol]
    %T = temperature [degC]
    %p = pressure [kPa]
    %
    %outputs:
    %ro = gas density [kg/m3]
    if nargin<3
        M=1.030388422.*28.02; %g/mol (assumed 78 v-% N2, 14 v-% O2, 4 v-% CO2, 4 v-% H2O
    end
    M=M./1000; %kg/mol
    R=8.314462; % J/(mol K) kaasuvakio
    T=T+273.15; %C->K
    p=p.*1000; %kPa->Pa
    
    ro = (p.*M)./(R.*T);
end

function mu=viscosityD(T)
    %mu = viscosity (Pa s)
    %T = temperature (C)
    T=T+273.15; %C->K
    % [T(K) mu e-5(kg/(ms)) v e-6 (m2/s)]
    %from engineering toolbox
    % http://www.engineeringtoolbox.com/dry-air-properties-d_973.html
    M=[ 100	    0.6924  1.923
        150	    1.0283	4.343
        200	    1.3289	7.490
        250	    1.488	9.49    
        300	    1.983	15.68
        350	    2.075	20.76
        400	    2.286	25.90
        450	    2.484	28.86
        500	    2.671	37.90
        550	    2.848	44.34
        600	    3.018	51.34
        650	    3.177	58.51
        700	    3.332	66.25
        750	    3.481	73.91
        800	    3.625	82.29
        850	    3.765	90.75
        900	    3.899	99.30
        950	    4.023	108.2
        1000	4.152	117.8
        1100	4.44	138.6
        1200	4.69	159.1
        1300	4.93	182.1
        1400	5.17	205.5
        1500	5.40	229.1
        1600	5.63	254.5];
    mu=interp1(M(:,1),M(:,2),T);
    mu=mu.*1e-5;
end