
function Jdp = get_deltaP_jacobian(DPFDiam, DPFLen,mf,exhtemp,pout,sootload,ashload,ws,cpsi,R_alpha)
% TUHKAMALLI PIELESSÄ!!!!

% DPF-params:
 
% DPF-dimensions
rim = 0.003875; % [m]
plugdepth=0.005; % [m]
D=DPFDiam*0.0254-2*rim; %inch->m
L=DPFLen*0.0254-2*plugdepth; %inch->m
V=pi*D^2/4*L; %filter total volume m3

% filter-properties
ws=ws*2.54e-5; %mil->m
cpsi=cpsi/0.0254^2; %1/inch^2 -> 1/m2
nopen=cpsi*pi*D^2/4/2;
alphaMean=sqrt(1/cpsi)-ws; %m mean side lenght of the filter opening
alphaIn = 2*alphaMean/(1+1/R_alpha); %m side length of the filter inlet opening 
alphaOut = 2*alphaMean/(1+R_alpha); %m side length of the filter inlet opening 

% Constants and fitted parameters
ro_soot=34.633; % density of the soot layer [kg/m3]
ro_ash=5.1832e+05; % density of ash [g/m3] !!!

F = 28.454;
%ksoot=4.399e-14;% permeability of soot layer [m2]
%kash=1.3204e-09; % permeability of ash layer [m2]
ksoot=9.14130161581979e-14;% permeability of soot layer [m2]
kash=1.65e-13; % permeability of ash layer [m2]. 

% Time variant parameters
mf=mf/(60*60); %mass flow kg/s
exhtemp = exhtemp + 0; % temperature degC
pout = pout * 1; % pressure kPa
mu = viscosityD(exhtemp); 
ro = gasdensity(exhtemp, pout);
Q = mf ./ ro;

% Input dependent parameters
msoot = sootload * V; % kg
mash = ashload; % g !!!
wash = 1/2 .* (alphaIn - sqrt(alphaIn.^2 - mash./(nopen .* L .* ro_ash)));
wsoot = 1/2 .* (alphaIn - 2.*wash - sqrt((alphaIn - 2.*wash).^2 - msoot./(nopen .* L .* ro_soot)));


X = max(alphaIn^2 - mash ./ (L .* nopen .* ro_ash), 0);
Y = max(       X - msoot ./ (L .* nopen .* ro_soot), 0);
S = alphaOut - 2*wsoot - 2*wash;
T = alphaIn - 2*wsoot - 2*wash;
U = alphaOut - 2*wash;
M = 4*L*nopen;
P = alphaIn + alphaOut + 2*ws;

A = Q .* mu ./ (M .* ksoot .* S);
B = 4 .* F .* L^2 .* P^2 .* Q .* mu ./ (3 .* V .* T.^5);
C = A + B;
D = M .* ro_soot .* sqrt(Y);
E = (M .* ro_ash .* sqrt(X)).^(-1);
Fv = (M .* ro_ash .* sqrt(Y) ).^(-1);
G = Q .* mu ./ (M .* kash .* U);
H = Q .* mu ./ (M .* ksoot .* U) .* (U./S + 1);

Jdp = [C./D, C.*(E-Fv)-E.*(B-G+H)];
Jdp = Jdp * 1e-3; 

end

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