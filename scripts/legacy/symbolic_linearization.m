%% Layers
clc; clear
syms mash msoot ain nopen L ro_soot ro_ash positive


wash = 1/2* (ain - sqrt(ain^2- mash/(nopen*L*ro_ash)));
wsoot = 1/2* (ain - 2*wash - sqrt((ain-2*wash)^2- msoot/(nopen*L*ro_soot)));

pretty(wash)
pretty(wsoot)

%%
dwash_msoot = diff(wash, msoot);
dwash_mash = diff(wash,mash);

dwsoot_mash = diff(wsoot,mash);
dwsoot_msoot= diff(wsoot,msoot);

pretty(dwsoot_msoot)
pretty(dwsoot_mash)

pretty(dwash_msoot)
pretty(dwash_mash)
%%
particulateJacobian = [dwsoot_msoot, dwsoot_mash;
                       dwash_msoot,  dwash_mash];

pretty(particulateJacobian)

%% 
syms mu aout F Q V ws kash ksoot q a d positive 
syms wsoot wash positive
clc

ast = ain -2*wsoot - 2*wash;
% PDin = (ain + aout + 2*ws)^2 * F * L^2 * Q * mu/ (6 * V * ast^4);
% PDsoot = Q*mu / (8 * nopen * L * ksoot) * log( (aout - 2*wash)/(aout - 2*wash - 2*wsoot) );
% PDash  = Q*mu / (8 * nopen * L * kash) * log( (aout)/(aout - 2*wash ) );
PDin = (d)^2 * F * L^2 * q/ (3 * V * ast^4);
PDsoot =  q / (a* ksoot) * log( (aout - 2*wash)/(aout - 2*wash - 2*wsoot) );
PDash  = q / (a*kash) * log( (aout)/(aout - 2*wash ) );

pretty(PDin)
pretty(PDsoot)
pretty(PDash)

%%
clc
dPDin_wsoot = diff(PDin, wsoot);
dPDin_wash = diff(PDin, wash);
pretty(dPDin_wsoot)
pretty(dPDin_wash)
%%
clc
dPDsoot_wsoot = diff(PDsoot, wsoot);
dPDsoot_wash = diff(PDsoot, wash);
pretty(dPDsoot_wsoot)
pretty(dPDsoot_wash)

%%
clc
dPDash_wsoot = diff(PDash, wsoot);
dPDash_wash = diff(PDash, wash);
pretty(dPDash_wsoot)
pretty(dPDash_wash)

%%
clc
pretty(dPDin_wsoot+dPDsoot_wsoot+dPDash_wsoot)
