%%
% Compute the Jacobian Jw of w = [wsoot, wash] with respect to
% m = [msoot, mash]. 
clc; clear
syms mash msoot ain nopen L ro_soot ro_ash positive

wash = 1/2* (ain - sqrt(ain^2- mash/(nopen*L*ro_ash)));
wsoot = 1/2* (ain - 2*wash - sqrt((ain-2*wash)^2- msoot/(nopen*L*ro_soot)));


dwsoot_mash = diff(wsoot,mash);
dwsoot_msoot= diff(wsoot,msoot);

dwash_msoot = diff(wash, msoot);
dwash_mash = diff(wash,mash);

Jw_m = [dwsoot_msoot, dwsoot_mash;
       dwash_msoot,  dwash_mash];

pretty(Jw_m)

%% 
% Compute the Jacobian JdP of deltaP with respect to w.
syms mu aout F Q V ws kash ksoot kwall q a d positive 
syms wsoot wash positive
clc

ast = ain -2*wsoot - 2*wash;
PDin = (ain + aout + 2*ws)^2 * F * L^2 * Q * mu/ (6 * V * ast^4);
PDsoot = Q*mu / (8 * nopen * L * ksoot) * log( (aout - 2*wash)/(aout - 2*wash - 2*wsoot) );
PDash  = Q*mu / (8 * nopen * L * kash) * log( (aout)/(aout - 2*wash ) );

PDwall = ws / (a * aout * kwall) * Q*mu;
PDout = (ain + aout + 2*ws)^2 * F * L^2 * Q * mu/ (6 * V * aout^4);


PD = PDin + PDsoot + PDash + PDout + PDwall;
PD = simplify(PD);

JdP_w = [diff(PD, wsoot), diff(PD, wash)];

pretty(JdP_w)

%%
% The Jacobian of dP with respect to m is given by the chain rule.
clc
JdP_m = JdP_w*Jw_m;
% JdP_m = simplify(JdP_m);
% pretty(JdP_m)

% pretty(JdP_m(1))
pretty(JdP_m(2))

