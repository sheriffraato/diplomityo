%%
% Compute the Jacobian Jw of w = [wsoot, wash] with respect to
% m = [msoot, mash]. 
clc; clear
syms mash msoot alphaIn nopen L ro_soot ro_ash positive

wash = 1/2 .* (alphaIn - sqrt(alphaIn.^2 - mash./(nopen .* L .* ro_ash)));
wsoot = 1/2 .* (alphaIn - 2.*wash - sqrt((alphaIn - 2.*wash).^2 - msoot./(nopen .* L .* ro_soot)));

dwsoot_mash = diff(wsoot,mash);
dwsoot_msoot= diff(wsoot,msoot);

dwash_msoot = diff(wash, msoot);
dwash_mash = diff(wash,mash);

% Jw_m = [dwsoot_msoot, dwsoot_mash;
%        dwash_msoot,  dwash_mash];
Jw_m = [dwash_mash, dwash_msoot;
        dwsoot_mash, dwsoot_msoot];

pretty(Jw_m)

%% 
% Compute the Jacobian JdP of deltaP with respect to w.
syms mu alphaOut F Q V ws kash ksoot kwall q a d positive 
syms wsoot wash positive
clc

ast = alphaIn -2*wsoot - 2*wash;
% PDin = (alphaIn + alphaOut + 2*ws)^2 * F * L^2 * Q * mu/ (6 * V * ast^4);
% PDsoot = Q*mu / (8 * nopen * L * ksoot) * log( (alphaOut - 2*wash)/(alphaOut - 2*wash - 2*wsoot) );
% PDash  = Q*mu / (8 * nopen * L * kash) * log( (alphaOut)/(alphaOut - 2*wash ) );
% PDwall = ws / (a * alphaOut * kwall) * Q*mu;
% PDout = (alphaIn + alphaOut + 2*ws)^2 * F * L^2 * Q * mu/ (6 * V * alphaOut^4);

cf = Q .* mu .* F .* L.^2 ./ (6 * V);
PDin = cf .* ((alphaIn + alphaOut + 2 .* ws).^2 ./ (alphaIn -2*wash -2*wsoot ).^4);
PDout = cf .* ((alphaIn + alphaOut + 2 .* ws).^2 ./ alphaOut.^4);
PDwall = Q .* mu .* ws ./ (4 .* nopen .* L .* alphaOut .* kwall);
PDash = Q.*mu ./ (8 * nopen .* L .* kash) .* log(alphaOut ./ (alphaOut - 2.*wash));
PDsoot = Q.*mu ./ (8 * nopen .* L .* ksoot) .* log( (alphaOut - 2.*wash) ./ (alphaOut - 2.*wash - 2.*wsoot));


PD = PDin + PDsoot + PDash + PDout + PDwall;
PD = simplify(PD);

JdP_w = [diff(PD, wsoot), diff(PD, wash)];

pretty(JdP_w)

%%
% The Jacobian of dP with respect to m is given by the chain rule.
clc
JdP_m = JdP_w*Jw_m;
JdP_m = simplify(JdP_m);
% pretty(JdP_m)

pretty(simplify(JdP_m(1)))
disp('AND \n')
pretty(JdP_m(2))

