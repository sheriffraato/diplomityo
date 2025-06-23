function params = DpfWM_params_HAC()
%DPFWM_PARAMS returns a struct that contains parameters of the DpfWM-model.

%Outputs:
%params: struct containing the parameters

params.DpfWM_cdr_P = 1.2803;
params.DpfWM_ExpnDflNo2AcvnERgn_P = 10792;
params.DpfWM_ExpnDflNo2AcvnERgnInh1_P = 0.5;
params.DpfWM_ExpnDflNo2AcvnERgnInh2_P = 3;
params.DpfWM_ExpnDflOxyAcvnERgn_P = 15000;
params.DpfWM_ExpnDflOxyAcvnERgnInh1_P = 1000;
params.DpfWM_ExpnDflOxyAcvnERgnInh2_P = 3;
params.DpfWM_ExpnNoxAcvnEOxd1_P = 47191;
params.DpfWM_ExpnNoxAcvnEOxd2_P = -2500;
params.DpfWM_ExpnNoxAcvnEOxd3_P = -1500;
params.DpfWM_ExpnNoxAcvnEOxd4_P = -2000;
params.DpfWM_ExpnSootCakeNo2AcvnERgn_P = 10792; 
params.DpfWM_ExpnSootCakeNo2AcvnERgnInh1_P = 0.5;
params.DpfWM_ExpnSootCakeNo2AcvnERgnInh2_P = 1;
params.DpfWM_ExpnSootCakeOxyAcvnERgn_P = 24654.699;
params.DpfWM_ExpnSootCakeOxyAcvnERgnInh1_P = 5000;
params.DpfWM_ExpnSootCakeOxyAcvnERgnInh2_P = 1;
params.DpfWM_GainDflNo2AcvnERgn_P = 3012830538;
params.DpfWM_GainDflNo2AcvnERgnInh_P = 44.4275;
params.DpfWM_GainDflOxyAcvnERgn_P = 547256;
params.DpfWM_GainDflOxyAcvnERgnInh_P = 540000;
params.DpfWM_GainNoxArfacOxd1_P = 2714658;
params.DpfWM_GainNoxArfacOxd2_P = 0;
params.DpfWM_GainNoxArfacOxd3_P = 0;
params.DpfWM_GainNoxArfacOxd4_P = 0;
params.DpfWM_GainSootCakeNo2AcvnERgn_P = 3012830538;
params.DpfWM_GainSootCakeNo2AcvnERgnInh_P = 44.4275;
params.DpfWM_GainSootCakeOxyAcvnERgn_P = 2072679936;
params.DpfWM_GainSootCakeOxyAcvnERgnInh_P = 540000;
params.DpfWM_m_P = 1;
params.DpfWM_MaxSootLoadingDFL_P = 0.75;
params.DpfWM_MmolOfAir_P = 2.886000e+01;
params.DpfWM_n_P = 1.700000e+00;
params.DpfWM_Rgas_P = 8.314;
params.DpfWM_RhoCell_P = 300;
params.DpfWM_RhoSoot_P = 100;
params.DpfWM_ThcknsWall_P = 2.2860e-4;
x = [0.087004958226877,0.836864352254616];
params.DpfWM_GainSootCakeNo2AcvnERgn_P = x(1)*params.DpfWM_GainSootCakeNo2AcvnERgn_P;
params.DpfWM_ExpnSootCakeNo2AcvnERgn_P = x(2)*params.DpfWM_ExpnSootCakeNo2AcvnERgn_P;
params.DpfWM_GainDflNo2AcvnERgn_P = params.DpfWM_GainSootCakeNo2AcvnERgn_P;
params.DpfWM_ExpnDflNo2AcvnERgn_P = params.DpfWM_ExpnSootCakeNo2AcvnERgn_P;

end