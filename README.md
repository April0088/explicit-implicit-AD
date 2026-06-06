# explicit-implicit-AD
This repository contains the code for the study: Dissociation between impaired explicit spatial remapping and preserved implicit neural dynamics in Alzheimer’s disease


pre-process: Data_geteventTS.m; Data_getvideo_ABBA.m; code_Process_Data_video.m
Get ratemap: do_all_ratemap_OF.m
Plot ratemap: do_all_plot_ratemap_v2.m
Get features related to remapping: do_all_remapping_ABBA_v2.m
Get tau of cell pairs: Pcorr_tau.m
Do SVM for tau of cell pairs: do_all_get_tau_for_svm.m ; tau_SVM_v6.m
Get position-tuning independent rate (PIR): do_posi_indep_rate.m
Get tau of PIR of cell pairs: Pcorr_tau_PIR.m
Get CCG for theta modulation: cell_pair_space_overlap_CCG.m
Detecte SWRs during awake rest: do_all_detect_SWR_OF_rest.m
Get reactivation in SWRs: do_all_spikesinSWR.m
