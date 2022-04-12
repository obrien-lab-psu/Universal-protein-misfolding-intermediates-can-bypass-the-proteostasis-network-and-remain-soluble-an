mol load psf multidom/1cli/output/1_r345_steps131680_1cli.psf dcd multidom/1cli/output/1_r345_trans_term1_1cli.dcd
mol addfile multidom/1cli/output/1_r345_trans_term2_1cli.dcd
mol addfile multidom/1cli/output/1_r345_trans_term3_1cli.dcd
mol addfile multidom/1cli/output/1_r345_trans_term4_1cli.dcd
mol addfile multidom/1cli/output/1_r345_trans_term5_1cli.dcd
mol addfile multidom/1cli/output/1_r345_trans_term6_1cli.dcd
animate write dcd multidom/1cli/output/1_termi_merged_1cli.dcd
exit