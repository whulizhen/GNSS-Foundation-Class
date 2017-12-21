#!/bin/sh
# Test crd_chk
cd test_files
../crd_chk 7080_lageos1_crd_20071202_16_00.npt > 7080_lageos1_crd_20071202_16_00.npt.chk
../crd_chk 7080_lageos1_crd_20071202_16_00.qlk > 7080_lageos1_crd_20071202_16_00.qlk.chk
../crd_chk 7080_lageos1_crd_20071202_16_00.frd > 7080_lageos1_crd_20071202_16_00.frd.chk
../crd_chk beaconc_20110111.npt > beaconc_20110111.npt.chk
../crd_cstg_np_cmp -1 7080_etalon2_crd_20080506_10_00.npt \
-2 7080_etalon2_crd_20080506_10_00.qld \
-o 7080_etalon2_crd_20080506_10_00.cmp
../crd_merit_fr_cmp -1 7080_etalon2_crd_20080506_10_00.frd \
-2 7080_etalon2_crd_20080506_10_00.mrt \
-o 7080_etalon2_crd_20080506_10_00.frcmp
echo Comparing 7080_lageos1_crd_20071202_16_00.npt.chk
diff 7080_lageos1_crd_20071202_16_00.npt.chk 7080_lageos1_crd_20071202_16_00.npt.chk.ref
echo Comparing 7080_lageos1_crd_20071202_16_00.qlk.chk
diff 7080_lageos1_crd_20071202_16_00.qlk.chk 7080_lageos1_crd_20071202_16_00.qlk.chk.ref
echo Comparing 7080_lageos1_crd_20071202_16_00.frd.chk
diff 7080_lageos1_crd_20071202_16_00.frd.chk 7080_lageos1_crd_20071202_16_00.frd.chk.ref
echo Comparing beaconc_20110111.npt.chk
diff beaconc_20110111.npt.chk beaconc_20110111.npt.chk.ref
echo Comparing 7080_etalon2_crd_20080506_10_00.cmp
diff 7080_etalon2_crd_20080506_10_00.cmp 7080_etalon2_crd_20080506_10_00.cmp.ref
echo Comparing 7080_etalon2_crd_20080506_10_00.frcmp
diff 7080_etalon2_crd_20080506_10_00.frcmp 7080_etalon2_crd_20080506_10_00.frcmp.ref
