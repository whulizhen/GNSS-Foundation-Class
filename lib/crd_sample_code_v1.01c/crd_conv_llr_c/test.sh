cd test_files
../cllr_to_crd -i moon_sum.9001 -o 7080_apollo15_1990_00.npt
../cllr_to_crd -i s25y07d333t0828#103.lu -o 7080_apollo15_20071129_08_00.crd
echo comparing .npt
diff 7080_apollo15_1990_00.npt 7080_apollo15_1990_00.npt.ref
echo comparing .crd
diff 7080_apollo15_20071129_08_00.crd 7080_apollo15_20071129_08_00.crd.ref

