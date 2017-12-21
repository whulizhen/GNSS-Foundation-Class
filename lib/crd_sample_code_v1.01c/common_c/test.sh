#!/bin/sh
# Test wr_test and crd_readf/crd_write modules
cd test_files
../wr_test 7080_giovea_crd_20080508_09_00.frd 7080_giovea_crd_20080508_09_00.frd_wr
../wr_test 7080_giovea_crd_20080508_09_00.npt 7080_giovea_crd_20080508_09_00.npt_wr
echo Comparing full rate files
diff 7080_giovea_crd_20080508_09_00.frd_wr 7080_giovea_crd_20080508_09_00.frd_wr.ref
echo comparing normal point files
diff 7080_giovea_crd_20080508_09_00.npt_wr 7080_giovea_crd_20080508_09_00.npt_wr.ref

