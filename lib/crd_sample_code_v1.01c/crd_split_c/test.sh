cd test_files
rm -f *.npt *.qlk *.frd
../crd_split moon.9312.crd
../crd_split 7080_crd_20071220_00.crd
../frd_strip s25y07d336t1600#1155.frfin 7080_lageos1_20071202_16_00.frd
for i in `ls *.npt`
do
  echo comparing $i
  diff $i $i.ref
done
for i in `ls *.qlk`
do
  echo comparing $i
  diff $i $i.ref
done
for i in `ls *.frd`
do
  echo comparing $i
  diff $i $i.ref
done
cd ..
./merge_crd_daily
cd merge_test_files/crd_merge
echo Differences on h1 headers are due to different creation dates.
for i in `ls *.npt`
do
  echo comparing $i
  diff $i $i.ref
done
for i in `ls *.qlk`
do
  echo comparing $i
  diff $i $i.ref
done
for i in `ls *.frd`
do
  echo comparing $i
  diff $i $i.ref
done

