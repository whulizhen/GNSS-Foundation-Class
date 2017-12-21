#!/bin/sh
# Test crd sample code 
# Now go through all directories...
for i in `ls`
  do
    if test -d $i
    then
      cd $i
      if test -e NO_TEST
      then
        echo "Bypassing $i"
      else
        echo "Doing $i"
	./test.sh
        echo $1
        echo "-c"
##	if test -z $1
##	then
##          make all
##	fi
      fi
      cd ..
    fi
  done
echo done
