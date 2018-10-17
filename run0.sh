# if [ -f "log" ];
# then
# rm log
# fi
#
# if [ -f "dump" ];
# then
# rm dump
# fi
#
# if [ -d "intermediate" ]
# then
# rm -r intermediate
# fi
#
# mkdir intermediate

qcc -O2 -Wall $file.c -o $file -lm
