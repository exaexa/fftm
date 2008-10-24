#!/bin/bash
make

a=''
b=''
program='./fftm'

step=0


echo
echo "Classical multiply demo NOW"
echo
sleep 1

while [ $step -lt 10 ] ; do
	a=$a$RANDOM
	b=$b$RANDOM
	echo "$a * $b = `echo -e $a '\n' $b | $program`"
	let ++step
done

sleep 1
echo
echo
echo "Now trying to compute some powers of 2"
echo
sleep 1

power=1
step=0
while [ $step -le 512 ] ; do
	echo "2^$step = $power"
	power=`echo -e $power '\n' 2 | $program`
	let ++step
done


sleep 1
echo
echo
echo "Now trying to compute some factorials"
echo
sleep 1

a=1
i=1
while [ $i -le 1000 ] ; do
	a=`echo -e $a '\n' $i | $program`
	echo "$i! = $a"
	let ++i
done


