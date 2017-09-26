#!/bin/bash


#===============================================================================
#  Simple script to kill executions
#===============================================================================

# This variable is key to defining your output
#VERSION="81rk4"

VERION=$1


echo ""
echo "==========================================="

if [ -z "$2" ]
  then
    echo "Two arguments are expected"
    echo "arg 1: string defining the version"
    echo "arg 2: string defining the site or all"
    echo ""
    echo "example 1:"
    echo "prompt>./kill_test.sh 81rk4rapid all"
    echo ""
    echo "example 2:"
    echo "prompt>./kill_test.sh 82rk5long m34"
    echo "==========================================="
    echo ""
    exit
fi

echo "Killing jobs tagged: "$1"for user "${USER}

if [ $2 == "all" ];then

    echo $1
    squeue -u ${USER} -o %20j%4t%10A
    squeue -u ${USER} -o %20j%4t%10A | grep -i $1 > killfiles.txt
    njobs=`wc -l killfiles.txt | awk '{print $1}'`
    echo "Killing "$njobs" jobs"
    echo ""

    if [ $njobs -gt 0 ];then

    while IFS= read -r line
    do
      pid=`echo $line | awk '{print $3}'`
      pname=`echo $line | awk '{print $1}'`
      echo "killing pid: "$pid"  AKA: "$pname
      scancel $pid
    done < killfiles.txt
    rm killfiles.txt
    else
	echo "No Jobs Found"
    fi

else

    squeue -u ${USER} -o %20j%4t%10A | grep -i $1 | grep -i $2 > killfiles.txt
    njobs=`wc -l killfiles.txt | awk '{print $1}'`
    echo "Killing "$njobs" jobs"
    echo ""
    if [ $njobs -gt 0 ];then
    while IFS= read -r line
    do
      pid=`echo $line | awk '{print $3}'`
      pname=`echo $line | awk '{print $1}'`
      echo "killing pid: "$pid"  AKA: "$pname
      scancel $pid
    done < killfiles.txt
    rm killfiles.txt
    fi
fi

echo ""
echo "======================================"
exit