#!/bin/bash

#===============================================================================
#  Simple script to check executions
#===============================================================================

# This variable is key to defining your output

echo ""
echo "==========================================="

if [ -z "$1" ]
  then
    echo "One arguments is expected"
    echo "arg 1: string defining the run version"
    echo ""
    echo "example 1:"
    echo "prompt>./check_runs.sh 81rk4rapid"
    echo "==========================================="
    echo ""
    exit
fi

version=$1
folder=$1

nmain=`ls ${folder}/main*out | wc -l`
ntest=`ls ${folder}/test*out | wc -l`
ndbug=`ls ${folder}/dbug*out | wc -l`

ls ${folder}/main*out > mainfiles.txt
ls ${folder}/test*out > testfiles.txt
ls ${folder}/dbug*out > dbugfiles.txt

# Loop through the longest list and pull the tags

if [ $nmain -gt $ntest ]; then
    if [ $nmain -gt $ndbug ]; then
	#USE NMAIN
	LIST=mainfiles.txt
	FLAG=main
    else
	#USE NDBUG
	LIST=dbugfiles.txt
	FLAG=dbug
    fi
else
    if [ $ntest -gt $ndbug ]; then
	#USE NTEST
	LIST=testfiles.txt
	FLAG=test
    else
	#USE NDBUG
	LIST=dbugfiles.txt
	FLAG=dbug
    fi
fi


echo "SITE DBUG TEST MAIN"
while IFS= read -r line
  do

  echo $line

  id1=`expr index "$line" '/'`
  token=${line:$id1}
  id2=`expr index "$token" '_'`
  tag=${token:$id2:3}

  # Check dbug
  DCHECK=`qstat -u ${USER} | grep -is $tag | grep -is dbug | wc -l`
  if [ $DCHECK -gt 0 ]; then
      DRES="RUNN"
  else
      lastfile=`ls ${folder}/dbug_${tag}*out | awk '{ f=$NF };END{ print f }'`
      if [ -f $lastfile ]; then
	  DCHECK=`grep -is 'ED execution ends' $lastfile | wc -l`
	  if [ $DCHECK -gt 0 ];then
	      DRES="COMP"
	  else
	      DRES="FAIL"
	  fi
      else
	  DRES="----"
      fi
  fi

  # Check test
  TCHECK=`qstat -u ${USER} | grep -is $tag | grep -is test | wc -l`
  if [ $TCHECK -gt 0 ]; then
      TRES="RUNN"
  else
      lastfile=`ls ${folder}/test_${tag}*out | awk '{ f=$NF };END{ print f }'`
      if [ -f $lastfile ]; then
          TCHECK=`grep -is 'ED execution ends' $lastfile | wc -l`
          if [ $TCHECK -gt 0 ];then
              TRES="COMP"
          else
              TRES="FAIL"
          fi
      else
          TRES="----"
      fi
  fi


  # Check main
  MCHECK=`qstat -u ${USER} | grep -is $tag | grep -is main | wc -l`
  if [ $MCHECK -gt 0 ]; then
      MRES="RUNN"
  else
      lastfile=`ls ${folder}/main_${tag}*out | awk '{ f=$NF };END{ print f }'`
      if [ -f $lastfile ]; then
          MCHECK=`grep -is 'ED execution ends' $lastfile | wc -l`
          if [ $MCHECK -gt 0 ];then
              MRES="COMP"
          else
              MRES="FAIL"
          fi
      else
          MRES="----"
      fi
  fi

  echo $tag"  "$DRES $TRES $MRES

done < $LIST

# clean-up

rm -f testfiles.txt
rm -f mainfiles.txt
rm -f dbugfiles.txt