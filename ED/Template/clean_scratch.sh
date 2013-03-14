#!/bin/sh
myself=`whoami`

#----- Check which queue we will clean. ---------------------------------------------------#
if [ 'x'${1} == 'x' ]
then
   echo -n 'Which queue do you want to clean?'
   read queue
else
   queue=${1}
fi
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Select the nodes to be deleted.                                                     #
#------------------------------------------------------------------------------------------#
case ${queue} in
moorcroft_6100a)
   nodes="moorcroft01 moorcroft02 moorcroft03 moorcroft04"
   ;;
moorcroft_6100b)
   nodes="moorcroft05 moorcroft06 moorcroft07 moorcroft08 moorcroft09
          moorcroft10 moorcroft11 moorcroft12 moorcroft13 moorcroft14
          moorcroft15 moorcroft16 moorcroft17 moorcroft18 moorcroft19
          moorcroft20 moorcroft21 moorcroft22 moorcroft23 moorcroft24
          moorcroft25 moorcroft26 moorcroft27 moorcroft28 moorcroft29
          moorcroft30 moorcroft31 moorcroft32 moorcroft33 moorcroft34
          moorcroft35 moorcroft36 moorcroft37 moorcroft38 moorcroft39
          moorcroft40 moorcroft41 moorcroft42 moorcroft43 moorcroft44"
   ;;
moorcroft2a)
   nodes="hero4001 hero4002 hero4003 hero4004 hero4005 hero4006
          hero4007 hero4008 hero4009 hero4010 hero4011"
   ;;
moorcroft2b)
   nodes="hero4014 hero4015 hero4016 hero4101 hero4102 hero4103
          hero4104 hero4105 hero4106 hero4107 hero4108"
   ;;
moorcroft2c)
   nodes="hero4109 hero4110 hero4111 hero4112 hero4113 hero4114 hero4115 hero4116"
   ;;
wofsy)
   nodes="wofsy011 wofsy012 wofsy013 wofsy014 wofsy021 wofsy022 wofsy023 wofsy024"
   ;;
camd)
   nodes="camd04 camd05 camd06 camd07 camd09 camd10 camd11 camd13"
   ;;
unrestricted_parallel)
   nodes="hero3102 hero3103 hero3104 hero3107 hero3108 hero3109 hero3110
          hero3111 hero3112 hero3113 hero3114 hero3115 hero3116 hero3201
          hero3202 hero3203 hero3204 hero3205 hero3206 hero3207 hero3208
          hero3209 hero3210 hero3211 hero3212 hero3213 hero3214 hero3215
          hero3216"
   ;;
unrestricted_serial)
   nodes="soph57 soph58 soph59 soph60 soph61 soph62 soph63 soph64"
   ;;
long_serial)
   nodes="soph03 soph04 soph05 soph13 soph14 soph16 soph17 soph18 soph19 soph24 
          soph25 soph26 soph27 soph28 soph29 soph30 soph31 soph32 soph33 soph34 
          soph35 soph36 soph37 soph38 soph39 soph40 soph41 soph42 soph43 soph44 
          soph45 soph46 soph47 soph48 soph49 soph50 soph51 soph52 soph53 soph54 
          soph55 soph56"
   ;;
*)
   echo ' I cannot recognise queue '${queue}'...'
   exit 39
   ;;
esac
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#     Delete the files in all nodes that the queue uses.                                   #
#------------------------------------------------------------------------------------------#
for node in ${nodes}
do
   echo -n ' Scheduling files for deletion - node '${node}'...'
   ssh ${node} /bin/mv /scratch/${myself} /scratch/goodbye-${myself} 1> /dev/null 2>&1
   ssh ${node} rm -fr /scratch/goodbye-${myself}                     1> /dev/null 2>&1 &
   echo 'Done!'
done
#------------------------------------------------------------------------------------------#
