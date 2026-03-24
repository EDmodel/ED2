#!/bin/bash

#===============================================================================
#  This script looks through the ED2IN files and generates
#  a matrix of parameters across sites and across simulation.
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

declare -a SITEGRP1=(M34 S67 HAR PDG TON CAX)
declare -a SITEGRP2=(GYF ATA PET HIM HIP RJG)
declare -a SITETAG1=(m34 s67 har pdg ton cax)
declare -a SITETAG2=(gyf ata pet him hip rjg)


declare -a VARGRP1=(IED_INIT_MODE INTEGRATION_SCHEME \
                      POI_LAT POI_LON ISOILFLG ISOILCOL \
                      NZG NZS ISOILBC IBIGLEAF IBRANCH_THERMO IALLOM IGRASS \
                      IPHEN_SCHEME CROWN_MOD H2O_PLANT_LIM)

declare -a VARGRP2=(IDDMORT_SCHEME THETACRIT INCLUDE_FIRE IANTH_DISTURB \
                   IPERCOL MAXSITE MAXPATCH MAXCOHORT TREEFALL_DISTURBANCE_RATE \
                   IMETAVG IMETRAD IDETAILED PATCH_KEEP)

declare -a VARGRP3=(IPHYSIOL H2O_PLANT_LIM ICANTURB ISFCLYR REPRO_SCHEME \
                    DECOMP_SCHEME ISTRUCT_GROWTH_SCHEME)


declare -a VNAME1=(INIT_MODE INTEGRATION LAT LON ISOILFLG ISOILCOL \
                      NZG NZS ISOILBC IBIGLEAF IBRANCH IALLOM IGRASS \
                      IPHEN CROWNMOD )

declare -a VNAME2=(IDDMORT_SCHEME THETACRIT INCLUDE_FIRE IANTH_DISTURB \
                   IPERCOL MAXSITE MAXPATCH MAXCOHORT TREEFALL \
                   IMETAVG IMETRAD IDETAILED PATCH_KEEP)

declare -a VNAME3=(IPHYSIOL H2O_LIM ICANTURB ISFCLYR IREPRO \
                   IDECOMP STRUCT_GROW)

#cat > $LATEXTABLE <<\EOF
#\documentclass{article}
#\usepackage{fullpage}
#\title{Runtime Parameter Table}
#\usepackage{graphicx}
#\usepackage[margin=1in]{geometry}% http://ctan.org/pkg/geometry
#\usepackage{color}
#\usepackage{amsmath,amsthm,amssymb,amsfonts}
#\usepackage{rotating}
#\usepackage[table]{xcolor}
#\usepackage[normalem]{ulem}
#\renewcommand{\arraystretch}{1.25}% Increase the row height of tabular/array
#\newcolumntype{C}{>{\centering\arraybackslash$}p{\linewidth}<{$}}
#\newcommand{\tagarray}{%
#\mbox{}\refstepcounter{equation}%
#$(\theequation)$%
#}
#%\setlength{\mathindent}{0pt}
#%\author{Ryan G. Knox}
#\begin{document}
#\maketitle
#\tableofcontents



# PAGE 1, GROUP 1 VARS1
# =======================================================
LATEXTABLE=${folder}/report/vartable1.tex
cat > $LATEXTABLE <<\EOF
\frame{
\tiny
\begin{columns}
\begin{column}{\textwidth}
\begin{center}
EOF

let NVARS=${#VARGRP1[@]}
let NVARSM1=NVARS-1
let NVARSP1=NVARS+1
let NSITE=${#SITEGRP1[@]}
let NSITEM1=NSITE-1

unset VARSTR
unset VNAME
unset SITEPFX
unset TVAL
unset MVAL
unset DVAL

for i in `seq 0 $NVARSM1`
do
    VARSTR[i]=${VARGRP1[i]}
    VNAME[i]=${VNAME1[i]}
done
for i in `seq 0 $NSITEM1`
do
    SITEPFX[i]=${SITEGRP1[i]}
done

TABHEAD='\begin{tabular}{| r c'
TABVAR="${VNAME[0]//_/\_}"
RSTR='\rowcolor{gray!35} & \rotatebox[origin=c]{90}{\texttt{'$TABVAR'}}'
for i in `seq 1 $NVARSM1`
do
    # CONVERT THE VARIABLE NAME SO LATEX DOESNT FREAK
    TABVAR="${VNAME[i]//_/\_}"
    RSTR=$RSTR' &  \rotatebox[origin=c]{90}{\texttt{'$TABVAR'}}'
    TABHEAD=$TABHEAD' c'
done
RSTR=$RSTR' \\'
TABHEAD=$TABHEAD' |}'
echo $TABHEAD >> $LATEXTABLE
echo '\hline' >> $LATEXTABLE
echo $RSTR >> $LATEXTABLE
for i in ${!SITEPFX[@]}
do
    echo "PROCESSING SITE: "${SITEPFX[i]}
    ED2INFT=${folder}/ED2IN-${SITEPFX[i]}-TEST
    ED2INFD=${folder}/ED2IN-${SITEPFX[i]}-DBUG
    ED2INFM=${folder}/ED2IN-${SITEPFX[i]}-MAIN
    let fcnt=0
    if [ -f $ED2INFT ]; then
	let fcnt=fcnt+1
    fi
    if [ -f $ED2INFD ]; then
	let fcnt=fcnt+1
    fi
    if [ -f $ED2INFM ]; then
	let fcnt=fcnt+1
    fi
    if [ fcnt==3 ]; then
	TABSITE='\multicolumn{'$NVARSP1'}{|>{\columncolor[gray]{0.9}}l|}{'${SITEPFX[i]}'}\\'
	echo $TABSITE >> $LATEXTABLE
	DROW='DBUG: '
	TROW='TEST: '
	MROW='MAIN: '
	# LOOP through the variables
	for j in ${!VARSTR[@]}
	do
	    # FOR THE TEST FILE
	    grep ${VARSTR[j]} $ED2INFT > tmp
	    while iline= read -r LINE
	    do
		SUBSTR=`echo $LINE | awk -v varn="NL%${VARSTR[j]}" '$0 ~ varn'`
		NFIELDS=`echo $SUBSTR | awk '{ print NF }'`
		if [ $NFIELDS == 3 ];then
		    TVAL[j]=`echo $SUBSTR | awk '{ print $3 }'`
		elif [ $NFIELDS == 1 ] || [ $NFIELDS == 2 ]; then
		    echo "" 
		    echo "========================="
		    echo " MAKE A SLIGHT MODIFICATION"
		    echo " TO $SUBSTR IN $ED2INFT "
		    echo " SO THE STRING IS 3 FIELDS"
		    echo " EXITING ..."
		    echo "========================="
		    exit
		fi
	    done < tmp

	    # FOR THE DBUG FILE
	    grep ${VARSTR[j]} $ED2INFD > tmp
	    while iline= read -r LINE
	    do
		SUBSTR=`echo $LINE | awk -v varn="NL%${VARSTR[j]}" '$0 ~ varn'`
		NFIELDS=`echo $SUBSTR | awk '{ print NF }'`
		if [ $NFIELDS == 3 ];then
		    DVAL[j]=`echo $SUBSTR | awk '{ print $3 }'`
		elif [ $NFIELDS == 1 ] || [ $NFIELDS == 2 ]; then
		    echo "" 
		    echo "========================="
		    echo " MAKE A SLIGHT MODIFICATION"
		    echo " TO $SUBSTR IN $ED2INFT "
		    echo " SO THE STRING IS 3 FIELDS"
		    echo " EXITING ..."
		    echo "========================="
		    exit
		fi
	    done < tmp

	    # FOR THE MAIN FILE
	    grep ${VARSTR[j]} $ED2INFM > tmp
	    while iline= read -r LINE
	    do
		SUBSTR=`echo $LINE | awk -v varn="NL%${VARSTR[j]}" '$0 ~ varn'`
		NFIELDS=`echo $SUBSTR | awk '{ print NF }'`
		if [ $NFIELDS == 3 ];then
		    MVAL[j]=`echo $SUBSTR | awk '{ print $3 }'`
		elif [ $NFIELDS == 1 ] || [ $NFIELDS == 2 ]; then
		    echo "" 
		    echo "========================="
		    echo " MAKE A SLIGHT MODIFICATION"
		    echo " TO $SUBSTR IN $ED2INFT "
		    echo " SO THE STRING IS 3 FIELDS"
		    echo " EXITING ..."
		    echo "========================="
		    exit
		fi
	    done < tmp

	    TROW=$TROW" & "${TVAL[j]}
	    DROW=$DROW" & "${DVAL[j]}
	    MROW=$MROW" & "${MVAL[j]}

	done
	
	echo $TROW' \\' >> $LATEXTABLE
	echo $DROW' \\' >> $LATEXTABLE
	echo $MROW' \\' >> $LATEXTABLE

	echo ${TVAL[@]}
	echo ${DVAL[@]}
	echo ${MVAL[@]}
    fi


done

cat >> $LATEXTABLE <<\EOF
  \hline
    \end{tabular}
  \end{center}
\end{column}
\end{columns}
}
EOF




# PAGE 2, GROUP 2 VARS1
# =======================================================
LATEXTABLE=${folder}/report/vartable2.tex
cat > $LATEXTABLE <<\EOF
\frame{
\tiny
\begin{columns}
\begin{column}{\textwidth}
\begin{center}
EOF

let NVARS=${#VARGRP1[@]}
let NVARSM1=NVARS-1
let NVARSP1=NVARS+1
let NSITE=${#SITEGRP2[@]}
let NSITEM1=NSITE-1

unset VARSTR
unset SITEPFX
unset TVAL
unset MVAL
unset DVAL

for i in `seq 0 $NVARSM1`
do
    VARSTR[i]=${VARGRP1[i]}
    VNAME[i]=${VNAME1[i]}
done
for i in `seq 0 $NSITEM1`
do
    SITEPFX[i]=${SITEGRP2[i]}
done

TABHEAD='\begin{tabular}{| r c'
TABVAR="${VNAME[0]//_/\_}"
RSTR='\rowcolor{gray!35} & \rotatebox[origin=c]{90}{\texttt{'$TABVAR'}}'
for i in `seq 1 $NVARSM1`
do
    # CONVERT THE VARIABLE NAME SO LATEX DOESNT FREAK
    TABVAR="${VNAME[i]//_/\_}"
    RSTR=$RSTR' &  \rotatebox[origin=c]{90}{\texttt{'$TABVAR'}}'
    TABHEAD=$TABHEAD' c'
done
RSTR=$RSTR' \\'
TABHEAD=$TABHEAD' |}'
echo $TABHEAD >> $LATEXTABLE
echo '\hline' >> $LATEXTABLE
echo $RSTR >> $LATEXTABLE
for i in ${!SITEPFX[@]}
do
    echo "PROCESSING SITE: "${SITEPFX[i]}
    ED2INFT=${folder}/ED2IN-${SITEPFX[i]}-TEST
    ED2INFD=${folder}/ED2IN-${SITEPFX[i]}-DBUG
    ED2INFM=${folder}/ED2IN-${SITEPFX[i]}-MAIN
    let fcnt=0
    if [ -f $ED2INFT ]; then
	let fcnt=fcnt+1
    fi
    if [ -f $ED2INFD ]; then
	let fcnt=fcnt+1
    fi
    if [ -f $ED2INFM ]; then
	let fcnt=fcnt+1
    fi
    if [ fcnt==3 ]; then
	TABSITE='\multicolumn{'$NVARSP1'}{|>{\columncolor[gray]{0.9}}l|}{'${SITEPFX[i]}'}\\'
	echo $TABSITE >> $LATEXTABLE
	DROW='DBUG: '
	TROW='TEST: '
	MROW='MAIN: '
	# LOOP through the variables
	for j in ${!VARSTR[@]}
	do
	    # FOR THE TEST FILE
	    grep ${VARSTR[j]} $ED2INFT > tmp
	    while iline= read -r LINE
	    do
		SUBSTR=`echo $LINE | awk -v varn="NL%${VARSTR[j]}" '$0 ~ varn'`
		NFIELDS=`echo $SUBSTR | awk '{ print NF }'`
		if [ $NFIELDS == 3 ];then
		    TVAL[j]=`echo $SUBSTR | awk '{ print $3 }'`
		elif [ $NFIELDS == 1 ] || [ $NFIELDS == 2 ]; then
		    echo "" 
		    echo "========================="
		    echo " MAKE A SLIGHT MODIFICATION"
		    echo " TO $SUBSTR IN $ED2INFT "
		    echo " SO THE STRING IS 3 FIELDS"
		    echo " EXITING ..."
		    echo "========================="
		    exit
		fi
	    done < tmp

	    # FOR THE DBUG FILE
	    grep ${VARSTR[j]} $ED2INFD > tmp
	    while iline= read -r LINE
	    do
		SUBSTR=`echo $LINE | awk -v varn="NL%${VARSTR[j]}" '$0 ~ varn'`
		NFIELDS=`echo $SUBSTR | awk '{ print NF }'`
		if [ $NFIELDS == 3 ];then
		    DVAL[j]=`echo $SUBSTR | awk '{ print $3 }'`
		elif [ $NFIELDS == 1 ] || [ $NFIELDS == 2 ]; then
		    echo "" 
		    echo "========================="
		    echo " MAKE A SLIGHT MODIFICATION"
		    echo " TO $SUBSTR IN $ED2INFT "
		    echo " SO THE STRING IS 3 FIELDS"
		    echo " EXITING ..."
		    echo "========================="
		    exit
		fi
	    done < tmp

	    # FOR THE MAIN FILE
	    grep ${VARSTR[j]} $ED2INFM > tmp
	    while iline= read -r LINE
	    do
		SUBSTR=`echo $LINE | awk -v varn="NL%${VARSTR[j]}" '$0 ~ varn'`
		NFIELDS=`echo $SUBSTR | awk '{ print NF }'`
		if [ $NFIELDS == 3 ];then
		    MVAL[j]=`echo $SUBSTR | awk '{ print $3 }'`
		elif [ $NFIELDS == 1 ] || [ $NFIELDS == 2 ]; then
		    echo "" 
		    echo "========================="
		    echo " MAKE A SLIGHT MODIFICATION"
		    echo " TO $SUBSTR IN $ED2INFT "
		    echo " SO THE STRING IS 3 FIELDS"
		    echo " EXITING ..."
		    echo "========================="
		    exit
		fi
	    done < tmp

	    TROW=$TROW" & "${TVAL[j]}
	    DROW=$DROW" & "${DVAL[j]}
	    MROW=$MROW" & "${MVAL[j]}

	done
	
	echo $TROW' \\' >> $LATEXTABLE
	echo $DROW' \\' >> $LATEXTABLE
	echo $MROW' \\' >> $LATEXTABLE

	echo ${TVAL[@]}
	echo ${DVAL[@]}
	echo ${MVAL[@]}
    fi


done

cat >> $LATEXTABLE <<\EOF
  \hline
    \end{tabular}
  \end{center}
\end{column}
\end{columns}
}
EOF


# PAGE 3, GROUP 1 VARS2
# =======================================================
LATEXTABLE=${folder}/report/vartable3.tex
cat > $LATEXTABLE <<\EOF
\frame{
\tiny
\begin{columns}
\begin{column}{\textwidth}
\begin{center}
EOF

let NVARS=${#VARGRP2[@]}
let NVARSM1=NVARS-1
let NVARSP1=NVARS+1
let NSITE=${#SITEGRP1[@]}
let NSITEM1=NSITE-1

unset VARSTR
unset VNAME
unset SITEPFX
unset TVAL
unset MVAL
unset DVAL

for i in `seq 0 $NVARSM1`
do
    VARSTR[i]=${VARGRP2[i]}
    VNAME[i]=${VNAME2[i]}
done
for i in `seq 0 $NSITEM1`
do
    SITEPFX[i]=${SITEGRP1[i]}
done

TABHEAD='\begin{tabular}{| r c'
TABVAR="${VNAME[0]//_/\_}"
RSTR='\rowcolor{gray!35} & \rotatebox[origin=c]{90}{\texttt{'$TABVAR'}}'
for i in `seq 1 $NVARSM1`
do
    # CONVERT THE VARIABLE NAME SO LATEX DOESNT FREAK
    TABVAR="${VNAME[i]//_/\_}"
    RSTR=$RSTR' &  \rotatebox[origin=c]{90}{\texttt{'$TABVAR'}}'
    TABHEAD=$TABHEAD' c'
done
RSTR=$RSTR' \\'
TABHEAD=$TABHEAD' |}'
echo $TABHEAD >> $LATEXTABLE
echo '\hline' >> $LATEXTABLE
echo $RSTR >> $LATEXTABLE
for i in ${!SITEPFX[@]}
do
    echo "PROCESSING SITE: "${SITEPFX[i]}
    ED2INFT=${folder}/ED2IN-${SITEPFX[i]}-TEST
    ED2INFD=${folder}/ED2IN-${SITEPFX[i]}-DBUG
    ED2INFM=${folder}/ED2IN-${SITEPFX[i]}-MAIN
    let fcnt=0
    if [ -f $ED2INFT ]; then
	let fcnt=fcnt+1
    fi
    if [ -f $ED2INFD ]; then
	let fcnt=fcnt+1
    fi
    if [ -f $ED2INFM ]; then
	let fcnt=fcnt+1
    fi
    if [ fcnt==3 ]; then
	TABSITE='\multicolumn{'$NVARSP1'}{|>{\columncolor[gray]{0.9}}l|}{'${SITEPFX[i]}'}\\'
	echo $TABSITE >> $LATEXTABLE
	DROW='DBUG: '
	TROW='TEST: '
	MROW='MAIN: '
	# LOOP through the variables
	for j in ${!VARSTR[@]}
	do
	    # FOR THE TEST FILE
	    grep ${VARSTR[j]} $ED2INFT > tmp
	    while iline= read -r LINE
	    do
		SUBSTR=`echo $LINE | awk -v varn="NL%${VARSTR[j]}" '$0 ~ varn'`
		NFIELDS=`echo $SUBSTR | awk '{ print NF }'`
		if [ $NFIELDS == 3 ];then
		    TVAL[j]=`echo $SUBSTR | awk '{ print $3 }'`
		elif [ $NFIELDS == 1 ] || [ $NFIELDS == 2 ]; then
		    echo "" 
		    echo "========================="
		    echo " MAKE A SLIGHT MODIFICATION"
		    echo " TO $SUBSTR IN $ED2INFT "
		    echo " SO THE STRING IS 3 FIELDS"
		    echo " EXITING ..."
		    echo "========================="
		    exit
		fi
	    done < tmp

	    # FOR THE DBUG FILE
	    grep ${VARSTR[j]} $ED2INFD > tmp
	    while iline= read -r LINE
	    do
		SUBSTR=`echo $LINE | awk -v varn="NL%${VARSTR[j]}" '$0 ~ varn'`
		NFIELDS=`echo $SUBSTR | awk '{ print NF }'`
		if [ $NFIELDS == 3 ];then
		    DVAL[j]=`echo $SUBSTR | awk '{ print $3 }'`
		elif [ $NFIELDS == 1 ] || [ $NFIELDS == 2 ]; then
		    echo "" 
		    echo "========================="
		    echo " MAKE A SLIGHT MODIFICATION"
		    echo " TO $SUBSTR IN $ED2INFT "
		    echo " SO THE STRING IS 3 FIELDS"
		    echo " EXITING ..."
		    echo "========================="
		    exit
		fi
	    done < tmp

	    # FOR THE MAIN FILE
	    grep ${VARSTR[j]} $ED2INFM > tmp
	    while iline= read -r LINE
	    do
		SUBSTR=`echo $LINE | awk -v varn="NL%${VARSTR[j]}" '$0 ~ varn'`
		NFIELDS=`echo $SUBSTR | awk '{ print NF }'`
		if [ $NFIELDS == 3 ];then
		    MVAL[j]=`echo $SUBSTR | awk '{ print $3 }'`
		elif [ $NFIELDS == 1 ] || [ $NFIELDS == 2 ]; then
		    echo "" 
		    echo "========================="
		    echo " MAKE A SLIGHT MODIFICATION"
		    echo " TO $SUBSTR IN $ED2INFT "
		    echo " SO THE STRING IS 3 FIELDS"
		    echo " EXITING ..."
		    echo "========================="
		    exit
		fi
	    done < tmp

	    TROW=$TROW" & "${TVAL[j]}
	    DROW=$DROW" & "${DVAL[j]}
	    MROW=$MROW" & "${MVAL[j]}

	done
	
	echo $TROW' \\' >> $LATEXTABLE
	echo $DROW' \\' >> $LATEXTABLE
	echo $MROW' \\' >> $LATEXTABLE

	echo ${TVAL[@]}
	echo ${DVAL[@]}
	echo ${MVAL[@]}
    fi


done

cat >> $LATEXTABLE <<\EOF
  \hline
    \end{tabular}
  \end{center}
\end{column}
\end{columns}
}
EOF


# PAGE 4, GROUP 2 VARS2
# =======================================================
LATEXTABLE=${folder}/report/vartable4.tex
cat > $LATEXTABLE <<\EOF
\frame{
\tiny
\begin{columns}
\begin{column}{\textwidth}
\begin{center}
EOF

let NVARS=${#VARGRP2[@]}
let NVARSM1=NVARS-1
let NVARSP1=NVARS+1
let NSITE=${#SITEGRP2[@]}
let NSITEM1=NSITE-1

unset VARSTR
unset VNAME
unset SITEPFX
unset TVAL
unset MVAL
unset DVAL

for i in `seq 0 $NVARSM1`
do
    VARSTR[i]=${VARGRP2[i]}
    VNAME[i]=${VNAME2[i]}
done
for i in `seq 0 $NSITEM1`
do
    SITEPFX[i]=${SITEGRP2[i]}
done

TABHEAD='\begin{tabular}{| r c'
TABVAR="${VNAME[0]//_/\_}"
RSTR='\rowcolor{gray!35} & \rotatebox[origin=c]{90}{\texttt{'$TABVAR'}}'
for i in `seq 1 $NVARSM1`
do
    # CONVERT THE VARIABLE NAME SO LATEX DOESNT FREAK
    TABVAR="${VNAME[i]//_/\_}"
    RSTR=$RSTR' &  \rotatebox[origin=c]{90}{\texttt{'$TABVAR'}}'
    TABHEAD=$TABHEAD' c'
done
RSTR=$RSTR' \\'
TABHEAD=$TABHEAD' |}'
echo $TABHEAD >> $LATEXTABLE
echo '\hline' >> $LATEXTABLE
echo $RSTR >> $LATEXTABLE
for i in ${!SITEPFX[@]}
do
    echo "PROCESSING SITE: "${SITEPFX[i]}
    ED2INFT=${folder}/ED2IN-${SITEPFX[i]}-TEST
    ED2INFD=${folder}/ED2IN-${SITEPFX[i]}-DBUG
    ED2INFM=${folder}/ED2IN-${SITEPFX[i]}-MAIN
    let fcnt=0
    if [ -f $ED2INFT ]; then
	let fcnt=fcnt+1
    fi
    if [ -f $ED2INFD ]; then
	let fcnt=fcnt+1
    fi
    if [ -f $ED2INFM ]; then
	let fcnt=fcnt+1
    fi
    if [ fcnt==3 ]; then
	TABSITE='\multicolumn{'$NVARSP1'}{|>{\columncolor[gray]{0.9}}l|}{'${SITEPFX[i]}'}\\'
	echo $TABSITE >> $LATEXTABLE
	DROW='DBUG: '
	TROW='TEST: '
	MROW='MAIN: '
	# LOOP through the variables
	for j in ${!VARSTR[@]}
	do
	    # FOR THE TEST FILE
	    grep ${VARSTR[j]} $ED2INFT > tmp
	    while iline= read -r LINE
	    do
		SUBSTR=`echo $LINE | awk -v varn="NL%${VARSTR[j]}" '$0 ~ varn'`
		NFIELDS=`echo $SUBSTR | awk '{ print NF }'`
		if [ $NFIELDS == 3 ];then
		    TVAL[j]=`echo $SUBSTR | awk '{ print $3 }'`
		elif [ $NFIELDS == 1 ] || [ $NFIELDS == 2 ]; then
		    echo "" 
		    echo "========================="
		    echo " MAKE A SLIGHT MODIFICATION"
		    echo " TO $SUBSTR IN $ED2INFT "
		    echo " SO THE STRING IS 3 FIELDS"
		    echo " EXITING ..."
		    echo "========================="
		    exit
		fi
	    done < tmp

	    # FOR THE DBUG FILE
	    grep ${VARSTR[j]} $ED2INFD > tmp
	    while iline= read -r LINE
	    do
		SUBSTR=`echo $LINE | awk -v varn="NL%${VARSTR[j]}" '$0 ~ varn'`
		NFIELDS=`echo $SUBSTR | awk '{ print NF }'`
		if [ $NFIELDS == 3 ];then
		    DVAL[j]=`echo $SUBSTR | awk '{ print $3 }'`
		elif [ $NFIELDS == 1 ] || [ $NFIELDS == 2 ]; then
		    echo "" 
		    echo "========================="
		    echo " MAKE A SLIGHT MODIFICATION"
		    echo " TO $SUBSTR IN $ED2INFT "
		    echo " SO THE STRING IS 3 FIELDS"
		    echo " EXITING ..."
		    echo "========================="
		    exit
		fi
	    done < tmp

	    # FOR THE MAIN FILE
	    grep ${VARSTR[j]} $ED2INFM > tmp
	    while iline= read -r LINE
	    do
		SUBSTR=`echo $LINE | awk -v varn="NL%${VARSTR[j]}" '$0 ~ varn'`
		NFIELDS=`echo $SUBSTR | awk '{ print NF }'`
		if [ $NFIELDS == 3 ];then
		    MVAL[j]=`echo $SUBSTR | awk '{ print $3 }'`
		elif [ $NFIELDS == 1 ] || [ $NFIELDS == 2 ]; then
		    echo "" 
		    echo "========================="
		    echo " MAKE A SLIGHT MODIFICATION"
		    echo " TO $SUBSTR IN $ED2INFT "
		    echo " SO THE STRING IS 3 FIELDS"
		    echo " EXITING ..."
		    echo "========================="
		    exit
		fi
	    done < tmp

	    TROW=$TROW" & "${TVAL[j]}
	    DROW=$DROW" & "${DVAL[j]}
	    MROW=$MROW" & "${MVAL[j]}

	done
	
	echo $TROW' \\' >> $LATEXTABLE
	echo $DROW' \\' >> $LATEXTABLE
	echo $MROW' \\' >> $LATEXTABLE

	echo ${TVAL[@]}
	echo ${DVAL[@]}
	echo ${MVAL[@]}
    fi


done

cat >> $LATEXTABLE <<\EOF
  \hline
    \end{tabular}
  \end{center}
\end{column}
\end{columns}
}
EOF



















# PAGE 5, GROUP 1 VARS3
# =======================================================
LATEXTABLE=${folder}/report/vartable5.tex
cat > $LATEXTABLE <<\EOF
\frame{
\tiny
\begin{columns}
\begin{column}{\textwidth}
\begin{center}
EOF

let NVARS=${#VARGRP3[@]}
let NVARSM1=NVARS-1
let NVARSP1=NVARS+1
let NSITE=${#SITEGRP1[@]}
let NSITEM1=NSITE-1

unset VARSTR
unset VNAME
unset SITEPFX
unset TVAL
unset MVAL
unset DVAL

for i in `seq 0 $NVARSM1`
do
    VARSTR[i]=${VARGRP3[i]}
    VNAME[i]=${VNAME3[i]}
done
for i in `seq 0 $NSITEM1`
do
    SITEPFX[i]=${SITEGRP1[i]}
done

TABHEAD='\begin{tabular}{| r c'
TABVAR="${VNAME[0]//_/\_}"
RSTR='\rowcolor{gray!35} & \rotatebox[origin=c]{90}{\texttt{'$TABVAR'}}'
for i in `seq 1 $NVARSM1`
do
    # CONVERT THE VARIABLE NAME SO LATEX DOESNT FREAK
    TABVAR="${VNAME[i]//_/\_}"
    RSTR=$RSTR' &  \rotatebox[origin=c]{90}{\texttt{'$TABVAR'}}'
    TABHEAD=$TABHEAD' c'
done
RSTR=$RSTR' \\'
TABHEAD=$TABHEAD' |}'
echo $TABHEAD >> $LATEXTABLE
echo '\hline' >> $LATEXTABLE
echo $RSTR >> $LATEXTABLE
for i in ${!SITEPFX[@]}
do
    echo "PROCESSING SITE: "${SITEPFX[i]}
    ED2INFT=${folder}/ED2IN-${SITEPFX[i]}-TEST
    ED2INFD=${folder}/ED2IN-${SITEPFX[i]}-DBUG
    ED2INFM=${folder}/ED2IN-${SITEPFX[i]}-MAIN
    let fcnt=0
    if [ -f $ED2INFT ]; then
	let fcnt=fcnt+1
    fi
    if [ -f $ED2INFD ]; then
	let fcnt=fcnt+1
    fi
    if [ -f $ED2INFM ]; then
	let fcnt=fcnt+1
    fi
    if [ fcnt==3 ]; then
	TABSITE='\multicolumn{'$NVARSP1'}{|>{\columncolor[gray]{0.9}}l|}{'${SITEPFX[i]}'}\\'
	echo $TABSITE >> $LATEXTABLE
	DROW='DBUG: '
	TROW='TEST: '
	MROW='MAIN: '
	# LOOP through the variables
	for j in ${!VARSTR[@]}
	do
	    # FOR THE TEST FILE
	    grep ${VARSTR[j]} $ED2INFT > tmp
	    while iline= read -r LINE
	    do
		SUBSTR=`echo $LINE | awk -v varn="NL%${VARSTR[j]}" '$0 ~ varn'`
		NFIELDS=`echo $SUBSTR | awk '{ print NF }'`
		if [ $NFIELDS == 3 ];then
		    TVAL[j]=`echo $SUBSTR | awk '{ print $3 }'`
		elif [ $NFIELDS == 1 ] || [ $NFIELDS == 2 ]; then
		    echo "" 
		    echo "========================="
		    echo " MAKE A SLIGHT MODIFICATION"
		    echo " TO $SUBSTR IN $ED2INFT "
		    echo " SO THE STRING IS 3 FIELDS"
		    echo " EXITING ..."
		    echo "========================="
		    exit
		fi
	    done < tmp

	    # FOR THE DBUG FILE
	    grep ${VARSTR[j]} $ED2INFD > tmp
	    while iline= read -r LINE
	    do
		SUBSTR=`echo $LINE | awk -v varn="NL%${VARSTR[j]}" '$0 ~ varn'`
		NFIELDS=`echo $SUBSTR | awk '{ print NF }'`
		if [ $NFIELDS == 3 ];then
		    DVAL[j]=`echo $SUBSTR | awk '{ print $3 }'`
		elif [ $NFIELDS == 1 ] || [ $NFIELDS == 2 ]; then
		    echo "" 
		    echo "========================="
		    echo " MAKE A SLIGHT MODIFICATION"
		    echo " TO $SUBSTR IN $ED2INFT "
		    echo " SO THE STRING IS 3 FIELDS"
		    echo " EXITING ..."
		    echo "========================="
		    exit
		fi
	    done < tmp

	    # FOR THE MAIN FILE
	    grep ${VARSTR[j]} $ED2INFM > tmp
	    while iline= read -r LINE
	    do
		SUBSTR=`echo $LINE | awk -v varn="NL%${VARSTR[j]}" '$0 ~ varn'`
		NFIELDS=`echo $SUBSTR | awk '{ print NF }'`
		if [ $NFIELDS == 3 ];then
		    MVAL[j]=`echo $SUBSTR | awk '{ print $3 }'`
		elif [ $NFIELDS == 1 ] || [ $NFIELDS == 2 ]; then
		    echo "" 
		    echo "========================="
		    echo " MAKE A SLIGHT MODIFICATION"
		    echo " TO $SUBSTR IN $ED2INFT "
		    echo " SO THE STRING IS 3 FIELDS"
		    echo " EXITING ..."
		    echo "========================="
		    exit
		fi
	    done < tmp

	    TROW=$TROW" & "${TVAL[j]}
	    DROW=$DROW" & "${DVAL[j]}
	    MROW=$MROW" & "${MVAL[j]}

	done
	
	echo $TROW' \\' >> $LATEXTABLE
	echo $DROW' \\' >> $LATEXTABLE
	echo $MROW' \\' >> $LATEXTABLE

	echo ${TVAL[@]}
	echo ${DVAL[@]}
	echo ${MVAL[@]}
    fi


done

cat >> $LATEXTABLE <<\EOF
  \hline
    \end{tabular}
  \end{center}
\end{column}
\end{columns}
}
EOF


# PAGE 6, GROUP 2 VARS3
# =======================================================
LATEXTABLE=${folder}/report/vartable6.tex
cat > $LATEXTABLE <<\EOF
\frame{
\tiny
\begin{columns}
\begin{column}{\textwidth}
\begin{center}
EOF

let NVARS=${#VARGRP3[@]}
let NVARSM1=NVARS-1
let NVARSP1=NVARS+1
let NSITE=${#SITEGRP2[@]}
let NSITEM1=NSITE-1

unset VARSTR
unset VNAME
unset SITEPFX
unset TVAL
unset MVAL
unset DVAL

for i in `seq 0 $NVARSM1`
do
    VARSTR[i]=${VARGRP3[i]}
    VNAME[i]=${VNAME3[i]}
done
for i in `seq 0 $NSITEM1`
do
    SITEPFX[i]=${SITEGRP2[i]}
done

TABHEAD='\begin{tabular}{| r c'
TABVAR="${VNAME[0]//_/\_}"
RSTR='\rowcolor{gray!35} & \rotatebox[origin=c]{90}{\texttt{'$TABVAR'}}'
for i in `seq 1 $NVARSM1`
do
    # CONVERT THE VARIABLE NAME SO LATEX DOESNT FREAK
    TABVAR="${VNAME[i]//_/\_}"
    RSTR=$RSTR' &  \rotatebox[origin=c]{90}{\texttt{'$TABVAR'}}'
    TABHEAD=$TABHEAD' c'
done
RSTR=$RSTR' \\'
TABHEAD=$TABHEAD' |}'
echo $TABHEAD >> $LATEXTABLE
echo '\hline' >> $LATEXTABLE
echo $RSTR >> $LATEXTABLE
for i in ${!SITEPFX[@]}
do
    echo "PROCESSING SITE: "${SITEPFX[i]}
    ED2INFT=${folder}/ED2IN-${SITEPFX[i]}-TEST
    ED2INFD=${folder}/ED2IN-${SITEPFX[i]}-DBUG
    ED2INFM=${folder}/ED2IN-${SITEPFX[i]}-MAIN
    let fcnt=0
    if [ -f $ED2INFT ]; then
	let fcnt=fcnt+1
    fi
    if [ -f $ED2INFD ]; then
	let fcnt=fcnt+1
    fi
    if [ -f $ED2INFM ]; then
	let fcnt=fcnt+1
    fi
    if [ fcnt==3 ]; then
	TABSITE='\multicolumn{'$NVARSP1'}{|>{\columncolor[gray]{0.9}}l|}{'${SITEPFX[i]}'}\\'
	echo $TABSITE >> $LATEXTABLE
	DROW='DBUG: '
	TROW='TEST: '
	MROW='MAIN: '
	# LOOP through the variables
	for j in ${!VARSTR[@]}
	do
	    # FOR THE TEST FILE
	    grep ${VARSTR[j]} $ED2INFT > tmp
	    while iline= read -r LINE
	    do
		SUBSTR=`echo $LINE | awk -v varn="NL%${VARSTR[j]}" '$0 ~ varn'`
		NFIELDS=`echo $SUBSTR | awk '{ print NF }'`
		if [ $NFIELDS == 3 ];then
		    TVAL[j]=`echo $SUBSTR | awk '{ print $3 }'`
		elif [ $NFIELDS == 1 ] || [ $NFIELDS == 2 ]; then
		    echo "" 
		    echo "========================="
		    echo " MAKE A SLIGHT MODIFICATION"
		    echo " TO $SUBSTR IN $ED2INFT "
		    echo " SO THE STRING IS 3 FIELDS"
		    echo " EXITING ..."
		    echo "========================="
		    exit
		fi
	    done < tmp

	    # FOR THE DBUG FILE
	    grep ${VARSTR[j]} $ED2INFD > tmp
	    while iline= read -r LINE
	    do
		SUBSTR=`echo $LINE | awk -v varn="NL%${VARSTR[j]}" '$0 ~ varn'`
		NFIELDS=`echo $SUBSTR | awk '{ print NF }'`
		if [ $NFIELDS == 3 ];then
		    DVAL[j]=`echo $SUBSTR | awk '{ print $3 }'`
		elif [ $NFIELDS == 1 ] || [ $NFIELDS == 2 ]; then
		    echo "" 
		    echo "========================="
		    echo " MAKE A SLIGHT MODIFICATION"
		    echo " TO $SUBSTR IN $ED2INFT "
		    echo " SO THE STRING IS 3 FIELDS"
		    echo " EXITING ..."
		    echo "========================="
		    exit
		fi
	    done < tmp

	    # FOR THE MAIN FILE
	    grep ${VARSTR[j]} $ED2INFM > tmp
	    while iline= read -r LINE
	    do
		SUBSTR=`echo $LINE | awk -v varn="NL%${VARSTR[j]}" '$0 ~ varn'`
		NFIELDS=`echo $SUBSTR | awk '{ print NF }'`
		if [ $NFIELDS == 3 ];then
		    MVAL[j]=`echo $SUBSTR | awk '{ print $3 }'`
		elif [ $NFIELDS == 1 ] || [ $NFIELDS == 2 ]; then
		    echo "" 
		    echo "========================="
		    echo " MAKE A SLIGHT MODIFICATION"
		    echo " TO $SUBSTR IN $ED2INFT "
		    echo " SO THE STRING IS 3 FIELDS"
		    echo " EXITING ..."
		    echo "========================="
		    exit
		fi
	    done < tmp

	    TROW=$TROW" & "${TVAL[j]}
	    DROW=$DROW" & "${DVAL[j]}
	    MROW=$MROW" & "${MVAL[j]}

	done
	
	echo $TROW' \\' >> $LATEXTABLE
	echo $DROW' \\' >> $LATEXTABLE
	echo $MROW' \\' >> $LATEXTABLE

	echo ${TVAL[@]}
	echo ${DVAL[@]}
	echo ${MVAL[@]}
    fi


done

cat >> $LATEXTABLE <<\EOF
  \hline
    \end{tabular}
  \end{center}
\end{column}
\end{columns}
}
EOF


