#!/bin/bash
here=`pwd`
email='manfredo.diporciaebrugnera@ugent.be'
variables="air.sig995 pres.sfc rhum.sig995 uwnd.sig995 vwnd.sig995 
           dlwrf.sfc.gauss nbdsf.sfc.gauss nddsf.sfc.gauss vbdsf.sfc.gauss
           vddsf.sfc.gauss prate.sfc.gauss pres.sfc.gauss"

years="1948 1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959
                 1960 1961 1962 1963 1964 1965 1966 1967 1968 1969
                 1970 1971 1972 1973 1974 1975 1976 1977 1978 1979
                 1980 1981 1982 1983 1984 1985 1986 1987 1988 1989
                 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999
                 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009
                 2010 2011 2012 2013 2014"

for var in ${variables}
do
   case ${var} in
     *.gauss ) dir="surface_gauss";;
     *       ) dir="surface";;     
   esac

   for year in ${years}
   do
     destination=${here}/${var}
     file=${var}'.'${year}'.nc'
     if [ ! -s ${destination} ]
     then
       echo 'Creating folder '${destination}'...'
       mkdir ${destination}
       cd ${destination}
     else
       cd ${destination}
     fi
     if [ ! -s ${destination}/${file} ]
     then

        echo -n '   ===> Downloading '${file}'...'
cat << Eof | ftp -in ftp.cdc.noaa.gov > /dev/null
user anonymous ${email}
bin
cd Datasets/ncep.reanalysis/${dir}
get ${file}
bye
Eof

       echo 'Done!'
     fi
   done
done

cd ${pwd}
