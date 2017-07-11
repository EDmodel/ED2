#!/bin/bash

# This script reads a list of files from a txt
# Then it moves all R files that are not in the list to a different folder
# contains function usage:
# contains aList anItem

#readarray -t KEEP < useful.txt
#set -x
while IFS=\= read var; do
    KEEP+=($var)
done < useful.txt

echo "${KEEP[@]}"

#contains() {
 #   [[ $1 =~ (^| )$2($| ) ]] && exit(0) || exit(1)}

for item in *.{r,R}
do
	if [[ ${KEEP[*]} =~ $item ]]
		then
		echo -e "skipping $item (to keep)\n"
	else
		echo -e "Moving item $item to not useful folder\n"
		mv $item maybe-not-useful/
	fi
done

echo '\n Done! \n'
exit 0