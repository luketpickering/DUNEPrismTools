#!/bin/sh

echo "" > gaus_manifest.txt

if [[ ! -e "${1}" ]]; then
  echo "No input file found."
  exit 1
fi

OUTFILE=D00NPrism_gaus.root
if [[ "${2}" ]]; then
  OUTFILE=${2}
fi

GLIMIT=4
if [[ "${3}" ]]; then
  GLIMIT=${3}
fi

for i in {1..7}; do
  if [ $i -eq 1 ]; then
    OUTCMD="o"
  else
    OUTCMD="a"
  fi

  PEAK=$(bc <<< "scale=4; ${i} * 0.5")
  WDTH=$(bc <<< "scale=4; ${PEAK} * 0.1")
  DNAME=$(echo "${PEAK}__${WDTH}" | sed "s/\./_/g")
  echo "./FitFluxes.exe -g ${PEAK},${WDTH} -gl ${GLIMIT} -${OUTCMD} ${OUTFILE} -d ${DNAME} -f ${1}, -n 100000 -c 30"
  ./FitFluxes.exe -g ${PEAK},${WDTH} -gl ${GLIMIT} -${OUTCMD} ${OUTFILE} -d ${DNAME} -f ${1} -n 100000 -c 30

  echo "${DNAME}" >> gaus_manifest.txt
done
