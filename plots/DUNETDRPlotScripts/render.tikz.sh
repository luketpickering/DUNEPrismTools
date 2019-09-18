INP=${1}
OUP=${2}

if [ -z ${INP} ] || [ ! -e ${INP} ]; then
  echo "Input directory: ${INP} does not exist."
  exit 1
fi

if [ -z ${OUP} ] || [ ! -e ${OUP} ]; then
  echo "Output directory: ${OUP} does not exist."
  exit 1
fi

for i in ${1}/*.tex; do

  IFNAME=${i##*/}
  echo "\documentclass[crop,tikz]{standalone}" > ${OUP}/${IFNAME}
  echo "\usepackage{tikz}" >> ${OUP}/${IFNAME}
  echo "\begin{document}" >> ${OUP}/${IFNAME}
  cat ${i} >> ${OUP}/${IFNAME}
  echo "\end{document}" >> ${OUP}/${IFNAME}
  cd ${OUP};
  pdflatex ${IFNAME}
  cd ../

done
