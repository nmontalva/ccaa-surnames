#!/bin/bash

for i in `seq 40 237`; do
  wget -O "$i".html "http://www.comunidadesagricolas.cl/index.php?option=com_content&task=view&id=$i"
  pdf=`grep 'Pinche Aqui' "$i".html |grep -o '[a-zA-Z0-9%/_]*\.pdf'`
  wget -O "$i".pdf "http://www.comunidadesagricolas.cl/$pdf" || echo "$i" >> error.txt
  rm -f "$i".html
  pdftohtml "$i".pdf "$i".html
  sed -ie 's/<br>//g' "$i"s.html
  sed -ie 's/&nbsp;/ /g' "$i"s.html
  sed -ie '/^<.*$/d' "$i"s.html
  mv "$i"s.html "$i".txt
  name=`head -n 1 "$i".txt |tr [A-Z] [a-z]`
  # the next line would be clojure -cp . --eval ... in ubuntu
  java -cp /usr/share/clojure/clojure.jar:. clojure.main --eval \
    "(do (load \"bienes-raices\") (convert \"$i.txt\" \"$i - $name.csv\"))" || echo "$i" >> error.txt
  # the backslash at the end of the line means that it continues on the next line
  rm -f "$i".html "$i"s.htmle "$i"_ind.html "$i"-1_1.jpg
done
