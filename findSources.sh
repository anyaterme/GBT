#!/bin/sh
echo "File 01"
list_srcs="3C286 DR21-OH G10_47 G10_62 G29_96 G31_41 G34_26 G34_43MM1 G34_43MM4 G45_07 G45_12 I20126"
for src in $list_srcs; do
   find . -name 'sb.**' -exec grep -H $src {} \;
done

echo "File 03"
list_srcs="3C286 G29_96 G30_97MM1 G31_41 G31_97MM1"
for src in $list_srcs; do
   find . -name 'sb.**' -exec grep -H $src {} \;
done

echo "File 04"
list_srcs="3C48 G34_43MM4"
for src in $list_srcs; do
   find . -name 'sb.**' -exec grep -H $src {} \;
done
