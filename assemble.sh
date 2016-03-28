megahit --presets meta -1 $1 -2 $2 -o out
cp out/*.fa .
rm out -R
