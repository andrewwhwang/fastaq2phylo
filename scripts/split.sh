while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line#>}.fa
        outfile=${outfile// /.}
        echo $line > $2/output/$outfile
    else
        echo $line >> $2/output/$outfile
    fi
done < $1
