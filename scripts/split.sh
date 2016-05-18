while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line#>}.fa
        outfile=${outfile// /.}
        echo $line > output/$outfile
    else
        echo $line >> output/$outfile
    fi
done < $1
