# BEGIN: Test
line_count=$(wc -l < /gpfs/home/jd5457/naiklink/jd5457/.hidden-libraries/multiome-processing/routes/filter_and_cluster.sh)
if [[ $line_count -gt 1 ]]; then
    echo "The code has more than one line."
else
    echo "The code has only one line."
fi
# END: Test