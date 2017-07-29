## folder made by date, in case we make major changes to exac-regions.py or resid-plot.py ##
date=2017_04_14 # default date value
if [ -s coveragecut.txt ]; then
    rm coveragecut.txt segdupscut.txt selfchaincut.txt #deletioncut.txt
fi
while getopts ":v:d:scnwf:" opt; do
    case $opt in
        v)
            echo "-version (date) was triggered, input: $OPTARG" >&2
            date=$OPTARG
            ;;
        s)
            echo "-synonymous variant density input into the model was triggered" >&2
            syn="-s"
            s="-synonymous"
            ;;
        c)
            echo "-CpG density input into the model was triggered" >&2
            cpg="-c"
            c="-cpg"
            ;;
        n)
            echo "-no singleton input into the model was triggered" >&2
            ns="-n"
            n="-nosingletons"
            ;;
        f)
            echo "-file input into the model was specified, input: $OPTARG" >&2
            file=$OPTARG
            ;;
        w)
            echo "-variant flags incorporated into the model" >&2
            var="-w"
            w="-novariant"
            ;;
        d)
            echo "depth and depth cutoff passed into the model" >&2
            depth+=("$OPTARG")
            echo "-depth input: '${depth[@]}'" >&2
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done
### for when something is absolutely required in this script
#if [ -z "$file" ]; then
#  echo "-f [option] is required"
#  exit
#fi
mkdir -p results/$date/
## generates regions and residuals files ##
python exac-regions.py $ns $var -c "data/exacv2.chr{chrom}.cov.txt.gz" -e data/Homo_sapiens.GRCh37.75.gtf.gz -x data/gnomad-vep-vt.vcf.gz -d "${depth[0]}" -l "${depth[1]}" > results/$date/exac-regions$n$w.txt # added $file as a placeholder for now, so we don't always hard code files
#python exac-regions.py $ns $var -c "data/Panel.chr{chrom}.coverage.txt.gz" -e data/Homo_sapiens.GRCh37.75.gtf.gz -x $DATA/ExAC.r1.vt.vep.vcf.gz -d "${depth[0]}" -l "${depth[1]}" > results/$date/exac-regions$n$w.txt # added $file as a placeholder for now, so we don't always hard code files
python resid-plot.py $ns $syn $cpg $var -f results/$date/exac-regions$n$w.txt > results/$date/resids$c$s$n$w.txt
cat <(head -1 results/$date/resids$c$s$n$w.txt) <(sed '1d' results/$date/resids$c$s$n$w.txt | sort -k12,12nr) > /tmp/residsort$c$s$n$w.txt
python weightpercentile.py /tmp/residsort$c$s$n$w.txt > results/$date/weightedresiduals$c$s$n$w.txt
