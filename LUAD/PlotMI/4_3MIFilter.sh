for i in $(ls *.sort)
do
	Rscript 4_3MIFilterAut.R $i
done
