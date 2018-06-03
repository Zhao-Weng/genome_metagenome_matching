truthTablePath="../strains2_training_metagenomes/truthTable.txt"
averageFile="../strains2_training_metagenomes/average.txt"
touch $truthTablePath
count=0
for genomeFile in [C-E]*.fna
do 
	genomeName=${genomeFile%.*}
	echo -n "$genomeName " >>  $truthTablePath
	# echo $genomeFile
	suffix="${genomeFile//[!0-9]/}"

	prefix="${genomeFile:0:2}"
	dirName=$prefix$suffix
	echo $dirName
	mkdir $dirName
	bowtie2-build $genomeFile $dirName
	mv *.bt2 $dirName
	(( count++ ))

	cd ../strains2_training_metagenomes

	for idx in 4
	do
		sum=0 
		for nextI in 1 2
		do
			metaFile=$(find . -name "*$idx\_$nextI.fastq")
			echo $metaFile
			outputFile=$dirName-$idx$nextI.sam
			#echo $outputFile
			bowtie2 -p 8 -x ../strains2_training_genomes/$dirName/$dirName $metaFile 2>stderr.log 1>$outputFile --local
			rate=$(cat stderr.log | tail -1 | grep -o '[0-9].[0-9][0-9]')
			rm stderr.log $outputFile
			sum=$(echo $sum $rate | awk '{ printf "%f", $1 + $2 }')
			echo $sum
		done
		average=$(echo $sum 2 | awk '{ printf "%f", $1 / $2 }')
		echo "$metaFile $genomeName $average" >> average.txt   # write to average file
		echo  >> average.txt  
		cmp=$(echo $average'>'5 | bc -l)
		echo -n "$cmp	" >>  truthTable.txt
		
	done
	echo >> truthTable.txt    # next line
	cd ../strains2_training_genomes
done
echo $count


