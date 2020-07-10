#!/usr/bin/bash
echo "GenProcessor Intiated"
echo "Make sure that GenProcessor.sh is located in the same directory as the annovar perl files: convert2annovar.pl AND table.pl"
echo "Make sure the following directories found in ALS_consortium are correctly labeled as shown below (M-D-Y means today's date in ex: 6/3/2020 format):"
echo "ALS_consortium"
echo "	M-D-Y"
echo " 	M-D-Y_Annotated"
echo "	AnnovarReady"
echo "		M-D-Y_annovar_ready"
echo "Would you like to proceed? [y/n]"
read proceed
yesholder="y"
if [[ "$proceed" == "$yesholder" ]]; then
	echo "Please provide the file path of ALS_consortium but not including ALS_consortium. (Start with /mnt/c/ where c is your c drive AND make sure to use forward slashes)"
	read ALSconsortiumFilePath
	dateHolder="$(date +"%-m-%-d-%Y")"
	echo "Please provide the file path of the AnnovarIDcontainer.txt in the same format as above, but include AnnovarIDcontainer.txt"
	read input
	echo "Please provide the file path of the QualityMetricsSingle.py in the same format as above, but include QualityMetricsSingle.py"
	read pythonfilepath
	a=1
	while IFS= read -r line
	do
		echo "$line"
		echo "$a"
		a=$((a+1))
		perl convert2annovar.pl -format vcf4 ""$ALSconsortiumFilePath"ALS_consortium/"$dateHolder"/CGND-HDA-"$line".filtered_acmg_als.vcf"> ""$ALSconsortiumFilePath"ALS_consortium/AnnovarReady/"$dateHolder"_annovar_ready/CGND-HDA-"$line".filtered_acmg_als.avinput"
		perl table_annovar.pl ""$ALSconsortiumFilePath"ALS_consortium/AnnovarReady/"$dateHolder"_annovar_ready/CGND-HDA-"$line".filtered_acmg_als.avinput" humandb/ -buildver hg38 -remove -protocol refGene,knownGene,exac03,gnomad211_genome,gnomad211_exome,abraom,dbscsnv11,revel,clinvar_20190305,mcap,dbnsfp35a,1000g2015aug_all,1000g2015aug_eur,abraom -operation g,g,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -csvout -out ""$ALSconsortiumFilePath"ALS_consortium/"$dateHolder"_Annotated/hg38_CGND-HDA-"$line".filtered_acmg_als.avinput_updated"
		python pythonfilepath $line $ALSconsortiumFilePath $dateHolder
	done < "$input"
	grep -H "[A-Za-z0-9]" "$ALSconsortiumFilePath"ALS_consortium/"$dateHolder"_Annotated/*.csv | sed "s|:|,|" > "$ALSconsortiumFilePath"ALS_consortium/"$dateHolder"_Annotated/"$dateHolder"_Annotated-Combined.csv
	grep -H "[A-Za-z0-9]" "$ALSconsortiumFilePath"ALS_consortium/"$dateHolder"_Annotated/Original_Annotation/*.csv | sed "s|:|,|" > "$ALSconsortiumFilePath"ALS_consortium/"$dateHolder"_Annotated/Original_Annotation/"$dateHolder"_Annotated-Combined.csv
fi
