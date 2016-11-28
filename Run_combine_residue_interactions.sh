# Run as $ bash script.sh *firstres* *lastres* 

echo "First residue is $1 and last is $2"
rep=1
while [ $rep -lt $1 ];
do
	python Combine_residue_interactions.py $rep $1 $2
	rep=$((rep+1))

done
