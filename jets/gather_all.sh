filenames=`ls *.root`
make gather
for filename in $filenames
do
	./gather $filename
done
