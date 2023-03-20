make clean
if [ "$1" = "fill" ]
then
    make fill
    ./fill
    make clean
elif [ "$1" = "gen" ]
then
    make gen
    ./gen
    make clean
else
    echo "invalid option"
    exit 1
fi
