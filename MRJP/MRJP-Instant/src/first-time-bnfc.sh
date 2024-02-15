# DON'T USE
# THIS WILL OVERRIDE MAKEFILE

if [ -z "$IKNOWWHATIAMDOING" ]
then
    echo "Please don't use this script. Use make instead."
    echo "If you really want to use this script, set the environment variable IKNOWWHATIAMDOING to any value."
    exit 1
fi


bnfc -m -d --haskell -o . Instant.cf --functor
