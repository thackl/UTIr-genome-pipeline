tool=sylph
db=databases/$tool

if [ -d "$db" ]; then
    echo "Existing $tool database detect at $db."
    echo "Skipping download. Remove folder to trigger redownload"
    exit;
fi;

echo "Downloading $tool database"
mkdir -p $db
cd $db
wget http://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r220-c200-dbv1.syldb
