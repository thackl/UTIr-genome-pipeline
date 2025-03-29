tool=gtdbtk
db=databases/$tool

if [ -d "$db" ]; then
    echo "Existing $tool database detect at $db."
    echo "Skipping download. Remove folder to trigger redownload"
    exit;
fi;

echo "Downloading $tool database"
download-db.sh $db
