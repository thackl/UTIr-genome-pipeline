tool=checkm2
db=resources/databases/$tool

if [ -d "$db" ]; then
    echo "Existing $tool database detect at $db."
    echo "Skipping download. Remove folder to trigger redownload"
    exit;
fi;

echo "Downloading $tool database"
checkm2 database --download --path $db
