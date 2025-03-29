tool=bakta
db=databases/$tool

if [ -d "$db" ]; then
    echo "Existing $tool database detect at $db."
    echo "Skipping download. Remove folder to trigger redownload"
    exit;
fi;

echo "Downloading $tool database"
bakta_db download --output $db --type full
