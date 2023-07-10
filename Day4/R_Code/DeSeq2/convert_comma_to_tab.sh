bash 

perl -wnlp -e 's/\t/,/g;' "$1" > "$2"
