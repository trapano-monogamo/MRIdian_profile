#!/bin/bash

for file in $1/**/*
do
		# remove non-utf-8 chars
		iconv -f "windows-1252" -t "UTF-8" -c "$file" -o "$file.utf8" > /dev/null
		cat "$file.utf8" > "$file"
		rm "$file.utf8"

		# remove trailing return chars (\r)
		tr -d '\15\32' < "$file" > "$file.temp"
		cat "$file.temp" > "$file"
		rm "$file.temp"
done
