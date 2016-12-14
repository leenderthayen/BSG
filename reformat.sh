#!/bin/bash
for file in "$@"
do
	echo "Reformatting $file"
	clang-format-3.6 -style=Google $file > $file.format && mv $file.format $file
done
