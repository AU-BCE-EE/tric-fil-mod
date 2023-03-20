
perl -pi -w -e 's/([^ ])\+([^ ])/$1 + $2/g' *.py
perl -pi -w -e 's/([^ ])\*([^ ])/$1 * $2/g' *.py
perl -pi -w -e 's/([^ ])\/([^ ])/$1 \/ $2/g' *.py
perl -pi -w -e 's/ \* \*/**/g' *.py
perl -pi -w -e 's/([^ ])=([^ ])/$1 = $2/g' *.py
perl -pi -w -e 's/([^ ]),([^ ])/$1, $2/g' *.py
perl -pi -w -e 's/([^>^<^ ^=^!])=([^ ^=])/$1 = $2/g' *.py
perl -pi -w -e 's/([^ ])==([^ ])/$1 == $2/g' *.py
perl -pi -w -e 's/([^ ])!=([^ ])/$1 != $2/g' *.py 
