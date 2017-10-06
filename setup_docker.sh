FUSTR_dir=$(pwd)
target=${1%/}

cp -r $target $FUSTR_dir
echo cp -r  $target $FUSTR_dir/$(basename $target)
echo $FUSTR_dir
# pwd
# work_dir=$FUSTR_dir/$(basename $target)
echo $work_dir
ls $work_dir
docker build -t  fustr --build-arg package=$(basename $target)  $FUSTR_dir --no-cache
