mkdir -p build
cd build

cmake .. && make

./test/run_tests

mkdir -p coverage

COV_FILE="coverage/coverage.info"
COV_DIR=".coverage"
COV_SITE=$COV_DIR/"index.html"
RM_PREFIX="${PWD%/*}/src"

lcov -c -d .. -o $COV_FILE
lcov -r $COV_FILE "*Oersted/lib/*" "*Oersted/test/*" "/usr/include/*" "/usr/lib/*" -o $COV_FILE
genhtml $COV_FILE -o $COV_DIR -p $RM_PREFIX
google-chrome $COV_SITE
