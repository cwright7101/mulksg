VERSION="$(cat assembler/VERSION)"
TARGET_DIR=SPAdes-$VERSION
rm -rf $TARGET_DIR
SRC_DIR=$TARGET_DIR/src
mkdir -p $SRC_DIR

cp -r assembler/src/common $SRC_DIR/
cp -r assembler/src/projects $SRC_DIR/
cp -r assembler/src/include $SRC_DIR/
cp -r assembler/src/cmake $SRC_DIR/
cp -r assembler/src/spades_pipeline $SRC_DIR/
cp assembler/src/CMakeLists.txt $SRC_DIR/

cp -r assembler/configs $TARGET_DIR/configs
cp -r assembler/ext $TARGET_DIR/ext
rm -r $TARGET_DIR/ext/include/cute
rm -r $TARGET_DIR/ext/include/teamcity_boost

# cleaning .pyc and .pyo
rm -f $SRC_DIR/spades_pipeline/*.pyc
rm -f $SRC_DIR/spades_pipeline/*.pyo
rm -f $TARGET_DIR/ext/include/python_libs/joblib/*.pyc
rm -f $TARGET_DIR/ext/include/python_libs/joblib/*.pyo

cp -r assembler/test_dataset $TARGET_DIR/test_dataset
cp -r assembler/test_dataset_truspades $TARGET_DIR/test_dataset_truspades
cp -r assembler/test_dataset_plasmid $TARGET_DIR/test_dataset_plasmid
cp assembler/LICENSE $TARGET_DIR/
cp README.md $TARGET_DIR/
cp assembler/VERSION $TARGET_DIR/
cp assembler/spades.py $TARGET_DIR/
cp assembler/mulksg.py $TARGET_DIR/
cp assembler/rnaspades.py $TARGET_DIR/
cp assembler/metaspades.py $TARGET_DIR/
cp assembler/plasmidspades.py $TARGET_DIR/
cp assembler/truspades.py $TARGET_DIR/
cp assembler/spades_compile.sh $TARGET_DIR/
cp assembler/spades_init.py $TARGET_DIR/
cp assembler/mulksg_init.py $TARGET_DIR/
cp assembler/manual.html $TARGET_DIR/
cp assembler/truspades_manual.html $TARGET_DIR/
cp assembler/rnaspades_manual.html $TARGET_DIR/
cp assembler/changelog.html $TARGET_DIR/
cp assembler/GPLv2.txt $TARGET_DIR/

cd $TARGET_DIR/
touch src/CMakeListsInternal.txt
find . -name ".?*" | xargs rm -r

cd ..

tar -pczf SPAdes-$VERSION.tar.gz SPAdes-$VERSION
rm -r SPAdes-$VERSION
