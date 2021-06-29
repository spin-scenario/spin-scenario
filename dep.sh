#! / bin/sh
rm -r ./3rd-party/linenoise-ng
rm -r ./3rd-party/yacas

git submodule sync
git submodule init
git submodule update

cd ./3rd-party/eigen
git checkout 3.3.5
cd ../../

cd ./3rd-party/lua
mkdir build
cd build
cmake ..
make
cd ../../../

cd ./3rd-party/linenoise-ng
git checkout v1.0.1
sed -i '80d' CMakeLists.txt # delete STATIC option.
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=ON ..
make
cd ../../../

cd ./3rd-party/yacas
git checkout v1.6.1
# COPY THE SCRIPTS FOR YACAS.
mkdir scripts ../../share/spin-scenario/config/yacas/
cp -r scripts ../../share/spin-scenario/config/yacas/scripts 
cd cyacas/libyacas
sed -i '1i cmake_minimum_required(VERSION 2.8)' CMakeLists.txt
sed -i '99c install (TARGETS libyacas LIBRARY DESTINATION lib ARCHIVE DESTINATION lib RUNTIME DESTINATION bin COMPONENT app)' CMakeLists.txt
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -std=c++11" -DBUILD_SHARED_LIBS=ON ..
make -j8
cd ../../../../../


# how to install tensorflow into system: https://tensorflow.google.cn/install/pip  sudo -H pip3 install tensorflow
# uncomment the following if TensorFlow v10 has been install in the system.
#cd ./share/spin-scenario/tf_files/matrix_exp_op
#cmake .
#make
#cd ../../../../

#cd ./3rd-party/tensorflow
#git checkout r1.10 # http://releases.bazel.build/0.15.0/release/index.html
#./configure
#bazel build -c opt --jobs=24 -k //tensorflow:libtensorflow_cc.so
#cd ../../
