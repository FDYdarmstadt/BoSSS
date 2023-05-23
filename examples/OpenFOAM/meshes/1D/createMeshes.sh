#!/usr/bin/env sh

cd small
rm -r constant/polyMesh/*
cp blockMeshDict constant/polyMesh/
blockMesh
cd ..

cd medium
rm -r constant/polyMesh/*
cp blockMeshDict constant/polyMesh/
blockMesh
cd ..

cd large
rm -r constant/polyMesh/*
cp blockMeshDict constant/polyMesh/
blockMesh
cd ..

cp -r small/constant/polyMesh ~/BoSSS-experimental/public/src/L4-application/ExternalBinding/meshes/big/small
cp -r medium/constant/polyMesh ~/BoSSS-experimental/public/src/L4-application/ExternalBinding/meshes/big/medium
cp -r large/constant/polyMesh ~/BoSSS-experimental/public/src/L4-application/ExternalBinding/meshes/big/large
