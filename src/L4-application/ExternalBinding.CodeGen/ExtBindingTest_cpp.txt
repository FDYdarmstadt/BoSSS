// ExtBindingTest.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "BoSSScpp.h"

using namespace BoSSS::Application::ExternalBinding;
using namespace BoSSS::Foundation::Grid;
using namespace BoSSS::Foundation::Grid::Classic;

void InitMeshTest() {
    int nCells = 9;

    
    int __faces[][4] = {
     {1, 5, 21, 17},
     {4, 20, 21, 5},
     {2, 6, 22, 18},
     {5, 21, 22, 6},
     {6, 22, 23, 7},
     {5, 9, 25, 21},
     {8, 24, 25, 9},
     {6, 10, 26, 22},
     {9, 25, 26, 10},
     {10, 26, 27, 11},
     {9, 13, 29, 25},
     {10, 14, 30, 26},
     {12, 28, 29, 13},
     {13, 29, 30, 14},
     {14, 30, 31, 15},
     {0, 16, 20, 4},
     {4, 20, 24, 8},
     {8, 24, 28, 12},
     {3, 7, 23, 19},
     {7, 11, 27, 23},
     {11, 15, 31, 27},
     {0, 1, 17, 16},
     {1, 2, 18, 17},
     {2, 3, 19, 18},
     {0, 4, 5, 1},
     {4, 8, 9, 5},
     {8, 12, 13, 9},
     {1, 5, 6, 2},
     {5, 9, 10, 6},
     {9, 13, 14, 10},
     {2, 6, 7, 3},
     {6, 10, 11, 7},
     {10, 14, 15, 11},
     {16, 17, 21, 20},
     {20, 21, 25, 24},
     {24, 25, 29, 28},
     {17, 18, 22, 21},
     {21, 22, 26, 25},
     {25, 26, 30, 29},
     {18, 19, 23, 22},
     {22, 23, 27, 26},
     {26, 27, 31, 30}
    };

    int nFaces = 42;
    int* faces[42];
    int vertices_per_face[42];
    for (int i = 0; i < nFaces; i++) {
        vertices_per_face[i] = 4;
        faces[i] = __faces[i];
    }

   

    int neighbour[] = {
            1,
            3,
            2,
            4,
            5,
            4,
            6,
            5,
            7,
            8,
            7,
            8
    };
    
    int owner[] = {
            0,
            0,
            1,
            1,
            2,
            3,
            3,
            4,
            4,
            5,
            6,
            7,
            6,
            7,
            8,
            0,
            3,
            6,
            2,
            5,
            8,
            0,
            1,
            2,
            0,
            3,
            6,
            1,
            4,
            7,
            2,
            5,
            8,
            0,
            3,
            6,
            1,
            4,
            7,
            2,
            5,
            8
    };

    double points[][3] = {
     {0,0,0},
     {0.03333333333,0,0},
     {0.06666666667,0,0},
     {0.1,0,0},
     {0,0.03333333333,0},
     {0.03333333333,0.03333333333,0},
     {0.06666666667,0.03333333333,0},
     {0.1,0.03333333333,0},
     {0,0.06666666667,0},
     {0.03333333333,0.06666666667,0},
     {0.06666666667,0.06666666667,0},
     {0.1,0.06666666667,0},
     {0,0.1,0},
     {0.03333333333,0.1,0},
     {0.06666666667,0.1,0},
     {0.1,0.1,0},
     {0,0,0.01},
     {0.03333333333,0,0.01},
     {0.06666666667,0,0.01},
     {0.1,0,0.01},
     {0,0.03333333333,0.01},
     {0.03333333333,0.03333333333,0.01},
     {0.06666666667,0.03333333333,0.01},
     {0.1,0.03333333333,0.01},
     {0,0.06666666667,0.01},
     {0.03333333333,0.06666666667,0.01},
     {0.06666666667,0.06666666667,0.01},
     {0.1,0.06666666667,0.01},
     {0,0.1,0.01},
     {0.03333333333,0.1,0.01},
     {0.06666666667,0.1,0.01},
     {0.1,0.1,0.01}
    };

    int nPoints = 32;
    
    int nInternalFaces = 12;
    int** names = (int**) malloc(sizeof(int*)*3);
    int name1[] {1,2,3,4};
    int name2[] {1,2,3,4,5};
    int name3[] {1,2,3,4,6};
    names[0] = name1;
    names[1] = name2;
    names[2] = name3;

        int patchIDs[] {0,
                                                       1,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2,
                                                       2
        };

        int nameLengths[3] {4, 5, 5};


    BoSSS::Foundation::Grid::OpenFOAMGrid
        grd(nPoints, nCells, nFaces, nInternalFaces, 3, nameLengths, -1,
            faces, vertices_per_face, neighbour, owner, (double*) points,
            names, patchIDs
            );

    int RR;
    RR = grd.TestMethod(7);
    printf("Received in cpp: %d\n", RR);


    GridData* grid1 = grd.GetGridData();
    
}

int main()
{
	printf("Hello cruel world!\n");

    mono_config_parse(NULL);

#ifdef _WIN32
    //char monodir[] = "C:\\Program Files\\Mono";
    char exedir[] = "C:\\Users\\florian\\Documents\\BoSSS-master\\public\\src\\L4-application\\ExternalBinding\\bin\\Debug\\";
#else
    char exedir[] = "./Debug/";
#endif 

    BoSSS::Globals::Init(exedir);
    Initializer MyInit;
    MyInit.BoSSSInitialize();
    printf("init done.\n");

    InitMeshTest();

    MyInit.BoSSSFinalize();
	return 0;
}

