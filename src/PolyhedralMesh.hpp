#pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary {

struct PolyhedralMesh
{
	// VERTICI
    vector<unsigned int> Cell0DsId = {}; // id celle 0D
    Eigen::MatrixXd Cell0DsCoordinates = {}; // coordinate celle 0D
	vector<vector<unsigned int>> Cell0DsFlag= {}; // flag celle 0D
	vector<unsigned int> Cell0DsMarker = {}; // marker celle 0D

    // LATI/SPIGOLI
    vector<unsigned int> Cell1DsId = {}; // id celle 1D
    MatrixXi Cell1DsExtrema = {}; // id dei vertici (partenza, arrivo) celle 1D
    vector<unsigned int> Cell1DsFlag= {}; // flag celle 1D
    vector<bool> Cell1DsOriginalFlag = {}; 
    vector<unsigned int> Cell1DsMarker = {}; // marker celle 1D
    
    // FACCE
    vector<unsigned int> Cell2DsId = {}; // id celle 2D
    vector<vector<unsigned int>> Cell2DsVertices = {}; // id dei vertici celle 2D
    vector<vector<unsigned int>> Cell2DsEdges = {}; // id dei lati celle 2D
    
    // POLIEDRI
    unsigned int Cell3DsId = 0; // id celle 3D
   	unsigned int NumCells0Ds = 0;
   	unsigned int NumCells1Ds = 0;
   	unsigned int NumCells2Ds = 0;
    vector<unsigned int> Cell3DsVertices = {}; // id dei vertici celle 3D
    vector<unsigned int> Cell3DsEdges = {}; // id dei lati celle 3D
    vector<unsigned int> Cell3DsFaces = {}; // id delle facce celle 3D   

};

}