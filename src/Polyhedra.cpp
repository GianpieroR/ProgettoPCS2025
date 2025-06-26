#include <iostream>
#include <fstream>
#include <sstream>
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{	
	void generaTetraedro(PolyhedralMesh& mesh) {
		
		// VERTICI DEL TETRAEDRO
		double r = sqrt(3.0) / 3.0;
		
		mesh.Cell0DsCoordinates = MatrixXd::Zero(3, 4);
	
		double coords[4][3] = {
			{r,  r,  r},
			{-r, -r,  r},
			{-r,  r, -r},
			{r,  -r, -r}
		};

		for (int i = 0; i < 4; ++i) {
			mesh.Cell0DsCoordinates(0, i) = coords[i][0];
			mesh.Cell0DsCoordinates(1, i) = coords[i][1];
			mesh.Cell0DsCoordinates(2, i) = coords[i][2];
		}
		
		mesh.Cell0DsId = {0,1,2,3};
		
		// LATI DEL TETRAEDRO
		mesh.Cell1DsId = {0,1,2,3,4,5};
	
		mesh.Cell1DsExtrema = MatrixXi::Zero(6, 2);
		mesh.Cell1DsExtrema <<
			0, 1, // lato 0
			1, 2, // lato 1
			2, 0, // lato 2
			0, 3, // lato 3 
			3, 1, // lato 4
			2, 3; // lato 5

		// FACCE DEL TETRAEDRO
		mesh.Cell2DsId = {0, 1, 2, 3};
	
		mesh.Cell2DsVertices = {
			{0, 1, 2}, // faccia 0
			{0, 3, 1}, // faccia 1
			{1, 3, 2}, // faccia 2
			{2, 3, 0}  // faccia 3
		};
	
		mesh.Cell2DsEdges = {
			{0, 1, 2}, // faccia 0 
			{3, 4, 0}, // faccia 1
			{4, 5, 1}, // faccia 2
			{5, 3, 2}  // faccia 3
		};
	
		// TETRAEDRO
		mesh.Cell3DsId = {0};
		mesh.Cell3DsVertices = {0, 1, 2, 3};
		mesh.Cell3DsEdges = {0, 1, 2, 3, 4, 5};
		mesh.Cell3DsFaces = {0, 1, 2, 3};
	
	}
				
	void generaOttaedro(PolyhedralMesh& mesh) {
		
		// VERTICI DELL'OTTAEDRO
		double r = 1.0; 
	
		mesh.Cell0DsCoordinates = MatrixXd::Zero(3, 6);
	
				double coords[6][3] = {
			{ r,  0.0, 0.0},
			{-r,  0.0, 0.0},
			{0.0,  r,  0.0},
			{0.0, -r,  0.0},
			{0.0, 0.0,  r},
			{0.0, 0.0, -r}
		};

		for (int i = 0; i < 6; ++i) {
			mesh.Cell0DsCoordinates(0, i) = coords[i][0];
			mesh.Cell0DsCoordinates(1, i) = coords[i][1];
			mesh.Cell0DsCoordinates(2, i) = coords[i][2];
		}
	
		mesh.Cell0DsId = {0, 1, 2, 3, 4, 5};
		
		// LATI DELL'OTTAEDRO
		mesh.Cell1DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
	
		mesh.Cell1DsExtrema = MatrixXi::Zero(12, 2);
		mesh.Cell1DsExtrema <<
			0, 2,  // Lato 0
			2, 1,  // Lato 1
			1, 3,  // Lato 2
			3, 0,  // Lato 3
			0, 4,  // Lato 4
			2, 4,  // Lato 5
			4, 1,  // Lato 6
			4, 3,  // Lato 7
			1, 5,  // Lato 8
			3, 5,  // Lato 9
			0, 5,  // Lato 10
			2, 5;  // Lato 11
	
		// FACCE DELL'OTTAEDRO
		mesh.Cell2DsId = {0, 1, 2, 3, 4, 5, 6, 7};

		mesh.Cell2DsVertices = {
			{0, 2, 4},  // Faccia 0
			{0, 4, 3},  // Faccia 1
			{3, 4, 1},  // Faccia 2
			{1, 2, 4},  // Faccia 3
			{2, 5, 0},  // Faccia 4
			{2, 1, 5},  // Faccia 5
			{1, 5, 3},  // Faccia 6
			{0, 3, 5}   // Faccia 7
		};
	
		mesh.Cell2DsEdges = {
			{0, 5, 4},  // Faccia 0
			{4, 7, 3},  // Faccia 1
			{7, 6, 2}, // Faccia 2
			{1, 5, 6},  // Faccia 3
			{11, 10, 0}, // Faccia 4
			{1, 8 , 11}, // Faccia 5
			{8, 9, 2},  // Faccia 6
			{3, 9, 10}  // Faccia 7
		};
	
		// OTTAEDRO	
		mesh.Cell3DsId = {0};
		mesh.Cell3DsVertices = {0, 1, 2, 3, 4, 5};
		mesh.Cell3DsEdges = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
		mesh.Cell3DsFaces = {0, 1, 2, 3, 4, 5, 6, 7};
		
	}
	
	void generaIcosaedro(PolyhedralMesh& mesh) {
		// VERTICI DELL'ICOSAEDRO
		double r = 1.0;
    	double phi = (1.0 + sqrt(5.0)) / 2.0;
	
		mesh.Cell0DsCoordinates = MatrixXd::Zero(3, 12);
		
		double norm = sqrt(10.0 + 2.0 * sqrt(5.0)) / 2.0;
				double coords[12][3] = {
			{0.0,       r / norm,  phi / norm},
			{0.0,      -r / norm,  phi / norm},
			{0.0,       r / norm, -phi / norm},
			{0.0,      -r / norm, -phi / norm},
			{r / norm,  phi / norm, 0.0},
			{-r / norm, phi / norm, 0.0},
			{r / norm, -phi / norm, 0.0},
			{-r / norm,-phi / norm, 0.0},
			{phi / norm, 0.0,       r / norm},
			{-phi / norm,0.0,       r / norm},
			{phi / norm, 0.0,      -r / norm},
			{-phi / norm,0.0,      -r / norm}
		};

		for (int i = 0; i < 12; ++i) {
			mesh.Cell0DsCoordinates(0, i) = coords[i][0];
			mesh.Cell0DsCoordinates(1, i) = coords[i][1];
			mesh.Cell0DsCoordinates(2, i) = coords[i][2];
        }
	
		mesh.Cell0DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
	
		// LATI DELL'ICOSAEDRO
		mesh.Cell1DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
					 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
		
		mesh.Cell1DsExtrema = MatrixXi::Zero(30, 2);
		mesh.Cell1DsExtrema << 0,1, 0,4, 0,5, 0,8, 0,9, 1,6, 1,7, 1,8, 1,9, 2,3,
							  2,10, 2,4, 3,6, 3,7, 3,10, 3,11, 4,5, 9,7, 4,8, 4,10,
							  5,2, 5,9, 5,11, 6,8, 6,7, 6,10, 7,11, 8,10, 9,11, 2,11;
		 
		// FACCE DELL'ICOSAEDRO
		mesh.Cell2DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
		
		mesh.Cell2DsVertices = {
			{0, 1, 9},{0, 9, 5},{0, 5, 4},{0, 4, 8},{0, 8, 1},  {1, 6, 8}, {8, 6, 10}, {10, 8, 4}, {10, 4, 2}, {4, 2, 5}, {5, 2, 11}, 
			{5, 11, 9}, {9, 7, 11}, {1, 9, 7}, {1, 7, 6}, {7, 11, 3}, {11, 3, 2}, {2, 3, 10}, {10, 6, 3},{3, 6, 7}
		};
		
		mesh.Cell2DsEdges = {
			{0, 8, 4},{4, 21, 2},{2, 16, 1},{1, 18, 3},{3, 7, 0},{5, 23, 7},{23, 25, 27},{27, 18, 19},{19, 11, 10},{11, 20, 16},  
			{20, 29, 22},{22, 28, 21},{17, 26, 28},{8, 17, 6},{6, 24, 5},{26, 15, 13},{15, 9, 29},{9, 14, 10},{25, 12, 14},  
			{12, 24, 13}    

		};
		
		// ICOSAEDRO
		mesh.Cell3DsId = {0};
		mesh.Cell3DsVertices = {0, 1, 2, 3, 4, 5,6 ,7 ,8, 9, 10, 11};
		mesh.Cell3DsEdges = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
			16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
		mesh.Cell3DsFaces = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
		
	}
}