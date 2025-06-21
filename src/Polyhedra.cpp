#include <iostream>
#include <fstream>
#include <sstream>
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{	
	void generateTetrahedron(PolyhedralMesh& mesh) {
		
		// VERTICI
		double r = sqrt(3.0) / 3.0;
		
		mesh.Cell0DsCoordinates = MatrixXd::Zero(3, 4);
	
		mesh.Cell0DsCoordinates <<
        r, -r, -r,  r,
        r, -r,  r, -r,
        r,  r, -r, -r;
		
		mesh.Cell0DsId = {0,1,2,3};
		
		// LATI
		mesh.Cell1DsId = {0,1,2,3,4,5};
	
		mesh.Cell1DsExtrema = MatrixXi::Zero(6, 2);
		mesh.Cell1DsExtrema <<
			0, 1, // lato 0
			1, 2, // lato 1
			2, 0, // lato 2
			0, 3, // lato 3 
			3, 1, // lato 4
			2, 3; // lato 5

		// FACCE
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
	
		// POLIEDRO
		mesh.Cell3DsId = {0};
		mesh.Cell3DsVertices = {0, 1, 2, 3};
		mesh.Cell3DsEdges = {0, 1, 2, 3, 4, 5};
		mesh.Cell3DsFaces = {0, 1, 2, 3};
	
	}
	
	void generateCube(PolyhedralMesh& mesh) {
		
		// VERTICI
		double r = sqrt(3.0) / 3.0;
	
		mesh.Cell0DsCoordinates = MatrixXd::Zero(3, 8); 
	
		mesh.Cell0DsCoordinates <<
		 r, -r, -r,  r,  r, -r, -r,  r,
		 r,  r, -r, -r,  r,  r, -r, -r,
		 r,  r,  r,  r, -r, -r, -r, -r;
		
		mesh.Cell0DsId = {0,1,2,3,4,5,6,7};
		
		// LATI
		mesh.Cell1DsId = {0,1,2,3,4,5,6,7,8,9,10,11};
	
		mesh.Cell1DsExtrema = MatrixXi::Zero(12, 2);
		mesh.Cell1DsExtrema <<
			0, 1,  // Lato 0
			1, 2,  // Lato 1
			2, 3,  // Lato 2
			3, 0,  // Lato 3
			4, 5,  // Lato 4
			5, 6,  // Lato 5
			6, 7,  // Lato 6
			7, 4,  // Lato 7
			0, 4,  // Lato 8
			1, 5,  // Lato 9
			2, 6,  // Lato 10
			3, 7;  // Lato 11
			
		// FACCE
		mesh.Cell2DsId = {0, 1, 2, 3, 4, 5};

		mesh.Cell2DsVertices = {
			{0, 1, 2, 3},  
			{4, 5, 6, 7},  
			{0, 1, 5, 4},  
			{1, 2, 6, 5},  
			{2, 3, 7, 6},
			{3, 0, 4, 7} 
		};
		
		mesh.Cell2DsEdges = {
			{0, 1, 2, 3},
			{4, 5, 6, 7},
			{0, 9, 4, 8},
			{1, 10, 5, 9},
			{2, 11, 6, 10},
			{3, 8, 7, 11}
		};
		
		// POLIEDRO 
		mesh.Cell3DsId = {0};
		mesh.Cell3DsVertices = {0, 1, 2, 3, 4, 5, 6, 7};
		mesh.Cell3DsEdges = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
		mesh.Cell3DsFaces = {0, 1, 2, 3, 4, 5};
		
	}
			
	void generateOctahedron(PolyhedralMesh& mesh) {
		
		// VERTICI
		double r = 1.0; 
	
		mesh.Cell0DsCoordinates = MatrixXd::Zero(3, 6);
	
		mesh.Cell0DsCoordinates <<
		 r, -r, -r,  r,  r, -r, -r,  r,
		 0,  0,  0,  0,  r,  r, -r, -r,
		 0,  0,  0,  0,  r, -r, -r, -r;
	
		mesh.Cell0DsId = {0, 1, 2, 3, 4, 5};
		
		// LATI
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
	
		// FACCE
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
	
		// POLIEDRO	
		mesh.Cell3DsId = {0};
		mesh.Cell3DsVertices = {0, 1, 2, 3, 4, 5};
		mesh.Cell3DsEdges = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
		mesh.Cell3DsFaces = {0, 1, 2, 3, 4, 5, 6, 7};
		
	}
	
	void generateDodecahedron(PolyhedralMesh& mesh) { // probabilmente da togliere
		
		// VERTICI
		mesh.Cell0DsCoordinates = MatrixXd::Zero(3, 20);
		double phi = (1.0 + sqrt(5.0)) / 2.0;
		double invPhi = 1.0 / phi;
		double invNorm = 1.0 / sqrt(3.0);
		
		Vector3d points[] = {
			{  1,  1,  1}, {  1,  1, -1}, {  1, -1,  1}, {  1, -1, -1},
			{ -1,  1,  1}, { -1,  1, -1}, { -1, -1,  1}, { -1, -1, -1},
			{  0,  invPhi,  phi}, {  0,  invPhi, -phi}, {  0, -invPhi,  phi}, {  0, -invPhi, -phi},
			{ invPhi,  phi,  0}, { -invPhi,  phi,  0}, { invPhi, -phi,  0}, { -invPhi, -phi,  0},
			{  phi,  0,  invPhi}, {  phi,  0, -invPhi}, { -phi,  0,  invPhi}, { -phi,  0, -invPhi}
		};

		for (int i = 0; i < 20; ++i) {
			mesh.Cell0DsCoordinates.col(i) = invNorm * points[i];
		}
	
		mesh.Cell0DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
		
		// LATI
		mesh.Cell1DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
					 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
	
		mesh.Cell1DsExtrema = MatrixXi::Zero(30, 2);
		mesh.Cell1DsExtrema <<
			0, 8,
			0, 12,
			0, 16,
			1, 9,
			1, 12,
			1, 17,
			2, 10,
			2, 14,
			2, 16,
			3, 11,
			3, 14,
			3, 17,
			4, 8,
			4, 13,
			4, 18,
			5, 9,
			5, 13,
			5, 19,
			6, 10,
			6, 15,
			6, 18,
			7, 11,
			7, 15,
			7, 19,
			8, 12,
			9, 12,
			10, 14,
			11, 14,
			13, 18,
			15, 18;
	
		// FACCE 
		mesh.Cell2DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
		
		 mesh.Cell2DsVertices = {
			{0, 1, 2, 3, 4},  // Faccia 0
			{1, 5, 6, 7, 2},  // Faccia 1
			{2, 7, 8, 19, 3},  // Faccia 2
			{10, 6, 7, 8, 9},  // Faccia 3
			{9, 8, 19, 17, 16},  // Faccia 4
			{3, 19, 17, 18, 4},  // Faccia 5
			{16, 17, 18, 14, 15},  // Faccia 6
			{15, 14, 13, 12, 11},   // Faccia 7
			{0, 4, 18, 14, 13},  // Faccia 8
			{13, 12, 5, 1, 0},  // Faccia 9
			{11, 12, 5, 6, 10},  // Faccia 10
			{9, 16, 15, 11, 10},  // Faccia 11
		};
	
		mesh.Cell2DsEdges = {
			{0, 1, 2, 3, 4},  // Faccia 0
			{5, 6, 7, 8, 1},  // Faccia 1
			{8, 29, 27, 28, 2},  // Faccia 2
			{9, 7, 29, 10, 11},  // Faccia 3
			{10, 27, 25, 24, 26},  // Faccia 4
			{28, 25, 23, 22, 3},  // Faccia 5
			{24, 23, 21, 17, 20},  // Faccia 6
			{17, 18, 15, 14, 16},   // Faccia 7
			{4, 22, 21, 18, 19},  // Faccia 8
			{15, 13, 5, 0, 19},  // Faccia 9
			{14, 13, 6, 9, 12},  // Faccia 10
			{26, 20, 16, 12, 11},  // Faccia 11
		};
	
		// POLIEDRO
		mesh.Cell3DsId = {0};
		mesh.Cell3DsVertices = {0, 1, 2, 3, 4, 5,6 ,7 ,8, 9, 
			10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
		mesh.Cell3DsEdges = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
			16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
		mesh.Cell3DsFaces = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
		
	}
	
	void generateIcosahedron(PolyhedralMesh& mesh) {
		// VERTICI
		double r = 1.0;
    	double phi = (1.0 + sqrt(5.0)) / 2.0;
	
		mesh.Cell0DsCoordinates = MatrixXd::Zero(3, 12);
		
		double norm = sqrt(10.0 + 2.0 * sqrt(5.0)) / 2.0;
		Vector3d points[] = {
			{0.0,  r,  phi}, {0.0, -r,  phi}, {0.0,  r, -phi}, {0.0, -r, -phi},
			{r,  phi, 0.0}, {-r,  phi, 0.0}, {r, -phi, 0.0}, {-r, -phi, 0.0},
			{ phi, 0.0,  r}, {-phi, 0.0,  r}, { phi, 0.0, -r}, {-phi, 0.0, -r}
		};

		for (int i = 0; i < 12; ++i) {
			mesh.Cell0DsCoordinates.col(i) = points[i] / norm;
		}
		mesh.Cell0DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
	
		// LATI
		mesh.Cell1DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
					 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
		
		mesh.Cell1DsExtrema = MatrixXi::Zero(30, 2);
		mesh.Cell1DsExtrema <<
		  0, 1,   // 0
		  0, 4,   // 1
		  0, 5,   // 2
		  0, 8,   // 3
		  0, 9,   // 4
		  1, 6,   // 5
		  1, 7,   // 6
		  1, 8,   // 7
		  1, 9,   // 8
		  2, 3,   // 9
		  2, 10,   //10
		  2, 4,   //11
		  3, 6,   //12
		  3, 7,   //13
		  3, 10,  //14
		  3, 11,  //15
		  4, 5,   //16
		  9, 7,   //17
		  4, 8,   //18
		  4, 10,  //19
		  5, 2,   //20
		  5, 9,   //21
		  5, 11,  //22
		  6, 8,   //23
		  6, 7,   //24
		  6, 10,  //25
		  7, 11,  //26
		  8, 10,  //27
		  9, 11,  //28
		  2, 11;  //29	
		 
		// FACCE
		mesh.Cell2DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
		
		mesh.Cell2DsVertices = {
			{0, 1, 9},
			{0, 9, 5},
			{0, 5, 4},
			{0, 4, 8},
			{0, 8, 1},  
			{1, 6, 8}, 
			{8, 6, 10}, 
			{10, 8, 4}, 
			{10, 4, 2}, 
			{4, 2, 5}, 
			{5, 2, 11}, 
			{5, 11, 9}, 
			{9, 7, 11}, 
			{1, 9, 7}, 
			{1, 7, 6}, 
			{7, 11, 3}, 
			{11, 3, 2}, 
			{2, 3, 10}, 
			{10, 6, 3},
			{3, 6, 7}
		};
		
		mesh.Cell2DsEdges = {
			{0, 8, 4},     
			{4, 21, 2},    
			{2, 16, 1},    
			{1, 18, 3},    
			{3, 7, 0},     
			{5, 23, 7},    
			{23, 25, 27},  
			{27, 18, 19}, 
			{19, 11, 10}, 
			{11, 20, 16},  
			{20, 29, 22},  
			{22, 28, 21},  
			{17, 26, 28},  
			{8, 17, 6},    
			{6, 24, 5},    
			{26, 15, 13},   
			{15, 9, 29},  
			{9, 14, 10},  
			{25, 12, 14},  
			{12, 24, 13}    

		};
		
		// POLIEDRO
		mesh.Cell3DsId = {0};
		mesh.Cell3DsVertices = {0, 1, 2, 3, 4, 5,6 ,7 ,8, 9, 10, 11};
		mesh.Cell3DsEdges = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 
			16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
		mesh.Cell3DsFaces = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
		
	}
}