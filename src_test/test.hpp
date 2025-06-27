#pragma once

#include <iostream>
#include <vector>
#include "Eigen/Eigen"

#include <gtest/gtest.h>
#include "PolyhedralMesh.hpp"
#include "Utils.hpp"

using namespace Eigen;
using namespace PolyhedralLibrary;

namespace PolyhedraTest {

// Dimension
TEST(TestPolyedra, TestComputePolyhedronVEF)
{
	int b = 2;
	int c = 0;
	
	// caso 1: p=3, q=3
	int q1 = 3;
	vector<int> expected1 = {10, 24, 16};
	vector<int> result1 = CalcoloVEFPoliedro(q1, b, c);
	EXPECT_EQ(expected1, result1);
	
	// caso 2: p=3, q=4
	int q2 = 4;
	vector<int> expected2 = {18, 48, 32};
	vector<int> result2 =CalcoloVEFPoliedro(q2, b, c);
	EXPECT_EQ(expected2, result2);
	
    // caso 3: p=3, q=5
	int q3 = 5;
	vector<int> expected3 = {42, 120, 80};
	vector<int> result3 =CalcoloVEFPoliedro(q3, b, c);
	EXPECT_EQ(expected3, result3);
}

TEST(TestPolyedra, TestCalculateDuplicated)
{
	int b = 2;
	int c = 0;
	
	// caso 1: p=3, q=3
	int q1 = 3;
	vector<int> dimension1 = {10, 24, 16};
	vector<int> expected1 = {24, 36, 16};
	vector<int> result1 = CalcoloDuplicato(q1, b, c, dimension1);
	EXPECT_EQ(expected1, result1);
	
	// caso 2: p=3, q=4
	int q2 = 4;
	vector<int> dimension2 = {18, 48, 32};
	vector<int> expected2 = {48, 72, 32};
	vector<int> result2 = CalcoloDuplicato(q2, b, c, dimension2);
	EXPECT_EQ(expected2, result2);
	
	// caso 3: p=3, q=5
	int q3 = 5;
	vector<int> dimension3 = {42, 120, 80};
	vector<int> expected3 = {120, 180 , 80};
	vector<int> result3 = CalcoloDuplicato(q3, b, c, dimension3);
	EXPECT_EQ(expected3, result3);
}


// Triangulation
TEST(TestPolyedra, TestTriangulationTetrahedron)
{
	PolyhedralMesh meshExpected;
	PolyhedralMesh meshFinal;

	// VERTICI
	meshExpected.Cell0DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

	meshExpected.Cell0DsCoordinates = MatrixXd(3, 10);
	
	meshExpected.Cell0DsCoordinates.col(0)  <<  0, 0, 0.57735;
	meshExpected.Cell0DsCoordinates.col(1)  <<  -0.57735, -0.57735,  0.57735;
	meshExpected.Cell0DsCoordinates.col(2)  <<  0, -0.57735, 0;
	meshExpected.Cell0DsCoordinates.col(3)  <<  -0.57735, 0, 0;
	meshExpected.Cell0DsCoordinates.col(4)  <<  -0.57735,  0.57735, -0.57735;
	meshExpected.Cell0DsCoordinates.col(5)  <<  0, 0, -0.57735;
	meshExpected.Cell0DsCoordinates.col(6)  <<  0, 0.57735, 0;
	meshExpected.Cell0DsCoordinates.col(7)  <<  0.57735,  -0.57735, -0.57735;
	meshExpected.Cell0DsCoordinates.col(8)  <<  0.57735, 0, 0;
	meshExpected.Cell0DsCoordinates.col(9)  <<  0.57735, 0.57735, 0.57735;

	meshExpected.Cell0DsFlag = {};

	// LATI/SPIGOLI
	meshExpected.Cell1DsId ={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};

	meshExpected.Cell1DsExtrema = MatrixXi(24, 2);
	meshExpected.Cell1DsExtrema << 
		0, 6,
		0, 3,
		3, 6,
		0, 8,
		0, 9,
		2, 8,
		0, 2,
		0, 1,
		1, 2,
		2, 3,
		1, 3,
		2, 7,
		2, 5,
		3, 5,
		3, 4,
		4, 5,	
		5, 6,
		4, 6,
		5, 7,
		7, 8,
		5, 8,
		6, 8,
		8, 9,
		6, 9;

	meshExpected.Cell1DsFlag = {};
	
	meshExpected.Cell1DsOriginalFlag = {1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0};
	
	// FACCE
	meshExpected.Cell2DsId = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

	meshExpected.Cell2DsVertices = {
		{9, 0, 6},
		{0, 1, 3},
		{0, 3, 6},
		{6, 3, 4},
		{9, 8, 0},
		{8, 7, 2},
		{8, 2, 0},
		{0, 2, 1},
		{1, 2, 3},
		{2, 7, 5},
		{2, 5, 3},
		{3, 5, 4},
		{4, 5, 6},
		{5, 7, 8},
		{5, 8, 6},
		{6, 8, 9}
	};

	meshExpected.Cell2DsEdges = {
		{4, 0, 23},
		{7, 10, 1},
		{1, 2, 0},
		{2, 14, 17},
		{22, 3, 4},
		{19, 11, 5},
		{5, 6, 3},
		{6, 8, 7},
		{8, 9, 10},
		{11, 18, 12},
		{12, 13, 9},
		{13, 15, 14},
		{15, 16, 17},
		{18, 19, 20},
		{20, 21, 16},
		{21, 22, 23}
	};

	// POLIEDRI
	meshExpected.Cell3DsId = {0}; 
	meshExpected.NumCells0Ds = {10};
	meshExpected.NumCells1Ds = {24};
	meshExpected.NumCells2Ds = {16};
	meshExpected.Cell3DsVertices = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	meshExpected.Cell3DsEdges = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};
	meshExpected.Cell3DsFaces = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

	
	// mesh ottenuta utilizzando la funzione di triangolazione
	PolyhedralMesh meshTriangulated;
	PolyhedralMesh mesh;
	generaTetraedro(mesh);
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	vector<int> dimension =CalcoloVEFPoliedro(q, b, c);
	vector<int> dimensionDuplicated = CalcoloDuplicato(q, b, c, dimension);
	triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
	RimuoviLatiDuplicati(meshTriangulated);
	RimuoviVerticiDuplicati(meshTriangulated);
	NewMesh(meshTriangulated, meshFinal, dimension);
    PopulateCell3D(meshFinal);
	
	
	
	// CONFRONTI TRA MESH TRIANGOLATA E ATTESA
	
	// Id vertici
    EXPECT_EQ(meshExpected.Cell0DsId, meshFinal.Cell0DsId);

    // Coordinate vertici
    EXPECT_TRUE(meshExpected.Cell0DsCoordinates.isApprox(meshFinal.Cell0DsCoordinates, 1e-6));

    // Id spigoli
    EXPECT_EQ(meshExpected.Cell1DsId, meshFinal.Cell1DsId);

    // Estremi spigoli
    EXPECT_TRUE(meshExpected.Cell1DsExtrema == meshFinal.Cell1DsExtrema);

    // Flag spigoli
    EXPECT_EQ(meshExpected.Cell1DsFlag, meshFinal.Cell1DsFlag);

    // Id facce
    EXPECT_EQ(meshExpected.Cell2DsId, meshFinal.Cell2DsId);

    // Vertici delle facce
    ASSERT_EQ(meshExpected.Cell2DsVertices.size(), meshFinal.Cell2DsVertices.size());
    for (size_t i = 0; i < meshExpected.Cell2DsVertices.size(); ++i) {
        EXPECT_EQ(meshExpected.Cell2DsVertices[i], meshFinal.Cell2DsVertices[i]) << "Mismatch in Cell2DsVertices at face " << i;
    }

    // Spigoli delle facce
    ASSERT_EQ(meshExpected.Cell2DsEdges.size(), meshFinal.Cell2DsEdges.size());
    for (size_t i = 0; i < meshExpected.Cell2DsEdges.size(); ++i) {
        EXPECT_EQ(meshExpected.Cell2DsEdges[i], meshFinal.Cell2DsEdges[i]) << "Mismatch in Cell2DsEdges at face " << i;
    }

    // POLIEDRI 
    EXPECT_EQ(meshExpected.Cell3DsId, meshFinal.Cell3DsId);
    EXPECT_EQ(meshExpected.Cell3DsVertices, meshFinal.Cell3DsVertices);
    EXPECT_EQ(meshExpected.Cell3DsEdges, meshFinal.Cell3DsEdges);
    EXPECT_EQ(meshExpected.Cell3DsFaces, meshFinal.Cell3DsFaces);
}


TEST(TestPolyedra, TestOrderedEdges)
{ 
   // mesh ottenuta utilizzando la funzione di triangolazione
	PolyhedralMesh meshTriangulated;
	PolyhedralMesh mesh;
	generaTetraedro(mesh);
	
	int q = 3;
	int b = 2;
	int c = 0; 
	unsigned int maxFlag = numeric_limits<unsigned int>::max();
	
	vector<int> dimension =CalcoloVEFPoliedro(q, b, c);
	vector<int> dimensionDuplicated = CalcoloDuplicato(q, b, c, dimension);
	triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
	RimuoviVerticiDuplicati(meshTriangulated);
	RimuoviLatiDuplicati(meshTriangulated);
	
	
	// per ogni faccia i lati sono ordinati in modo che la fine dell'arco e coincida con l'inizio dell'arco successivo (e+1)%E
	// il vertice e della faccia deve corrispondere all'origine dell'arco e

    // ciclo su tutte le facce della mesh triangolata 
    for (size_t f = 0; f < meshTriangulated.Cell2DsId.size(); ++f) {
		const auto& edges = meshTriangulated.Cell2DsEdges[f]; // lista dei lati di una faccia
        const auto& vertices = meshTriangulated.Cell2DsVertices[f]; // lista dei vertici di una faccia 
        size_t E = edges.size(); // numero di vertici della faccia
        
		
		ASSERT_EQ(vertices.size(), E) << "Numero di vertici e di lati non corrispondono per faccia " << f;
		
		// itero su ogni lato e della faccia 
		 for (size_t e = 0; e < E; ++e) {
			unsigned int currentEdge = edges[e]; // lato corrente
            unsigned int nextEdge = edges[(e + 1) % E]; // lato successivo 
            
            unsigned int currentEdgeOrigin = maxFlag;
			unsigned int currentEdgeEnd = maxFlag;
			unsigned int nextEdgeOrigin = maxFlag;
			unsigned int nextEdgeEnd = maxFlag;
            
            bool foundCurrent = false;
            for (unsigned int i = 0; i < meshTriangulated.Cell1DsId.size(); i++) {
                if (currentEdge == meshTriangulated.Cell1DsId[i]) {
                    currentEdgeOrigin = meshTriangulated.Cell1DsExtrema(i, 0);
                    currentEdgeEnd = meshTriangulated.Cell1DsExtrema(i, 1);
                    foundCurrent = true;
                    break; 
                }
            }
           
            ASSERT_TRUE(foundCurrent) << "Master Edge ID " << currentEdge << " not found in Cell1DsId!";

            // trovo gli estremi per nextEdgeMasterId
            bool foundNext = false;
            for (unsigned int i = 0; i < meshTriangulated.Cell1DsId.size(); i++) {
                if (nextEdge == meshTriangulated.Cell1DsId[i]) {
                    nextEdgeOrigin = meshTriangulated.Cell1DsExtrema(i, 0);
                    nextEdgeEnd = meshTriangulated.Cell1DsExtrema(i, 1);
                    foundNext = true;
                    break; 
                }
            }
            ASSERT_TRUE(foundNext) << "Master Edge ID " << nextEdge << " not found in Cell1DsId!";

            // recupero il vertice della faccia
            unsigned int faceVertex = vertices[e]; 

            bool condition = (currentEdgeOrigin == nextEdgeOrigin) ||
                             (currentEdgeOrigin == nextEdgeEnd) ||
                             (currentEdgeEnd == nextEdgeOrigin) ||
                             (currentEdgeEnd == nextEdgeEnd);
            EXPECT_TRUE(condition);

            // Controllo che il vertice e-esimo della faccia coincida con l'origine del lato e-esimo
            
            bool condition2 = (faceVertex == currentEdgeOrigin) ||
            					(faceVertex == currentEdgeEnd);
            EXPECT_TRUE(condition2);

		}		
	}
}

TEST(TestPolyedra, TestNotNullArea)
{
	double eps = numeric_limits<double>::epsilon();

	// mesh ottenuta utilizzando la funzione di triangolazione
	PolyhedralMesh meshTriangulated;
	PolyhedralMesh mesh;
	generaTetraedro(mesh);
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	vector<int> dimension =CalcoloVEFPoliedro(q, b, c);
	vector<int> dimensionDuplicated = CalcoloDuplicato(q, b, c, dimension);
	triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
	RimuoviLatiDuplicati(meshTriangulated);
	RimuoviVerticiDuplicati(meshTriangulated);
	
	// ciclo su tutti i triangoli
	for (size_t i = 0; i < meshTriangulated.Cell2DsVertices.size(); ++i) {
		const auto& tri = meshTriangulated.Cell2DsVertices[i];
		ASSERT_EQ(tri.size(), 3); // controllo se il triangolo ha 3 vertici
		
		// accedo alle coordinate dei vertici
		Vector3d A = meshTriangulated.Cell0DsCoordinates.col(tri[0]); // per ogni vertice prendo la colonna che contiene le coordinate
        Vector3d B = meshTriangulated.Cell0DsCoordinates.col(tri[1]);
        Vector3d C = meshTriangulated.Cell0DsCoordinates.col(tri[2]);
		
		// calcolo l'area del triangolo
		double area = 0.5 * ((B - A).cross(C - A)).norm(); // prodotto vettoriale
		EXPECT_GT(area, eps) << "Triangolo con area nulla o quasi nulla al triangolo " << i;
	}	
	
}

TEST(TestPolyedra, TestNotNullEdges){
	
	double eps = numeric_limits<double>::epsilon();

	// mesh ottenuta utilizzando la funzione di triangolazione
	PolyhedralMesh meshTriangulated;
	PolyhedralMesh mesh;
	generaTetraedro(mesh);
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	vector<int> dimension =CalcoloVEFPoliedro(q, b, c);
	vector<int> dimensionDuplicated = CalcoloDuplicato(q, b, c, dimension);
	triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
	RimuoviLatiDuplicati(meshTriangulated);
	RimuoviVerticiDuplicati(meshTriangulated);
	
	// itero su ogni lato della mesh triangolata 
	for (unsigned int i = 0; i < meshTriangulated.Cell1DsExtrema.rows(); ++i) {
		// prendo gli indici dei due vertici che definisco il lato i
		int vStart = meshTriangulated.Cell1DsExtrema(i, 0); // indice del vertice di partenza 
        int vEnd = meshTriangulated.Cell1DsExtrema(i, 1); // indice del vertice di arrivo
		
		// prendo le coordinate dei due vertici
		Vector3d startPoint = meshTriangulated.Cell0DsCoordinates.col(vStart);
        Vector3d endPoint = meshTriangulated.Cell0DsCoordinates.col(vEnd);
		
		// calcolo la norma della lunghezza del lato
		double length = (endPoint - startPoint).norm();
		
		EXPECT_GT(length, eps) << "Lato con lunghezza nulla o quasi nulla all'edge " << i;
		
	}
}



TEST(TestPolyedra, TestNotNullEdgesDual){
	
	double eps = numeric_limits<double>::epsilon();
	PolyhedralMesh mesh;
	PolyhedralMesh meshDual;
	generaTetraedro(mesh);
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	TriangolazioneDuale(q, b, c, mesh, meshDual);
	
	
	// itero su ogni lato della mesh triangolata 
	for (unsigned int i = 0; i < meshDual.Cell1DsExtrema.rows(); ++i) {
		// prendo gli indici dei due vertici che definisco il lato i
		int vStart = meshDual.Cell1DsExtrema(i, 0); // indice del vertice di partenza 
        int vEnd = meshDual.Cell1DsExtrema(i, 1); // indice del vertice di arrivo
		
		// prendo le coordinate dei due vertici
		Vector3d startPoint = meshDual.Cell0DsCoordinates.col(vStart);
        Vector3d endPoint = meshDual.Cell0DsCoordinates.col(vEnd);
		
		// calcolo la norma della lunghezza del lato
		double length = (endPoint - startPoint).norm();
		
		EXPECT_GT(length, eps) << "Lato con lunghezza nulla o quasi nulla all'edge " << i;
		
	}
}

TEST(TestPolyedra, TestNotNullEdgesTri2){
	
	double eps = numeric_limits<double>::epsilon();
	PolyhedralMesh mesh;
	PolyhedralMesh meshTriangulated;
	generaTetraedro(mesh);
	
	int q = 3;
	int b = 2;
	
	Triangolazione2(q, b, mesh, meshTriangulated);
	
	
	// itero su ogni lato della mesh triangolata 
	for (unsigned int i = 0; i < meshTriangulated.Cell1DsExtrema.rows(); ++i) {
		// prendo gli indici dei due vertici che definisco il lato i
		int vStart = meshTriangulated.Cell1DsExtrema(i, 0); // indice del vertice di partenza 
        int vEnd = meshTriangulated.Cell1DsExtrema(i, 1); // indice del vertice di arrivo
		
		// prendo le coordinate dei due vertici
		Vector3d startPoint = meshTriangulated.Cell0DsCoordinates.col(vStart);
        Vector3d endPoint = meshTriangulated.Cell0DsCoordinates.col(vEnd);
		
		// calcolo la norma della lunghezza del lato
		double length = (endPoint - startPoint).norm();
		
		EXPECT_GT(length, eps) << "Lato con lunghezza nulla o quasi nulla all'edge " << i;
		
	}
}


TEST(TestPolyedra, TestNotNullAreaTri2)
{
	double eps = numeric_limits<double>::epsilon();

	PolyhedralMesh mesh;
	generaTetraedro(mesh);
	
	PolyhedralMesh meshTriangulated;
	PolyhedralMesh meshFinal;
	PolyhedralMesh meshTriangulated2;
	
	int q = 3;
	int b = 2;
	
	vector<int> dimension =CalcoloVEFPoliedro(q, b, 0);
	vector<int> dimensionDuplicated = CalcoloDuplicato(q, b, 0, dimension);
	triangulateAndStore(mesh, meshTriangulated, b, 0,  dimensionDuplicated);
	RimuoviVerticiDuplicati(meshTriangulated);
	RimuoviLatiDuplicati(meshTriangulated);
	NewMesh(meshTriangulated, meshFinal, dimension);
	
	vector<int> dimension2 = CalcoloDimensione2(b, q);
	map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap = buildEdgeToFacesMap(meshFinal);
	triangulateAndStore2(meshFinal, meshTriangulated2, dimension2, edgeToFacesMap);
	
	// ciclo su tutti i triangoli
	for (size_t i = 0; i < meshTriangulated2.Cell2DsVertices.size(); ++i) {
		const auto& tri = meshTriangulated2.Cell2DsVertices[i];
		ASSERT_EQ(tri.size(), 3); // controllo se il triangolo ha 3 vertici
		
		// accedo alle coordinate dei vertici
		Vector3d A = meshTriangulated2.Cell0DsCoordinates.col(tri[0]); // per ogni vertice prendo la colonna che contiene le coordinate
        Vector3d B = meshTriangulated2.Cell0DsCoordinates.col(tri[1]);
        Vector3d C = meshTriangulated2.Cell0DsCoordinates.col(tri[2]);
		
		// calcolo l'area del triangolo
		double area = 0.5 * ((B - A).cross(C - A)).norm(); // prodotto vettoriale
		EXPECT_GT(area, eps) << "Triangolo con area nulla o quasi nulla al triangolo " << i;
	}	
	
}


TEST(Polyedra, DualTest){
	
	// mesh ottenuta utilizzando la funzione di triangolazioneMore actions
	PolyhedralMesh meshTriangulated;

	PolyhedralMesh mesh;
	PolyhedralMesh meshFinal;
	PolyhedralMesh meshDual;
	generaTetraedro(mesh);
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	vector<int> dimension =CalcoloVEFPoliedro(q, b, c);
	vector<int> dimensionDuplicated = CalcoloDuplicato(q, b, c, dimension);
	triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
	RimuoviLatiDuplicati(meshTriangulated);
	RimuoviVerticiDuplicati(meshTriangulated);
	NewMesh(meshTriangulated, meshFinal, dimension);
	map <pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap = buildEdgeToFacesMap(meshFinal);
	
	CalculateDual(meshFinal, meshDual, edgeToFacesMap);
	
	double eps = numeric_limits<double>::epsilon();
	unsigned int maxFlag = numeric_limits<unsigned int>::max();

	size_t expectedVerticesDual   = meshFinal.Cell2DsId.size(); // 16 facce triangolate → 16 vertici nel duale
    size_t expectedEdgesDual      = meshFinal.Cell1DsId.size(); // Calcoloto dinamicamente, può essere 24 per subdivisionLevel=2
   
    // numero di vertici (id)
	// numero di lati(id)
   
    EXPECT_EQ(meshDual.Cell0DsId.size(), expectedVerticesDual);
    EXPECT_EQ(meshDual.Cell1DsId.size(), expectedEdgesDual);
	
	// ogni vertice del poliedro originale genera una faccia nel duale
	EXPECT_GE(meshDual.Cell2DsId.size(), 4);  // Minimo 4 se è un tetraedro chiuso
	
	// verifico se ogni faccia duale ha almento 3 vertici
	for (const auto& dualFace : meshDual.Cell2DsVertices) {
    	EXPECT_GE(dualFace.size(), 3);
    }
	
	// controllo che uno spigolo non connetta un vertice a se stesso
	for (int i = 0; i < meshDual.Cell1DsExtrema.rows(); ++i) {
		int a = meshDual.Cell1DsExtrema(i, 0);
		int b = meshDual.Cell1DsExtrema(i, 1);
		EXPECT_NE(a, b); // Mi verifica che a e b siano diversi 
    }
	
	// verifico che il calcolo del baricentro sia corretto
	// Per ogni faccia della mesh triangolata
    for (size_t faceId = 0; faceId < meshFinal.Cell2DsId.size(); ++faceId) {
        double sumX = 0.0;
        double sumY = 0.0;
        double sumZ = 0.0;
		
        // accedo alla lista dei vertici che compongono la faccia
        const vector<unsigned int>& faceVertices = meshFinal.Cell2DsVertices[faceId];

        // sommo le coordinate dei vertici della faccia
        for (unsigned int v_id : faceVertices) {
			sumX += meshFinal.Cell0DsCoordinates(0, v_id); // coordinata x del vertice v_idMore actions
            sumY += meshFinal.Cell0DsCoordinates(1, v_id); // coordinata y del vertice v_id
            sumZ += meshFinal.Cell0DsCoordinates(2, v_id); // coordinata z del vertice v_id
        }

        // calcolo il baricentro
        double n = (double) faceVertices.size(); // perché faceVertices.size() restituisce un tipo size_t quindi evito una conversione implicita
        double baryX = sumX / n;
        double baryY = sumY / n;
        double baryZ = sumZ / n;

        // estraggo la coordinata del vertice duale corrispondente (baricentro)
        double X = meshDual.Cell0DsCoordinates(0, faceId);
        double Y = meshDual.Cell0DsCoordinates(1, faceId);
        double Z = meshDual.Cell0DsCoordinates(2, faceId);

        EXPECT_NEAR(baryX, X, eps);
        EXPECT_NEAR(baryY, Y, eps);
        EXPECT_NEAR(baryZ, Z, eps);
    }

    // verifico che i lati siano costruiti nel modo corretto ---> devono connettere i baricentri
	for (int i = 0; i < meshDual.Cell1DsExtrema.rows(); ++i) {
		Vector2i edge = meshDual.Cell1DsExtrema.row(i);
		// id delle facce del poliedro originale (vertici duali sono baricentri di facce originali)
        int f1 = edge(0); 
        int f2 = edge(1);
		 // controlla che le due facce originali condividano un edge
    bool foundCommonEdge = false;
    for (int e1 : meshFinal.Cell2DsEdges[f1]) {
        for (int e2 : meshFinal.Cell2DsEdges[f2]) {
            if (e1 == e2) {
                foundCommonEdge = true;
                break;
            }
        }
        if (foundCommonEdge) break;
    }
    EXPECT_TRUE(foundCommonEdge);
	}
	
	// verifico che i lati e le facce siano ordinate correttamente 
	
	for (size_t f = 0; f < meshDual.Cell2DsId.size(); ++f) {
		const auto& edges = meshDual.Cell2DsEdges[f]; // lista dei lati di una faccia
        const auto& vertices = meshDual.Cell2DsVertices[f]; // lista dei vertici di una faccia 
        size_t E = edges.size(); // numero di vertici dela faccia
       
		// itero su ogni lato e della faccia 
		 for (size_t e = 0; e < E; ++e) {
			unsigned int currentEdge = edges[e]; // lato corrente
            unsigned int nextEdge = edges[(e + 1) % E]; // lato successivo 
            
            unsigned int currentEdgeOrigin = maxFlag;
			unsigned int currentEdgeEnd = maxFlag;
			unsigned int nextEdgeOrigin = maxFlag;
			unsigned int nextEdgeEnd = maxFlag;
            
            bool foundCurrent = false;
            for (unsigned int i = 0; i < meshDual.Cell1DsId.size(); i++) {
                if (currentEdge == meshDual.Cell1DsId[i]) {
                    currentEdgeOrigin = meshDual.Cell1DsExtrema(i, 0);
                    currentEdgeEnd = meshDual.Cell1DsExtrema(i, 1);
                    foundCurrent = true;
                    break; 
                }
            }
            
            ASSERT_TRUE(foundCurrent) << "Master Edge ID " << currentEdge << " not found in Cell1DsId!";

            // trovo gli estremi per nextEdgeMasterId
            bool foundNext = false;
            for (unsigned int i = 0; i < meshDual.Cell1DsId.size(); i++) {
                if (nextEdge == meshDual.Cell1DsId[i]) {
                    nextEdgeOrigin = meshDual.Cell1DsExtrema(i, 0);
                    nextEdgeEnd = meshDual.Cell1DsExtrema(i, 1);
                    foundNext = true;
                    break; 
                }
            }
            ASSERT_TRUE(foundNext);

            // Recupera il vertice della faccia
            unsigned int faceVertex = vertices[e]; // ID master del vertice della faccia

            bool condition = (currentEdgeOrigin == nextEdgeOrigin) ||
                             (currentEdgeOrigin == nextEdgeEnd) ||
                             (currentEdgeEnd == nextEdgeOrigin) ||
                             (currentEdgeEnd == nextEdgeEnd);
            EXPECT_TRUE(condition);
              
            // Controllo che il vertice e-esimo della faccia coincida con l'origine del lato e-esimo
            bool condition2 = (faceVertex == currentEdgeOrigin) ||
            					(faceVertex == currentEdgeEnd);
            EXPECT_TRUE(condition2);
		}		
	}	 
}


//TRIANGULATION2

TEST(TestPolyedra, TestOrderedEdges2)
{ 
	PolyhedralLibrary::PolyhedralMesh mesh;
	PolyhedralLibrary::generaTetraedro(mesh);
	
	PolyhedralMesh meshTriangulated;
	PolyhedralMesh meshFinal;
	PolyhedralMesh meshTriangulated2;
	
	int q = 3;
	int b = 2;
	unsigned int maxFlag = numeric_limits<unsigned int>::max();
	
	vector<int> dimension =CalcoloVEFPoliedro(q, b, 0);
	vector<int> dimensionDuplicated = CalcoloDuplicato(q, b, 0, dimension);
	triangulateAndStore(mesh, meshTriangulated, b, 0,  dimensionDuplicated);
	RimuoviVerticiDuplicati(meshTriangulated);
	RimuoviLatiDuplicati(meshTriangulated);
	NewMesh(meshTriangulated, meshFinal, dimension);
	
	vector<int> dimension2 = CalcoloDimensione2(b, q);
	map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap = buildEdgeToFacesMap(meshFinal);
	triangulateAndStore2(meshFinal, meshTriangulated2, dimension2, edgeToFacesMap);
	
	// per ogni faccia i lati sono ordinati in modo che la fine dell'arco e coincida con l'inizio dell'arco successivo (e+1)%E
	// il vertice e della faccia deve corrispondere all'origine dell'arco e

    // ciclo su tutte le facce della mesh triangolata 
    for (size_t f = 0; f < meshFinal.Cell2DsId.size(); ++f) {
		const auto& edges = meshFinal.Cell2DsEdges[f]; // lista dei lati di una faccia
        const auto& vertices = meshFinal.Cell2DsVertices[f]; // lista dei vertici di una faccia 
        size_t E = edges.size(); // numero di vertici dela faccia
        
		
		ASSERT_EQ(vertices.size(), E) << "Numero di vertici e di lati non corrispondono per faccia " << f;
		
		// itero su ogni lato e della faccia 
		 for (size_t e = 0; e < E; ++e) {
			unsigned int currentEdge = edges[e]; // lato corrente
            unsigned int nextEdge = edges[(e + 1) % E]; // lato successivo 
            
            unsigned int currentEdgeOrigin = maxFlag;
			unsigned int currentEdgeEnd = maxFlag;
			unsigned int nextEdgeOrigin = maxFlag;
			unsigned int nextEdgeEnd = maxFlag;
            
            bool foundCurrent = false;
            for (unsigned int i = 0; i < meshFinal.Cell1DsId.size(); i++) {
                if (currentEdge == meshFinal.Cell1DsId[i]) {
                    currentEdgeOrigin = meshFinal.Cell1DsExtrema(i, 0);
                    currentEdgeEnd = meshFinal.Cell1DsExtrema(i, 1);
                    foundCurrent = true;
                    break; 
                }
            }
           
            ASSERT_TRUE(foundCurrent) << "Master Edge ID " << currentEdge << " not found in Cell1DsId!";

            // trovo gli estremi per nextEdgeMasterId
            bool foundNext = false;
            for (unsigned int i = 0; i < meshFinal.Cell1DsId.size(); i++) {
                if (nextEdge == meshFinal.Cell1DsId[i]) {
                    nextEdgeOrigin = meshFinal.Cell1DsExtrema(i, 0);
                    nextEdgeEnd = meshFinal.Cell1DsExtrema(i, 1);
                    foundNext = true;
                    break; 
                }
            }
            ASSERT_TRUE(foundNext) << "Master Edge ID " << nextEdge << " not found in Cell1DsId!";

            // recupero il vertice della faccia
            unsigned int faceVertex = vertices[e]; 

            bool condition = (currentEdgeOrigin == nextEdgeOrigin) ||
                             (currentEdgeOrigin == nextEdgeEnd) ||
                             (currentEdgeEnd == nextEdgeOrigin) ||
                             (currentEdgeEnd == nextEdgeEnd);
            EXPECT_TRUE(condition);

            // Controllo che il vertice e-esimo della faccia coincida con l'origine del lato e-esimo
            
            bool condition2 = (faceVertex == currentEdgeOrigin) ||
            					(faceVertex == currentEdgeEnd);
            EXPECT_TRUE(condition2);

		}		
	}
}

// CAMMINO MINIMO
TEST(TestPolyedra, ShortestPath){
	
	PolyhedralMesh mesh;
	PolyhedralMesh meshTriangulated;
	PolyhedralMesh meshFinal;
	generaTetraedro(mesh);
	
	int q = 3;
	int b = 2;
	int c = 0;
	
	int startVertexId = 0;
	int endVertexId = 7;	
	
	vector<int> dimension =CalcoloVEFPoliedro(q, b, c);
	vector<int> dimensionDuplicated = CalcoloDuplicato(q, b, c, dimension);
	triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
	RimuoviLatiDuplicati(meshTriangulated);
	RimuoviVerticiDuplicati(meshTriangulated);
	NewMesh(meshTriangulated, meshFinal, dimension);
	
	ShortestPathResult expected(2, 1.63299);
	MatrixXi adjMatrix = calculateAdjacencyMatrix(meshFinal);
	ShortestPathResult result = findShortestPathDijkstra(meshFinal, adjMatrix, startVertexId, endVertexId);
	
	EXPECT_NEAR(result.numEdges, expected.numEdges, 1e-5) 
		<< "La lunghezza del percorso più breve non corrisponde all'atteso."
		<< " Atteso: " << expected.numEdges << ", Trovato: " << result.numEdges;
		
	EXPECT_NEAR(result.totalLength, expected.totalLength, 1e-5) 
        << "La lunghezza del percorso più breve non corrisponde all'atteso."
        << " Atteso: " << expected.totalLength << ", Trovato: " << result.totalLength;
	 
	}

}