#pragma once
#include <iostream>
#include "PolyhedralMesh.hpp"
#include <Eigen/Dense>


using namespace PolyhedralLibrary;
using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{
	// DIMENSION
	
	// Inverte i valori di p e q
	void invertiValori(int& p, int& q);
	
	// Calcola il numero di vertici, lati e facce del poliedro triangolato con la triangolazione I
	// q, b, c : parametri passati dall'utente che identificano il poliedro e la sua triangolazione
	vector<int> ComputePolyhedronVEF(int q, int b, int c);
	
	// Calcola il numero di vertici, lati e facce del poliedro triangolato (considerati i duplicati)
	// q, b, c : parametri passati dall'utente che identificano il poliedro e la sua triangolazione
	// dimension : vettore che contiene il numero di vertici, lati e facce del poliedro triangolato (senza duplicati)
	vector<int> CalculateDuplicated(int q, int b, int c, const vector<int>& dimension);
	
	// Calcola il numero di vertici, lati e facce del poliedro triangolato con la triangolazione II
	// q, b : parametri passati dall'utente che identificano il poliedro e la sua triangolazione
	vector<int> CalculateDimension2(int b, int q);
	
	// Assegna un flag ai vertici che devono essere tenuti (non duplicati)
	// meshTriangulated : una struct PolyhedralMesh triangolata
	void RemoveDuplicatedVertices(PolyhedralMesh& meshTriangulated);
	
	// Assegna un flag ai lati che devono essere tenuti (non duplicati)
	// meshTriangulated : una struct PolyhedralMesh triangolata
	void RemoveDuplicatedEdges(PolyhedralMesh& meshTriangulated);
	
	// Crea una nuova mesh senza duplicati
	// meshTriangulated : una struct PolyhedralMesh (con i duplicati)
	// meshFinal : una struct PolyhedralMesh (senza i duplicati)
	// dimension : vettore che contiene il numero di vertici, lati e facce del poliedro triangolato (senza duplicati)
	void NewMesh(PolyhedralMesh& meshTriangulated, PolyhedralMesh& meshFinal, const vector<int>& dimension);
	
	// ----------------------------------------------------------------------------------------------------------- // 
	
	// POLYHEDRA 
	
	// Riempiono la struct PolyhedralMesh con i dati dei poliedri non triangolati
	// mesh: una struct PolyhedralMesh
	void generateTetrahedron(PolyhedralMesh& mesh);
	void generateCube(PolyhedralMesh& mesh);
	void generateOctahedron(PolyhedralMesh& mesh);
	void generateDodecahedron(PolyhedralMesh& mesh);
	void generateIcosahedron(PolyhedralMesh& mesh);
	
	// ----------------------------------------------------------------------------------------------------------- //
	