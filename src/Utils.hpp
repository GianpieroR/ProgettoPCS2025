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
	vector<int> CalcoloVEFPoliedro(const int q, const int b, const int c);
	
	// Calcola il numero di vertici, lati e facce del poliedro triangolato (considerati i duplicati)
	// dove q, b, c sono i parametri passati dall'utente che identificano il poliedro e la sua triangolazione
	//  e dimension è il vettore che contiene il numero di vertici, lati e facce del poliedro triangolato (senza duplicati)
	vector<int> CalcoloDuplicato(const int q, const int b, const int c, const vector<int>& dimension);
	
	// Calcola il numero di vertici, lati e facce del poliedro triangolato con la triangolazione II
	// dove q, b  sono i parametri passati dall'utente che identificano il poliedro e la sua triangolazione
	vector<int> CalcoloDimensione2(const int b, const int q);
	
	// Assegna un flag ai vertici che devono essere tenuti (non duplicati)
	// meshTriangulated : una struct MeshPoliedrale triangolata
	void RimuoviVerticiDuplicati(PolyhedralMesh& meshTriangulated);
	
	// Assegna un flag ai lati che devono essere tenuti (non duplicati)
	// meshTriangulated : una struct PolyhedralMesh triangolata
	void RimuoviLatiDuplicati(PolyhedralMesh& meshTriangulated);
	
	// Crea una nuova mesh senza duplicati dove
	// meshTriangulated è una struct MeshPoliedrale (con i duplicati)
	// meshFinal è una struct MeshPoliedrale (senza i duplicati)
	// dimension è un vettore che contiene il numero di vertici, lati e facce del poliedro triangolato (senza duplicati)
	void NewMesh(const PolyhedralMesh& meshTriangulated, PolyhedralMesh& meshFinal, const vector<int>& dimension);
	
	// ----------------------------------------------------------------------------------------------------------- // 
	
	// POLYHEDRA 
	
	// Riempiono la struct MeshPoliedrale con i dati dei poliedri non triangolati
	// mesh: una struct PolyhedralMesh
	void generaTetraedro(PolyhedralMesh& mesh);
	void generaOttaedro(PolyhedralMesh& mesh);
	void generaIcosaedro(PolyhedralMesh& mesh);
	
	// ----------------------------------------------------------------------------------------------------------- //
	
	// TRIANGOLAZIONE
	
	// Racchiude le funzioni per la triangolazione di classe I se p=3 dove
	// q, b, c sono i parametri passati dall'utente che identificano il poliedro e la sua triangolazione
	// mesh è una struct PolyhedralMesh
	void Triangolazione(const int q, const int b, const int c, PolyhedralMesh& mesh, PolyhedralMesh& meshFinal);
	
	// Racchiude le funzioni per la triangolazione di classe I se q=3 dove
	// q, b, c sono i parametri passati dall'utente che identificano il poliedro e la sua triangolazione
	// mesh è una una struct PolyhedralMesh
	void TriangolazioneDuale(const int q, const int b, const int c, PolyhedralMesh& mesh,  PolyhedralMesh& meshFinal);
	
	// Racchiude le funzioni per la triangolazione di classe II se p=3
	// q, b sono i parametri passati dall'utente che identificano il poliedro e la sua triangolazione
	// mesh: una struct PolyhedralMesh
	void Triangolazione2(const int q, const int b, PolyhedralMesh& mesh,  PolyhedralMesh& meshFinal);
	
	// Riempie le Celle3d dopo la triangolazione
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	void PopulateCell3D(PolyhedralMesh& meshTriangulated);
	
	// ----------------------------------------------------------------------------------------------------------- //
	
	// TRIANGOLAZIONE1
	
	// Triangola il poliedro di classe I
	// mesh : una struct MeshPoliedrale, quella di partenza non triangolata
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	// b,c sono i parametri passati dall'utente che identificano il poliedro
	// dimension è il vettore che contiene il numero di vertici, lati e facce del poliedro triangolato
	void triangulateAndStore(PolyhedralMesh& mesh, PolyhedralMesh& meshTriangulated,
							  const unsigned int b, const unsigned int c, const vector<int>& dimension);
	
	 // Aggiunge un lato alla mesh triangolata se non è già presente
	 // a : id del primo estremo (vertice) del lato
	 // b : id del secondo estremo (vertice) del lato
	 // meshTriangulated : una struct PolyhedralMesh, quella triangolata
	 // edgeId : id del nuovo lato
	 // triangleId : id della faccia a cui appartiene
	void FindAddEdge(const unsigned int a, const unsigned int b, PolyhedralMesh& meshTriangulated,
					 unsigned int& edgeID, const unsigned int triangleID);
	
	// ----------------------------------------------------------------------------------------------------------- //
	
	// TRIANGOLAZIONE2
					 
	// Triangola il poliedro di classe II
	 // mesh è una struct MeshPoliedrale, quella di partenza non triangolata
	 // meshTriangulated è una struct PolyhedralMesh, quella triangolata
	 // b,c sono i parametri passati dall'utente che identificano il poliedro
	 // dimension è un vettore che contiene il numero di vertici, lati e facce del poliedro triangolato
	 // edgeToFacesMap è mappa che associa ad ogni spigolo del poliedro originale l'elenco di tutte le facce che contengono quello spigolo
	void triangulateAndStore2(PolyhedralMesh& mesh, PolyhedralMesh& meshTriangulated, const vector<int>& dimension, map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap);
							  
	// Aggiunge un vertice alla mesh triangolata se non è già presente e restituisce il suo id
	// coord sono le coordinate del punto che vorremmo aggiungere
	// meshTriangulated è una struct PolyhedralMesh, quella triangolata
	// k1 è prossimo id per i vertici disponibile 
	unsigned int FindAddVertice(const Vector3d& coord, PolyhedralMesh& meshTriangulated, unsigned int& k1);
	
	// Aggiunge un lato alla mesh triangolata se non è già presente dove
	// a è id del primo estremo (vertice) del lato
	// b è id del secondo estremo (vertice) del lato
	// meshTriangulated è una struct PolyhedralMesh, quella triangolata
	// k2 è prossimo id per i lati disponibile
	unsigned int FindAddEdge2(const unsigned int a, const unsigned int b, PolyhedralMesh& meshTriangulated, unsigned int& k2);
	
	// A partire da un lato e una faccia che stiamo considerando, trova l'altra faccia adiacente
	// al lato e restituisce il suo baricentro
	// meshTriangulated è una struct MeshPoliedrale
	// edgeId è id del lato
	// currentFacdeId è faccia che stiamo considerando
	// edgeToFacesMap è mappa che associa ad ogni spigolo del poliedro originale l'elenco di tutte le facce che contengono quello spigolo
	Vector3d FindNearBarycenter(const PolyhedralMesh& meshTriangulated, const unsigned int edgeId, const unsigned int currentFaceId, map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap);
	
	// Aggiunge una faccia alla mesh triangolata se non è già presente
	// new_face_vertices : id dei vertici della faccia che vorremmo aggiungere
	// new_face_edges : id dei lati della faccia che vorremmo aggiungere
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	// k3 : prossimo id per le facce disponibile
	void FindAddFace(const vector<unsigned int>& new_face_vertices, const vector<unsigned int>& new_face_edges, PolyhedralMesh& meshTriangulated, unsigned int& k3);
	
	// Ruota la sequenza degli id dei lati in modo che inizi con il più piccolo (8,5,10) --> (5,10,8)
	// current_edges : id dei lati che vogliamo riordinare
	vector<unsigned int> get_cyclic_normalized(const vector<unsigned int>& current_edges);
	
	// Normalizza l'ordine combinando la normalizzazione ciclica con un confronto delle forme invertite (0,1,2) e (0,2,1) sono uguali
	// face_edges : id dei lati che vogliamo riordinare
	vector<unsigned int> NormalizeFaceEdges(const vector<unsigned int>& face_edges);
	
	// ----------------------------------------------------------------------------------------------------------- //
	
	// DUALE
	
	// Calcola il duale di un poliedro
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	// meshDual : una struct PolyhedralMesh, quella duale
	// edgeToFacesMap : mappa che associa ad ogni spigolo del poliedro originale l'elenco di tutte le facce che contengono quello spigolo
	void CalculateDual(PolyhedralMesh& meshTriangulated, PolyhedralMesh& meshDual, map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap);
	
	// Calcola il baricentro di una faccia
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	// faceId : id della faccia di cui vogliamo calcolare il baricentro
	Vector3d getFaceBarycenter(const PolyhedralMesh& meshTriangulated, const unsigned int faceId);
	
	// Crea una mappa che associa ad ogni spigolo del poliedro originale l'elenco di tutte le facce che contengono quello spigolo
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
        map <pair<unsigned int, unsigned int>, vector<unsigned int>> buildEdgeToFacesMap(const PolyhedralMesh& meshTriangulated);
	
	// Crea una mappa che associa ad ogni vertice originale l'elenco di tutte le facce che contengono quel vertice
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	map<unsigned int, vector<unsigned int>> buildVertexToFacesMap(const PolyhedralMesh& meshTriangulated);
	
	// Crea una mappa che associa ad ogni vertice originale l'elenco di tutti gli spigoli che vi incidono
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	map<unsigned int, vector<unsigned int>> buildVertexToEdgesMap(const PolyhedralMesh& meshTriangulated);
	
	// Proietta i vertici del poliedro sulla sfera unitaria
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	void ProjectMeshToUnitSphere(PolyhedralMesh& meshTriangulated);
	
	// ----------------------------------------------------------------------------------------------------------- //
	
	// CAMMINO MINIMO
	 
	struct ShortestPathResult {
		unsigned int numEdges; // numero di lati del cammino
		double totalLength; // lunghezza totale del cammino
		
		// Costruttore con parametri per la dimensione dei vettori.
		ShortestPathResult(unsigned int nEdges, double j);
	};
	
	// Calcola la distanza euclidea tra due punti
	// mesh : una struct PolyhedralMesh
	// id1, id2 : id dei vertici di cui vogliamo calcolare la distanza 
	double calculateDistanceById(const PolyhedralMesh& mesh, const unsigned int id1, const unsigned int id2);
	
	// Calcola la matrice di adiancenza
	// djacencyMatrix[i,j] = 1 se esiste un lato tra il vertice i e il vertice j
	// mesh : una struct PolyhedralMesh
	MatrixXi calculateAdjacencyMatrix(const PolyhedralMesh& mesh);
	
	// Calcola il cammino minimo con l'algoritmo Dijkstra
	// Resituisce un oggetto ShortestPathResult che contiene le informazioni sul cammino minimo trovato
	// mesh : una struct PolyhedralMesh
	// adjMatrix : matrice di adiacenza che indica la connettività tra i vertici
	// startVertexId_real : id del vertice di partenza
	// endVertexId_real : id del vertice di arrivo
	ShortestPathResult findShortestPathDijkstra(
	    PolyhedralMesh& mesh,
	    const MatrixXi& adjMatrix,
	    const unsigned int startVertexId_real,
	    const unsigned int endVertexId_real
	);
	
	
	// ----------------------------------------------------------------------------------------------------------- //
	
	// EXPORT PARAVIEW
	
	// Esporta la mesh triangolata
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	void ExportParaview(const PolyhedralMesh& meshTriangulated);
	
	// Stampa la mesh triangolata
	// meshTriangulated : una struct PolyhedralMesh, quella triangolata
	void printMeshTriangulated(const PolyhedralMesh& meshTriangulated);
	
	// ----------------------------------------------------------------------------------------------------------- //
	
	// MESH EXPORT

	// Scrivono sui file TXT
	// mesh: una struct PolyhedralMesh
	void WriteCell0Ds(const PolyhedralMesh& mesh);
        void WriteCell1Ds(const PolyhedralMesh& mesh);
        void WriteCell2Ds(const PolyhedralMesh& mesh);
        void WriteCell3Ds(const PolyhedralMesh& mesh);
    
       // ----------------------------------------------------------------------------------------------------------- //

}

	