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

        // TRIANGULATION
        // Le seguenti funzioni gestiscono la triangolazione di poliedri secondo diverse classi (Classe I e Classe II).
        // I parametri `q`, `b`, `c` determinano la struttura del poliedro e il tipo di triangolazione da eseguire.
        // Tutte le funzioni operano su oggetti di tipo PolyhedralMesh, aggiornando la mesh triangolata.

        // Triangola un poliedro di Classe I nel caso in cui p = 3.
  
        // Questa funzione costruisce la triangolazione di un poliedro regolare o semiregolare sulla base
        // dei parametri forniti. È utilizzata per ottenere una mesh finale composta da soli triangoli.
 
        @param q Parametro che rappresenta il numero di lati delle facce del poliedro.
        @param b Parametro strutturale che influenza la forma della mesh.
        @param c Parametro aggiuntivo per la configurazione geometrica.
        @param mesh Mesh originale di input (non triangolata).
        @param meshFinal Mesh finale triangolata in uscita.

        void Triangulation(int q, int b, int c, PolyhedralMesh& mesh, PolyhedralMesh& meshFinal);


        // Calcola il duale della triangolazione di Classe I nel caso in cui q = 3.
  
        // Dopo aver triangolato il poliedro iniziale, la funzione costruisce la mesh duale associata.

 
        @param q Parametro che rappresenta il numero di lati delle facce (qui fisso a 3).
        @param b Parametro strutturale del poliedro.
        @param c Parametro configurativo.
        @param mesh Mesh originale di input.

        void TriangulationDual(int q, int b, int c, PolyhedralMesh& mesh, PolyhedralMesh& meshFinal);


        // Triangola un poliedro di Classe II nel caso in cui p = 3.

        // Questa funzione è analoga a Triangulation, ma utilizzata nel contesto della Classe II.
        // Utilizza una diversa logica strutturale nella generazione della mesh triangolata.
 
        @param q Parametro che identifica la forma del poliedro.
        @param b Parametro che ne controlla la struttura.
        @param mesh Mesh originale di input.
        @param meshFinal Mesh triangolata in uscita.

        void Triangulation2(int q, int b, PolyhedralMesh& mesh, PolyhedralMesh& meshFinal);


        // Calcola il duale della triangolazione di Classe II nel caso in cui q = 3.
 
        // Dopo la triangolazione secondo la logica della Classe II, viene calcolato il duale
        // della mesh ottenuta. Restituisce una rappresentazione alternativa utile per analisi topologiche.
 
        @param q Parametro geometrico del poliedro (fissato a 3).
        @param b Parametro che influisce sulla struttura.
        @param mesh Mesh di input (non triangolata).
        @param meshFinal Mesh duale risultante.
 
        void Triangulation2Dual(int q, int b, PolyhedralMesh& mesh, PolyhedralMesh& meshFinal);


        // Popola la struttura delle celle 3D a partire dalla mesh triangolata.
 
        // Questa funzione aggiorna la mesh specificando le relazioni tra celle 0D (vertici),
        // celle 1D (spigoli), celle 2D (facce) e celle 3D (volumi). 
  
        @param meshTriangulated Mesh già triangolata.
        @param dimension Vettore contenente il numero di vertici, spigoli e facce della mesh.
 
        void PopulateCell3D(PolyhedralMesh& meshTriangulated, const vector<int>& dimension);

/ ----------------------------------------------------------------------------------------------------------- //

// TRIANGOLAZIONE POLIEDRI DI CLASSE I

// Triangola un poliedro di classe I
// mesh : poliedro di partenza non triangolato (struct PolyhedralMesh)
// meshTriangulated : poliedro triangolato (struct PolyhedralMesh)
// b, c : parametri utente che identificano il poliedro
// dimension : vettore contenente il numero di vertici, lati e facce del poliedro triangolato
void triangulateAndStore(PolyhedralMesh& mesh, PolyhedralMesh& meshTriangulated,
                         const unsigned int b, const unsigned int c, const vector<int>& dimension);

// Aggiunge un lato alla mesh triangolata, se non è già presente
// a : id del primo vertice del lato
// b : id del secondo vertice del lato
// meshTriangulated : mesh triangolata (struct PolyhedralMesh)
// edgeID : id assegnato al nuovo lato (output)
// triangleID : id della faccia a cui il lato appartiene
void FindAddEdge(const unsigned int a, const unsigned int b, PolyhedralMesh& meshTriangulated,
                 unsigned int& edgeID, const unsigned int triangleID);

// ----------------------------------------------------------------------------------------------------------- //

// TRIANGOLAZIONE POLIEDRI DI CLASSE II

// Triangola un poliedro di classe II
// mesh : poliedro di partenza non triangolato (struct PolyhedralMesh)
// meshTriangulated : poliedro triangolato (struct PolyhedralMesh)
// dimension : vettore contenente numero di vertici, lati e facce triangolati
// edgeToFacesMap : mappa che associa ad ogni spigolo del poliedro originale l'elenco delle facce che lo contengono
void triangulateAndStore2(PolyhedralMesh& mesh, PolyhedralMesh& meshTriangulated,
                          const vector<int>& dimension, map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap);

// Aggiunge un vertice alla mesh triangolata se non è già presente e restituisce il suo id
// coord : coordinate del vertice da aggiungere
// meshTriangulated : mesh triangolata (struct PolyhedralMesh)
// k1 : prossimo id disponibile per i vertici (modificato in uscita)
unsigned int FindAddVertice(const Vector3d& coord, PolyhedralMesh& meshTriangulated, unsigned int& k1);

// Aggiunge un lato alla mesh triangolata se non è già presente e restituisce il suo id
// a, b : id dei vertici estremi del lato
// meshTriangulated : mesh triangolata (struct PolyhedralMesh)
// k2 : prossimo id disponibile per i lati (modificato in uscita)
unsigned int FindAddEdge2(const unsigned int a, const unsigned int b, PolyhedralMesh& meshTriangulated, unsigned int& k2);

// Dato un lato e una faccia, trova la faccia adiacente lungo quel lato e ne restituisce il baricentro
// meshTriangulated : mesh triangolata (struct PolyhedralMesh)
// edgeId : id del lato considerato
// currentFaceId : id della faccia di riferimento
// edgeToFacesMap : mappa spigoli → facce originali
Vector3d FindNearBarycenter(const PolyhedralMesh& meshTriangulated, const unsigned int edgeId,
                           const unsigned int currentFaceId, map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap);

// Aggiunge una faccia alla mesh triangolata se non è già presente
// new_face_vertices : id dei vertici della nuova faccia
// new_face_edges : id dei lati della nuova faccia
// meshTriangulated : mesh triangolata (struct PolyhedralMesh)
// k3 : prossimo id disponibile per le facce (modificato in uscita)
void FindAddFace(const vector<unsigned int>& new_face_vertices, const vector<unsigned int>& new_face_edges,
                 PolyhedralMesh& meshTriangulated, unsigned int& k3);

// Riordina ciclicamente una sequenza di id di lati in modo che inizi con il più piccolo
// current_edges : sequenza di id dei lati da normalizzare
vector<unsigned int> get_cyclic_normalized(const vector<unsigned int>& current_edges);

// Normalizza l’ordine degli id dei lati combinando normalizzazione ciclica e confronto tra ordine diretto e invertito
// face_edges : sequenza di id dei lati da normalizzare
vector<unsigned int> NormalizeFaceEdges(const vector<unsigned int>& face_edges);

// ----------------------------------------------------------------------------------------------------------- //

// DUALITA' DEL POLIEDRO

// Calcola il poliedro duale a partire da quello triangolato
// meshTriangulated : poliedro triangolato (struct PolyhedralMesh)
// meshDual : poliedro duale risultante (struct PolyhedralMesh)
// edgeToFacesMap : mappa spigoli → facce originali
void CalculateDual(PolyhedralMesh& meshTriangulated, PolyhedralMesh& meshDual,
                   map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap);

// Calcola il baricentro di una faccia specifica
// meshTriangulated : poliedro triangolato (struct PolyhedralMesh)
// faceId : id della faccia di cui calcolare il baricentro
Vector3d getFaceBarycenter(const PolyhedralMesh& meshTriangulated, const unsigned int faceId);

// Costruisce la mappa spigoli → facce per il poliedro triangolato
// meshTriangulated : poliedro triangolato (struct PolyhedralMesh)
map<pair<unsigned int, unsigned int>, vector<unsigned int>> buildEdgeToFacesMap(const PolyhedralMesh& meshTriangulated);

// Costruisce la mappa vertici → facce per il poliedro triangolato
// meshTriangulated : poliedro triangolato (struct PolyhedralMesh)
map<unsigned int, vector<unsigned int>> buildVertexToFacesMap(const PolyhedralMesh& meshTriangulated);

// Costruisce la mappa vertici → lati per il poliedro triangolato
// meshTriangulated : poliedro triangolato (struct PolyhedralMesh)
map<unsigned int, vector<unsigned int>> buildVertexToEdgesMap(const PolyhedralMesh& meshTriangulated);

// Proietta i vertici del poliedro sulla sfera unitaria
// meshTriangulated : poliedro triangolato (struct PolyhedralMesh)
void ProjectMeshToUnitSphere(PolyhedralMesh& meshTriangulated);

// ----------------------------------------------------------------------------------------------------------- //

// CAMMINO MINIMO

// Struttura per il risultato del cammino minimo
struct ShortestPathResult {
    unsigned int numEdges;   // Numero di lati nel cammino
    double totalLength;      // Lunghezza totale del cammino
    
    // Costruttore con parametri per dimensione e lunghezza
    ShortestPathResult(unsigned int nEdges, double length);
};

// Calcola la distanza euclidea tra due vertici identificati dai loro id
// mesh : poliedro (struct PolyhedralMesh)
// id1, id2 : id dei vertici
double calculateDistanceById(const PolyhedralMesh& mesh, const unsigned int id1, const unsigned int id2);

// Calcola la matrice di adiacenza del poliedro
// adjacencyMatrix[i,j] = 1 se esiste un lato tra i vertici i e j, 0 altrimenti
// mesh : poliedro (struct PolyhedralMesh)
MatrixXi calculateAdjacencyMatrix(const PolyhedralMesh& mesh);

// Calcola il cammino minimo tra due vertici usando l'algoritmo di Dijkstra
// Restituisce un ShortestPathResult con dettagli sul cammino trovato
// mesh : poliedro (struct PolyhedralMesh)
// adjMatrix : matrice di adiacenza tra i vertici
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

	