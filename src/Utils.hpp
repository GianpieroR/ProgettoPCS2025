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

	