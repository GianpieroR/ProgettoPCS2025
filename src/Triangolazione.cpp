#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <limits>

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{
	// La funzione Triangolazione esegue la generazione di una mesh triangolata a partire da una mesh poliedrica iniziale,
	// applicando un certo livello di suddivisione controllato dai parametri `b` e `c`, su un poliedro identificato da `q`.
	//
	// Il processo include:
	// 1. Calcolo teorico del numero di vertici, spigoli e facce previste (`dimension`).
	// 2. Calcolo del numero "duplicato" di vertici/spigoli necessario per la triangolazione iniziale (`dimensionDuplicated`),
	//    che considera la possibile presenza di duplicati prima della pulizia della mesh.
	// 3. Esecuzione della triangolazione e memorizzazione temporanea nella mesh `meshTriangulated`.
	// 4. Rimozione dei vertici e spigoli duplicati, per evitare ridondanze geometriche e topologiche.
	// 5. Creazione della mesh finale (`meshFinal`) a partire da quella triangolata e pulita.
	// 6. Costruzione della connettività delle celle 3D (cellule volumetriche) nella mesh finale.
	//
	// Il risultato è una mesh triangolata, strutturalmente coerente, pronta per l’uso in analisi numeriche o operazioni successive.
	void Triangolazione(const int q, const int b, const int c, PolyhedralMesh& mesh, PolyhedralMesh& meshFinal)
	{
		PolyhedralMesh meshTriangulated;

		vector<int> dimension = CalcoloVEFPoliedro(q, b, c);
		vector<int> dimensionDuplicated = CalcoloDuplicato(q, b, c, dimension);
		triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
		RimuoviVerticiDuplicati(meshTriangulated);
    	RimuoviLatiDuplicati(meshTriangulated);
    	NewMesh(meshTriangulated, meshFinal, dimension);
    	PopulateCell3D(meshFinal);
    }
    
	// La funzione TriangolazioneDuale esegue la costruzione della *dual mesh* (mesh duale)
	// a partire da una mesh poliedrica iniziale, tramite triangolazione e trasformazione duale.
	//
	// Il processo completo include:
	// 1. Calcolo delle dimensioni previste per la mesh (con e senza duplicati).
	// 2. Triangolazione della mesh di partenza e salvataggio dei dati nella mesh temporanea `meshTriangulated`.
	// 3. Rimozione di vertici e spigoli duplicati per ripulire la mesh triangolata.
	// 4. Costruzione della mesh finale triangolata `meshFinal` (come nella funzione Triangolazione).
	// 5. Costruzione della mappa spigolo → facce (edgeToFacesMap), che permette di individuare la topologia necessaria per la dualizzazione.
	// 6. Calcolo della mesh duale `meshDual` a partire dalla mesh triangolata e la mappa degli spigoli.
	// 7. Costruzione delle celle 3D della mesh duale per completare la struttura topologica.
	//

    void TriangolazioneDuale(const int q, const int b, const int c, PolyhedralMesh& mesh, PolyhedralMesh& meshDual)
	{
		PolyhedralMesh meshTriangulated;
		PolyhedralMesh meshFinal;
		
		vector<int> dimension = CalcoloVEFPoliedro(q, b, c);
		vector<int> dimensionDuplicated = CalcoloDuplicato(q, b, c, dimension);
		triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
		RimuoviVerticiDuplicati(meshTriangulated);
    	RimuoviLatiDuplicati(meshTriangulated);
    	NewMesh(meshTriangulated, meshFinal, dimension);
    	map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap = buildEdgeToFacesMap(meshFinal);
    	CalculateDual(meshFinal, meshDual, edgeToFacesMap);
		PopulateCell3D(meshDual);
    }
	
	// La funzione Triangolazione2 serve per fare una triangolazione della mesh in due passaggi,
	// in modo da ottenere una mesh finale più ordinata e ben strutturata.
	//
	// Ecco cosa fa passo dopo passo:
	// 1. Prima fa una triangolazione iniziale della mesh di partenza (`mesh`) usando il livello di suddivisione `b`.
	//    - Calcola quante facce, spigoli e vertici teorici dovrebbero esserci (con `CalcoloVEFPoliedro`),
	//      e anche quelli che servono considerando eventuali duplicati (con `CalcoloDuplicato`).
	//    - Poi triangola e salva il risultato in `meshTriangulated`.
	//    - Rimuove i vertici e i lati duplicati e crea una mesh "pulita" chiamata `meshFinal`.
	//
	// 2. Poi parte da `meshFinal` e fa una seconda triangolazione migliorata:
	//    - Calcola un'altra volta le dimensioni teoriche della mesh finale (con `CalcoloDimensione2`).
	//    - Costruisce una mappa che dice per ogni lato a quali facce appartiene (serve per la connettività).
	//    - Infine triangola di nuovo (con `triangulateAndStore2`) e salva tutto in `meshTriangulated2`.
	//
	// Alla fine si ottiene una mesh triangolata bene, senza duplicati e pronta per analisi o costruzione di duali.
    
    void Triangolazione2(const int q, const int b, PolyhedralMesh& mesh, PolyhedralMesh& meshTriangulated2)
	{
		PolyhedralMesh meshTriangulated;
		PolyhedralMesh meshFinal;
		
		vector<int> dimension = CalcoloVEFPoliedro(q, b, 0);
		vector<int> dimensionDuplicated = CalcoloDuplicato(q, b, 0, dimension);
		triangulateAndStore(mesh, meshTriangulated, b, 0,  dimensionDuplicated);
		RimuoviVerticiDuplicati(meshTriangulated);
    	RimuoviLatiDuplicati(meshTriangulated);
    	NewMesh(meshTriangulated, meshFinal, dimension);
    	
    	vector<int> dimension2 = CalcoloDimensione2(b, q);
    	map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap = buildEdgeToFacesMap(meshFinal);
		triangulateAndStore2(meshFinal, meshTriangulated2, dimension2, edgeToFacesMap);
		PopulateCell3D(meshTriangulated2);
    }
	
	// popolo cella 3d con valori gia esistenti provenineti da MeshTriangulated

	void PopulateCell3D(PolyhedralMesh& meshTriangulated){
		
		meshTriangulated.Cell3DsId = {0};
		meshTriangulated.NumCells0Ds = meshTriangulated.Cell0DsId.size();
		meshTriangulated.NumCells1Ds = meshTriangulated.Cell1DsId.size();
		meshTriangulated.NumCells2Ds = meshTriangulated.Cell2DsId.size();

		for (unsigned int i = 0; i < meshTriangulated.Cell0DsId.size(); i++){
			meshTriangulated.Cell3DsVertices.push_back(meshTriangulated.Cell0DsId[i]);
		}
		for (unsigned int i = 0; i < meshTriangulated.Cell1DsId.size(); i++){
			meshTriangulated.Cell3DsEdges.push_back(meshTriangulated.Cell1DsId[i]);

		}
		for (unsigned int i = 0; i < meshTriangulated.Cell2DsId.size(); i++){
			meshTriangulated.Cell3DsFaces.push_back(i);
		}
	}
	

}