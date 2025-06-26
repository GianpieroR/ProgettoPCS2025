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
	
	// In qeusta funzipone cerco un lato (edge) tra i vertici 'a' e 'b' nella mesh triangolata.
	// Se il lato esiste già (in un ordine qualsiasi), lo riutilizza e lo associa alla faccia y3.
	// Altrimenti, lo crea:
	//  - Aggiunge gli estremi a e b alla lista dei lati (Cell1DsExtrema)
	//  - Imposta l'id del nuovo lato
	//  - Cerca un flag comune tra i due vertici e lo assegna al lato, se presente
	//  - Altrimenti, assegna un flag massimo per indicare un lato "neutro"
	//  - Aggiunge il nuovo lato alla lista dei lati della faccia y3
	//  - Aggiorna l'indice y2 per il prossimo lato da inserire
	
	void FindAddEdge(
		const int a, const int b,
		PolyhedralMesh& meshTriangulated,
		unsigned int& y2,
		const unsigned int y3)
	{
		bool found = false;
	
		for (unsigned int i = 0; i <= y2; ++i) {  // Scorri solo fino all'ultimo edge inserito
			if ((meshTriangulated.Cell1DsExtrema(i, 0) == a && meshTriangulated.Cell1DsExtrema(i, 1) == b) ||
				(meshTriangulated.Cell1DsExtrema(i, 0) == b && meshTriangulated.Cell1DsExtrema(i, 1) == a)) {
				
				// Edge già presente ⇒ usalo
				meshTriangulated.Cell2DsEdges[y3].push_back(i);
				found = true;
				break;
			}
		}
		
		if (!found) {
			// Edge non esiste ⇒ lo creiamo
			meshTriangulated.Cell1DsExtrema.row(y2) << a, b;
			meshTriangulated.Cell1DsId[y2] = y2;
		
			const auto& c = meshTriangulated.Cell0DsFlag[a];
			const auto& d = meshTriangulated.Cell0DsFlag[b];
		
			bool common = false;
			for (size_t s = 0; s < c.size(); ++s) {
				for (size_t t = 0; t < d.size(); ++t) {
					if (c[s] == d[t]) {
						meshTriangulated.Cell1DsFlag[y2] = c[s];
						common = true;
						break;
					}
				}
				if (common) break;
			}
			
			if (!common) {
				meshTriangulated.Cell1DsFlag[y2] = numeric_limits<unsigned int>::max();
			}
		
			meshTriangulated.Cell2DsEdges[y3].push_back(y2);
			++y2;
		}
	}
	
	// Funzione principale che prende ogni faccia triangolare della mesh originale
	// e la suddivide in tanti triangolini più piccoli, in base al livello di suddivisione dato (b + c).
	// Ogni nuovo vertice calcolato viene salvato con le sue coordinate e con un "flag" che lo lega agli spigoli originali.
	// Inoltre, per ogni triangolino creato, si controlla se i lati esistono già (per non duplicarli inutilmente) usando FindAddEdge,
	// altrimenti li si crea. Tutti i nuovi vertici, lati e triangoli vengono salvati dentro 'meshTriangulated'.

	// Quindi in maniera sintetica:
	// - Si inizializzano le strutture dati della mesh triangolata con le dimensioni giuste (dimensionDuplicated)
	// - Per ogni faccia triangolare della mesh originale:
	//     - Si costruisce una griglia di punti che suddivide la faccia in tanti triangolini regolari
	//     - Ogni punto viene salvato con un "flag" che indica a quali spigoli originali appartiene
	//     - Poi si costruiscono i triangolini veri e propri, assegnando a ciascuno:
	//         - i suoi 3 vertici
	//         - e i 3 lati (riutilizzando quelli già trovati, se possibile)
	// - Alla fine, la meshTriangulated contiene tutti i vertici, lati e triangoli risultanti dalla suddivisione.
	void triangulateAndStore(PolyhedralMesh& mesh, PolyhedralMesh& meshTriangulated, const unsigned int b, const unsigned int c, const vector<int>& dimensionDuplicated) {
		
		unsigned int subdivisionLevel = b+c;
		
		meshTriangulated.Cell0DsId.resize(dimensionDuplicated[0]);
		meshTriangulated.Cell0DsCoordinates = MatrixXd::Zero(3, dimensionDuplicated[0]);
		meshTriangulated.Cell0DsFlag.resize(dimensionDuplicated[0]);
		
		meshTriangulated.Cell1DsId.resize(dimensionDuplicated[1]);
		meshTriangulated.Cell1DsExtrema = MatrixXi::Zero(dimensionDuplicated[1], 2);
		meshTriangulated.Cell1DsFlag.resize(dimensionDuplicated[1]);
		
		meshTriangulated.Cell2DsId.resize(dimensionDuplicated[2]);
		meshTriangulated.Cell2DsVertices.resize(dimensionDuplicated[2]);
		meshTriangulated.Cell2DsEdges.resize(dimensionDuplicated[2]);

		unsigned int y1=0; unsigned int y2=0; unsigned int y3=0;
		
		// Costruiamo una griglia triangolare di punti all’interno di ogni faccia originale della mesh.
		// 
		// Per ogni faccia:
		//
		// - Iteriamo su `i` da 0 a subdivisionLevel, creando `i+1` punti per ogni riga (forma triangolare).
		// - Per ogni riga calcoliamo due punti `start` e `end`, che sono interpolati linearmente lungo i lati V0-V1 e V0-V2.
		// - Poi, per ogni punto `j` nella riga, calcoliamo la posizione interpolata tra `start` e `end`.
		// - Questo crea una griglia di punti distribuiti uniformemente sulla faccia.
		//
		// Per ogni punto calcolato:
		//
		// - Salviamo la coordinata del punto nella mesh triangolata.
		// - Assegniamo un ID univoco a ciascun vertice.
		// - Assegniamo un flag per indicare a quali spigoli originali appartiene il vertice:
		//   * Vertici agli angoli prendono i flag di due spigoli.
		//   * Vertici sui lati prendono il flag di uno spigolo.
		//   * Vertici interni hanno un flag speciale (max) che indica che non appartengono a spigoli.
		// - Aggiungiamo l’ID del vertice alla riga corrente della griglia.
		// - Incrementiamo il contatore dei vertici.
		
		for (unsigned int faceId = 0; faceId < mesh.Cell2DsId.size() ; ++faceId){
			const auto& face = mesh.Cell2DsVertices[faceId];
			Vector3d V0 = mesh.Cell0DsCoordinates.col(face[0]); Vector3d V1 = mesh.Cell0DsCoordinates.col(face[1]); Vector3d V2 = mesh.Cell0DsCoordinates.col(face[2]);
	
			vector<vector<int>> vertexGrid;
	
			for (unsigned int i = 0; i <= subdivisionLevel; ++i) {
				vector<int> row;
				Vector3d start = ((double)i / subdivisionLevel) * V1 + ((double)(subdivisionLevel - i) / subdivisionLevel) * V0;
				Vector3d end   = ((double)i / subdivisionLevel) * V2 + ((double)(subdivisionLevel - i) / subdivisionLevel) * V0;
	
				for (unsigned int j = 0; j <= i; ++j) {
					Vector3d point;
					if (i == 0) {
						point = V0;
					} else {
						point = ((double)j / i) * end + ((double)(i - j) / i) * start;
					}
					
					meshTriangulated.Cell0DsCoordinates.col(y1) = point;
					meshTriangulated.Cell0DsId[y1]=y1;

					if (i == 0) {
					meshTriangulated.Cell0DsFlag[y1] = {mesh.Cell2DsEdges[faceId][0], mesh.Cell2DsEdges[faceId][2]};
					}
					else if (i == subdivisionLevel) {
						if (j == 0)
							meshTriangulated.Cell0DsFlag[y1] = {mesh.Cell2DsEdges[faceId][0], mesh.Cell2DsEdges[faceId][1]};
						else if (j == subdivisionLevel)
							meshTriangulated.Cell0DsFlag[y1] = {mesh.Cell2DsEdges[faceId][1], mesh.Cell2DsEdges[faceId][2]};
						else
							meshTriangulated.Cell0DsFlag[y1] = {mesh.Cell2DsEdges[faceId][1]};
					}
					else if (j == 0) {
						meshTriangulated.Cell0DsFlag[y1] = {mesh.Cell2DsEdges[faceId][0]};
					}
					else if (j == i) {
						meshTriangulated.Cell0DsFlag[y1] = {mesh.Cell2DsEdges[faceId][2]};
					}
					else {
						meshTriangulated.Cell0DsFlag[y1] = {numeric_limits<unsigned int>::max()};
					}		
	
					row.push_back(y1);
					y1++;
				}
				vertexGrid.push_back(row);
			}
			for (unsigned int i = 0; i < subdivisionLevel; ++i) {
				for (unsigned int j = 0; j < i; ++j) {
					unsigned int v1 = vertexGrid[i][j];
					unsigned int v2 = vertexGrid[i + 1][j];
					unsigned int v3 = vertexGrid[i + 1][j + 1];
		
					vector<unsigned int> verts1 = {v1, v2, v3};
					meshTriangulated.Cell2DsVertices[y3]=verts1;
						
					meshTriangulated.Cell2DsId[y3]=y3;
						
					for (unsigned int e = 0; e < 3; ++e) {
						int a = verts1[e];
						int b = verts1[(e + 1) % 3];
						FindAddEdge(a, b, meshTriangulated, y2, y3);
					}
	
					y3++;
							
					unsigned int v4 = vertexGrid[i][j];
					unsigned int v5 = vertexGrid[i + 1][j + 1];
					unsigned int v6 = vertexGrid[i][j + 1];
		
					vector<unsigned int> verts2 = {v4, v5, v6};

					meshTriangulated.Cell2DsVertices[y3]=verts2;
					meshTriangulated.Cell2DsId[y3]=y3;
						
					for (unsigned int e = 0; e < 3; ++e) {
						int a = verts2[e];
						int b = verts2[(e + 1) % 3];
						FindAddEdge(a, b, meshTriangulated, y2, y3);
					}
					y3++;
				}
		
				// Triangolo finale in basso a sinistra
				unsigned int v1 = vertexGrid[i][i];
				unsigned int v2 = vertexGrid[i + 1][i];
				unsigned int v3 = vertexGrid[i + 1][i + 1];
	
				vector<unsigned int> verts = {v1, v2, v3};
				meshTriangulated.Cell2DsVertices[y3]=verts;

				meshTriangulated.Cell2DsId[y3]=y3;
					
				for (unsigned int e = 0; e < 3; ++e) {
					int a = verts[e];
					int b = verts[(e + 1) % 3];
					FindAddEdge(a, b, meshTriangulated, y2, y3);
				}
				y3++;
		}
			
		}

	}
	

}