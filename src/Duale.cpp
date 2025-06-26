#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <fstream>
#include <sstream>
#include <Eigen/Dense> 
#include <cmath>

using namespace std;
using namespace Eigen;
namespace PolyhedralLibrary {
	
// Funzione per ottenere il baricentro di una faccia
Vector3d getFaceBarycenter(const PolyhedralMesh& meshTriangulated, const unsigned int faceId) {
    Vector3d barycenter = Vector3d::Zero();
    const auto& faceVertices = meshTriangulated.Cell2DsVertices[faceId];
    for (unsigned int v_id : faceVertices) {
        barycenter += meshTriangulated.Cell0DsCoordinates.col(v_id);
    }
    return barycenter /= (3.0);
}

// Mappa Spigolo Originale -> Facce che lo Contengono
map <pair<unsigned int, unsigned int>, vector<unsigned int>> buildEdgeToFacesMap(const PolyhedralMesh& meshTriangulated) {
    map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFaces;

    for (unsigned int faceId = 0; faceId < meshTriangulated.Cell2DsId.size(); ++faceId) {
        const vector<unsigned int>& faceEdges = meshTriangulated.Cell2DsEdges[faceId]; // spigoli della faccia corrente
        
        for (unsigned int edgeOriginalId : faceEdges) {
			unsigned int v1_id = meshTriangulated.Cell1DsExtrema(edgeOriginalId, 0);
			unsigned int v2_id = meshTriangulated.Cell1DsExtrema(edgeOriginalId, 1);
			pair<unsigned int, unsigned int> sortedEdgeVertices = {min(v1_id, v2_id), max(v1_id, v2_id)};
			// Ordiniamo i vertici dello spigolo per avere una chiave univoca nella mappa
			edgeToFaces[sortedEdgeVertices].push_back(faceId);
        }
    }
    return edgeToFaces;
}

// Mappa Vertice Originale -> Facce che lo Contengono
map<unsigned int, vector<unsigned int>> buildVertexToFacesMap(const PolyhedralMesh& meshTriangulated) {
    map<unsigned int, vector<unsigned int>> vertexToFaces;
    for (unsigned int faceId = 0; faceId < meshTriangulated.Cell2DsId.size(); ++faceId) {
        for (unsigned int vertexOriginalId : meshTriangulated.Cell2DsVertices[faceId]) {
            vertexToFaces[vertexOriginalId].push_back(faceId);
        }
    }
    return vertexToFaces;
}

// Mappa Vertice Originale -> Spigoli che vi incidono
map<unsigned int, vector<unsigned int>> buildVertexToEdgesMap(const PolyhedralMesh& meshTriangulated) {
    map<unsigned int, vector<unsigned int>> vertexToEdges;
    for (unsigned int edgeId = 0; edgeId < meshTriangulated.Cell1DsExtrema.rows(); ++edgeId) {
        unsigned int v1 = meshTriangulated.Cell1DsExtrema(edgeId, 0);
        unsigned int v2 = meshTriangulated.Cell1DsExtrema(edgeId, 1);
        vertexToEdges[v1].push_back(edgeId);
        vertexToEdges[v2].push_back(edgeId);
    }
    return vertexToEdges;
}


void CalculateDual(PolyhedralMesh& meshTriangulated, PolyhedralMesh& meshDual, map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap)
{
	// VERTICI
	// i vertici del duale sono i baricentri delle facce originali
	meshDual.Cell0DsId.resize(meshTriangulated.Cell2DsId.size());
	meshDual.Cell0DsCoordinates = MatrixXd::Zero(3, meshTriangulated.Cell2DsId.size());
	
	for (unsigned int faceId = 0; faceId < meshTriangulated.Cell2DsId.size() ; ++faceId){
		meshDual.Cell0DsCoordinates.col(faceId) = getFaceBarycenter(meshTriangulated, faceId);
        meshDual.Cell0DsId[faceId] = faceId; // L'ID del vertice duale è l'ID della faccia originale
    }
    
    // SPIGOLI
    
    vector<pair<unsigned int, unsigned int>> dualEdgesExtremaVector;
    // Ogni spigolo interno del poliedro originale (condiviso da due facce) genera uno spigolo nel poliedro duale che connette i baricentri di quelle due facce.
    // Useremo un vettore temporaneo e poi lo convertiremo in Eigen::MatrixXi.
        
    for (const auto& entry : edgeToFacesMap) {
        const vector<unsigned int>& facesSharingEdge = entry.second; // facce
        if (facesSharingEdge.size() == 2) {
			unsigned int faceId1 = facesSharingEdge[0]; // faccia condivisa 1
			unsigned int faceId2 = facesSharingEdge[1]; // faccia condivisa 2
			pair<unsigned int, unsigned int> dualEdge = {min(faceId1, faceId2), max(faceId1, faceId2)};
			// Gli ID delle facce originali diventano gli ID dei vertici del duale
			dualEdgesExtremaVector.push_back(dualEdge);
		}
	}

	// Questa struttura intermedia mi serve perchè non so a prescindere le dimensioni di meshDual.Cell1DsExtrema
	// il numero di spigoli duali è pari al numero di spigoli interni della mesh originale
		
	meshDual.Cell1DsId.resize(dualEdgesExtremaVector.size());
	meshDual.Cell1DsExtrema = MatrixXi::Zero(dualEdgesExtremaVector.size(), 2);
	
	for (unsigned int i = 0; i < dualEdgesExtremaVector.size(); ++i) {
		meshDual.Cell1DsId[i] = i;
		meshDual.Cell1DsExtrema(i, 0) = dualEdgesExtremaVector[i].first;
		meshDual.Cell1DsExtrema(i, 1) = dualEdgesExtremaVector[i].second;
    }
    
    // FACCE
    
    map<pair<unsigned int, unsigned int>, unsigned int> dualEdgeToIdMap;
    // mappa che per ogni coppia di id di facce originali mi associa l'id dello spigolo duale
    for (unsigned int i = 0; i < meshDual.Cell1DsId.size(); ++i) {
        unsigned int v1 = meshDual.Cell1DsExtrema(i, 0);
        unsigned int v2 = meshDual.Cell1DsExtrema(i, 1);
        dualEdgeToIdMap[{min(v1, v2), max(v1, v2)}] = i;
    }
    	 
    map<unsigned int, vector<unsigned int>> vertexToFacesMap = buildVertexToFacesMap(meshTriangulated); 
    map<unsigned int, vector<unsigned int>> vertexToEdgesMap = buildVertexToEdgesMap(meshTriangulated);

    meshDual.Cell2DsVertices.resize(vertexToFacesMap.size());
    meshDual.Cell2DsId.resize(vertexToFacesMap.size());
    meshDual.Cell2DsEdges.resize(vertexToFacesMap.size());

    unsigned int dualFaceIdCounter = 0; // contatore facce duali
    for (const auto& entry : vertexToFacesMap) {
        unsigned int vertexOriginalId = entry.first; // Il vertice originale che "genera" questa faccia duale
        vector<unsigned int> incidentFaces = entry.second;  // Le facce originali incidenti a questo vertice
        vector<unsigned int> incidentEdges = vertexToEdgesMap[vertexOriginalId]; // Gli spigoli originali incidenti a questo vertice

        // Un vertice deve essere incidente ad almeno 3 facce/spigoli per formare una faccia duale chiusa.
        if (incidentFaces.size() < 3 || incidentEdges.size() < 3) {
            cerr << "Warning: Vertice originale " << vertexOriginalId << " incidente a meno di 3 facce o spigoli. Non può formare una faccia duale chiusa." << endl;
            continue; // Salta la creazione della faccia duale
        }

        vector<unsigned int> orderedDualFaceVertices; // Vertici della faccia duale, ordinati
        vector<unsigned int> dualFaceEdges;           // Spigoli della faccia duale, ordinati

        // Iniziamo il percorso topologico:
        // Scegliamo uno spigolo incidente qualsiasi e una delle facce che lo contiene come punto di partenza.
        unsigned int currentEdgeOriginalId = incidentEdges[0]; // lato di partenza

        // Troviamo le due facce che condividono questo spigolo originale.
        pair<unsigned int, unsigned int> currentEdgeKey = {
            min(meshTriangulated.Cell1DsExtrema(currentEdgeOriginalId, 0), meshTriangulated.Cell1DsExtrema(currentEdgeOriginalId, 1)),
            max(meshTriangulated.Cell1DsExtrema(currentEdgeOriginalId, 0), meshTriangulated.Cell1DsExtrema(currentEdgeOriginalId, 1))
        };
        vector<unsigned int> facesSharingStartEdge = edgeToFacesMap.at(currentEdgeKey); 
        // Se la currentEdgeKey non esiste nella edgeToFacesMap, at() lancerà un'eccezione di tipo std::out_of_range
        // L'operatore [], se la chiave non esiste, inserirà automaticamente una nuova coppia chiave-valore nella mappa, 
        // con il valore predefinito per il tipo di valore

        // Ci assicuriamo che lo spigolo sia condiviso da esattamente due facce
        if (facesSharingStartEdge.size() != 2) {
             cerr << "Errore: Spigolo originale " << currentEdgeOriginalId << " non condiviso da esattamente due facce" << endl;
             continue; // Salta, non possiamo costruire la faccia duale correttamente
        }

        // Scelgiamo una delle due facce come punto di partenza per il ciclo
        unsigned int startFaceId = facesSharingStartEdge[0];
        unsigned int previousDualVertexId = startFaceId; // Il primo vertice della faccia duale (id della faccia originale)

        // Aggiungiamo il primo vertice duale della faccia che stiamo costruendo
        orderedDualFaceVertices.push_back(previousDualVertexId);

        // Percorriamo il ciclo attorno al vertice originale
        for (size_t i = 0; i < incidentEdges.size(); ++i) { // Il numero di spigoli incidenti è il numero di lati della faccia duale
            // Troviamo le due facce che condividono l'attuale spigolo `currentEdgeOriginalId`
            vector<unsigned int> currentEdgeFaces = edgeToFacesMap.at(currentEdgeKey);

            unsigned int face1 = currentEdgeFaces[0];
            unsigned int face2 = currentEdgeFaces[1];

            // La faccia "precedente" è previousDualVertexId. Troviamo la faccia "successiva".
            unsigned int nextDualVertexId;
            if (face1 == previousDualVertexId) {
                nextDualVertexId = face2;
            } else if (face2 == previousDualVertexId) {
                nextDualVertexId = face1;
            } else {
                // Questa situazione non dovrebbe accadere se la logica di navigazione è corretta
                cerr << "Errore logico: Faccia corrente non trovata tra le facce che condividono lo spigolo." << endl;
                break;
            }
            
            // Troviamo e aggiungiamo lo spigolo duale tra il vertice duale corrente e il successivo
            pair<unsigned int, unsigned int> dualEdgeToFind = {min(previousDualVertexId, nextDualVertexId), max(previousDualVertexId, nextDualVertexId)};
            auto it_dual_edge = dualEdgeToIdMap.find(dualEdgeToFind);

            if (it_dual_edge != dualEdgeToIdMap.end()) {
                dualFaceEdges.push_back(it_dual_edge->second); 
                // lo spigolo viene trovato e quindi aggiunto a dualFaceEdges
            } else {
                cerr << "Errore: Spigolo duale non trovato per facce (vertici duali) "
                          << previousDualVertexId << " e " << nextDualVertexId << " attorno al vertice originale "
                          << vertexOriginalId << endl;
                          // Significa che non esiste uno spigolo nella mesh originale che colleghi le due facce currentDualVertexId e
                          // nextDualVertexId (che sono i vertici della nostra faccia duale) e che passi attraverso vertexOriginalId
                break;
            }

            // Se abbiamo raggiunto il punto di partenza (per chiudere il poligono)
            if (nextDualVertexId == startFaceId) {
                if (orderedDualFaceVertices.size() == incidentFaces.size()) {
                    // Abbiamo aggiunto tutti i vertici unici, possiamo chiudere il ciclo.
                    break;
                } else {
                    // Siamo tornati all'inizio ma non abbiamo visitato tutte le facce incidenti.
                    cerr << "Errore topologico: Ciclo incompleto per vertice " << vertexOriginalId << endl;
                    break;
                }
            }
			
			// Aggiungiamo il prossimo vertice duale (id della faccia originale)
			bool alreadyAdded = false;
            for(unsigned int v : orderedDualFaceVertices) {
                if (v == nextDualVertexId) {
                    alreadyAdded = true;
                    break;
                }
            }
            if (!alreadyAdded) {
                orderedDualFaceVertices.push_back(nextDualVertexId);
            }
            
            // Aggiorniamo previousDualVertexId per la prossima iterazione
            previousDualVertexId = nextDualVertexId;
			
            // Troviamo il prossimo spigolo originale incidente a vertexOriginalId (vertice originale attorno al quale stiamo costruendo la faccia duale)
            // che connette currentDualVertexId (la faccia precedente) con nextDualVertexId (la faccia successiva).
            // Dobbiamo trovare lo spigolo della faccia `nextDualVertexId` che è incidente a `vertexOriginalId`
            // e che non è `currentEdgeOriginalId`.
            bool foundNextEdgeForCycle = false;
            const auto& nextFaceEdges = meshTriangulated.Cell2DsEdges[nextDualVertexId]; // spigoli della faccia originale che abbiamoa appena raggiunto
            for (unsigned int edgeOfNextFace : nextFaceEdges) {
                if (edgeOfNextFace == currentEdgeOriginalId) continue; // Non prendere lo spigolo da cui siamo venuti

                // Controlliamo se edgeOfNextFace contiene vertexOriginalId
                unsigned int v_e1 = meshTriangulated.Cell1DsExtrema(edgeOfNextFace, 0);
                unsigned int v_e2 = meshTriangulated.Cell1DsExtrema(edgeOfNextFace, 1);

                if (v_e1 == vertexOriginalId || v_e2 == vertexOriginalId) {
                    currentEdgeOriginalId = edgeOfNextFace; // Questo è il prossimo spigolo da percorrere
                    currentEdgeKey = {min(v_e1, v_e2), max(v_e1, v_e2)};
                    foundNextEdgeForCycle = true;
                    break;
                }
            }
            if (!foundNextEdgeForCycle) {
                cerr << "Errore: Impossibile trovare il prossimo spigolo per il ciclo attorno al vertice " << vertexOriginalId << endl;
                break;
            }
        } // Fine del ciclo di costruzione della singola faccia duale

        // Assegniamo la faccia duale costruita
        if (!orderedDualFaceVertices.empty()) {
            meshDual.Cell2DsId[dualFaceIdCounter] = dualFaceIdCounter;
            meshDual.Cell2DsVertices[dualFaceIdCounter] = orderedDualFaceVertices;
            meshDual.Cell2DsEdges[dualFaceIdCounter] = dualFaceEdges;
            dualFaceIdCounter++;
        }
	}
}

void ProjectMeshToUnitSphere(PolyhedralMesh& meshTriangulated) {
    for (int i = 0; i < meshTriangulated.Cell0DsCoordinates.cols(); ++i) {
        Vector3d vertexCoords = meshTriangulated.Cell0DsCoordinates.col(i);
        double norm = vertexCoords.norm(); // Equivalente a std::sqrt(vertexCoords.squaredNorm());
        if (norm < 1e-12) { 
            cerr << "Warning: Vertice " << i << " troppo vicino all'origine. Non proiettato." << endl;
            continue; // Salta la proiezione per questo vertice
        }
        meshTriangulated.Cell0DsCoordinates.col(i) = vertexCoords / norm;
    }
}

}



