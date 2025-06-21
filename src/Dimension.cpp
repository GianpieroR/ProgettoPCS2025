#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <limits> 
#include <Eigen/Dense>
#include <set> 
#include <cmath>

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{
	void invertiValori(int& a, int& b) {
		int temp = a;
		a = b;
		b = temp;
	}
	
	vector<int> ComputePolyhedronVEF(int q, int b, int c) {
    // Calcola il valore T in base ai parametri b e c
    int T = b * b + b * c + c * c;

    int vertices = 0;
    int edges = 0;
    int faces = 0;

    // Determina il numero di vertici, spigoli e facce a seconda del valore di q
    switch (q) {
        case 3:
            vertices = 2 * T + 2;
            edges = 6 * T;
            faces = 4 * T;
            break;
        case 4:
            vertices = 4 * T + 2;
            edges = 12 * T;
            faces = 8 * T;
            break;
        default: // per q >= 5 o altri valori
            vertices = 10 * T + 2;
            edges = 30 * T;
            faces = 20 * T;
            break;
    }

    // Inserisce i risultati in un vettore e lo restituisce
    return {vertices, edges, faces};
}
	
	vector<int> CalculateDimension2(int b, int q)
	{
		vector<int> result(3);
		int V = 0;
		int E = 0;
		int F = 0;
		
		if (q == 3) {
			V = 4 + 6 * (2 * b - 1) + static_cast<int>(round(4 * ((3.0 * b * b) / 2.0 - (3.0 * b) / 2.0 + 1)));
			E = 6 * 2 * b + static_cast<int>(round(4 * ((9.0 * b * b) / 2.0 + (3.0 * b) / 2.0)));
			F = 4 * (3 * b * b + 3 * b);
		}
		else if (q == 4) {
			V = 6 + 12 * (2 * b - 1) + static_cast<int>(round(8 * ((3.0 * b * b) / 2.0 - (3.0 * b) / 2.0 + 1)));
			E = 12 * 2 * b + static_cast<int>(round(8 * ((9.0 * b * b) / 2.0 + (3.0 * b) / 2.0)));
			F = 8 * (3 * b * b + 3 * b);
		}
		else {
			V = 12 + 30 * (2 * b - 1) + static_cast<int>(round(20 * ((3.0 * b * b) / 2.0 - (3.0 * b) / 2.0 + 1)));
			E = 30 * 2 * b + static_cast<int>(round(20 * ((9.0 * b * b) / 2.0 + (3.0 * b) / 2.0)));
			F = 20 * (3 * b * b + 3 * b);
		}

		result[0] = V;  
		result[1] = E;  
		result[2] = F;
		
		return result;
	}

	vector<int> CalculateDuplicated(int q, int b, int c, const vector<int>& dimension)
	{
		vector<int> result(3);
		int subdivisionLevel = 0;
		subdivisionLevel = b + c;
		int V = dimension[0];
		int E = dimension[1];
		
		if (q == 3) {
			V += 2*4 + 6*(subdivisionLevel - 1);
			E += 6*subdivisionLevel;
		}
		else if (q == 4) {
			V += 3*6 + 12*(subdivisionLevel - 1);
			E += 12*subdivisionLevel;
		}
		else {
			V += 4*12 + 30*(subdivisionLevel - 1);
			E += 30* subdivisionLevel;
		}
		result[0] = V;  
		result[1] = E;  
		result[2] = dimension[2]; //il numero delle facce non cambia
		
		return result;
	}
	
	void RemoveDuplicatedVertices(PolyhedralMesh& meshTriangulated) {
    const double tolerance = 1e-12;
    const unsigned int INVALID_FLAG = numeric_limits<unsigned int>::max();
    const size_t vertexCount = meshTriangulated.Cell0DsCoordinates.cols();

    // Array per mappare ogni vertice al suo "vertice master" finale.
    // Inizialmente ogni vertice si mappa su se stesso.
    vector<unsigned int> vertexMapping(vertexCount);
    for (unsigned int idx = 0; idx < vertexCount; ++idx) {
        vertexMapping[idx] = idx;
    }

    // Booleano per identificare se un vertice può ancora essere considerato come master
    // I duplicati verranno segnati come false.
    vector<bool> canBeMaster(vertexCount, true);

    // Iteriamo su ogni vertice per verificare duplicati
    for (size_t i = 0; i < vertexCount; ++i) {
        // Se il vertice i è già stato confermato come master (INVALID_FLAG)
        // o non è più un candidato master, saltalo.
        if (meshTriangulated.Cell0DsFlag[i][0] == INVALID_FLAG || !canBeMaster[i])
            continue;

        // Confronta il vertice i con tutti quelli successivi
        for (size_t j = i + 1; j < vertexCount; ++j) {
            if (meshTriangulated.Cell0DsFlag[j][0] == INVALID_FLAG || !canBeMaster[j])
                continue;

            // Verifica se i due vertici condividono almeno una faccia comune
            bool sharedFace = false;
            for (unsigned int faceI : meshTriangulated.Cell0DsFlag[i]) {
                for (unsigned int faceJ : meshTriangulated.Cell0DsFlag[j]) {
                    if (faceI == faceJ) {
                        sharedFace = true;
                        break;
                    }
                }
                if (sharedFace) break;
            }

            // Se hanno una faccia in comune e sono praticamente sovrapposti (distanza < tolleranza)
            if (sharedFace) {
                double distance = (meshTriangulated.Cell0DsCoordinates.col(i) - meshTriangulated.Cell0DsCoordinates.col(j)).norm();
                if (distance < tolerance) {
                    // j diventa il master, i è duplicato di j
                    canBeMaster[i] = false;

                    // Trova la radice del mapping per i e j
                    unsigned int root_i = i;
                    while (vertexMapping[root_i] != root_i) {
                        root_i = vertexMapping[root_i];
                    }
                    unsigned int root_j = j;
                    while (vertexMapping[root_j] != root_j) {
                        root_j = vertexMapping[root_j];
                    }

                    // Unisci i due gruppi di vertici duplicati se sono diversi
                    if (root_i != root_j) {
                        vertexMapping[root_i] = root_j;
                    }
                    // Assicura che i punti a j
                    vertexMapping[i] = root_j;

                    // Sincronizza la posizione di i con quella di j per coerenza
                    meshTriangulated.Cell0DsCoordinates.col(i) = meshTriangulated.Cell0DsCoordinates.col(j);
                }
            }
        }
    }

    // Passo finale di compressione: Assicura che tutti i vertici puntino alla loro radice finale
    for (unsigned int k = 0; k < vertexCount; ++k) {
        unsigned int current = k;
        while (vertexMapping[current] != current) {
            current = vertexMapping[current];
        }
        vertexMapping[k] = current;
    }

    // Aggiorna gli spigoli con i nuovi indici dei vertici
    for (int edgeId = 0; edgeId < meshTriangulated.Cell1DsExtrema.rows(); ++edgeId) {
        meshTriangulated.Cell1DsExtrema(edgeId, 0) = vertexMapping[meshTriangulated.Cell1DsExtrema(edgeId, 0)];
        meshTriangulated.Cell1DsExtrema(edgeId, 1) = vertexMapping[meshTriangulated.Cell1DsExtrema(edgeId, 1)];
    }

    // Aggiorna i vertici di ogni faccia secondo il nuovo mapping
    for (unsigned int faceId = 0; faceId < meshTriangulated.Cell2DsId.size(); ++faceId) {
        for (unsigned int& vertexId : meshTriangulated.Cell2DsVertices[faceId]) {
            vertexId = vertexMapping[vertexId];
        }
    }

    // Infine aggiorna gli ID e i flag dei vertici nel meshTriangulated
    for (unsigned int idx = 0; idx < vertexCount; ++idx) {
        meshTriangulated.Cell0DsId[idx] = vertexMapping[idx];
        if (vertexMapping[idx] == idx) {
            // Vertice master: assegna il flag speciale
            meshTriangulated.Cell0DsFlag[idx] = {INVALID_FLAG};
        }
    }
}

	void RemoveDuplicatedEdges(PolyhedralMesh& meshTriangulated) {
    const double tolerance = 1e-12;
    const unsigned int MAX_FLAG = numeric_limits<unsigned int>::max();
    const size_t edgeCount = meshTriangulated.Cell1DsExtrema.rows();  // Numero totale di lati

    // Controlla e inizializza il vettore dei flag originali per i lati,
    // usato per conservare lo stato originale (bordo o interno).
    if (meshTriangulated.Cell1DsOriginalFlag.empty() || meshTriangulated.Cell1DsOriginalFlag.size() != edgeCount) {
        meshTriangulated.Cell1DsOriginalFlag.resize(edgeCount);
    }

    // Salva lo stato originale di ogni lato: true se era flagged con maxFlag (ad esempio bordo), false altrimenti
    for (size_t k = 0; k < edgeCount; ++k) {
        meshTriangulated.Cell1DsOriginalFlag[k] = (meshTriangulated.Cell1DsFlag[k] == MAX_FLAG);
    }

    // --- FASE 1: Setup struttura di mappatura ---
    // Ogni lato inizialmente si mappa a se stesso (non ancora duplicato)
    vector<unsigned int> edgeMapping(edgeCount);
    for (unsigned int idx = 0; idx < edgeCount; ++idx) {
        edgeMapping[idx] = idx;
    }

    // Array di booleani per tracciare se un lato può essere ancora considerato "master"
    vector<bool> candidateMaster(edgeCount, true);

    // --- FASE 2: Ricerca e mappatura dei duplicati ---
    for (size_t i = 0; i < edgeCount; ++i) {
        // Ignora lati già marcati come master o già reindirizzati (duplicati)
        if (meshTriangulated.Cell1DsFlag[i] == MAX_FLAG || !candidateMaster[i]) {
            continue;
        }

        // Confronta con ogni lato successivo per trovare duplicati
        for (size_t j = i + 1; j < edgeCount; ++j) {
            if (meshTriangulated.Cell1DsFlag[j] == MAX_FLAG || !candidateMaster[j]) {
                continue;
            }

            // Due lati sono duplicati se:
            // 1) Hanno lo stesso flag
            // 2) Gli estremi corrispondono (in ordine diretto o invertito) entro la tolleranza

            if (meshTriangulated.Cell1DsFlag[i] == meshTriangulated.Cell1DsFlag[j]) {
                int i_start = meshTriangulated.Cell1DsExtrema(i, 0);
                int i_end = meshTriangulated.Cell1DsExtrema(i, 1);
                int j_start = meshTriangulated.Cell1DsExtrema(j, 0);
                int j_end = meshTriangulated.Cell1DsExtrema(j, 1);

                bool directMatch = ((meshTriangulated.Cell0DsCoordinates.col(i_start) - meshTriangulated.Cell0DsCoordinates.col(j_start)).norm() < tolerance) &&
                                   ((meshTriangulated.Cell0DsCoordinates.col(i_end) - meshTriangulated.Cell0DsCoordinates.col(j_end)).norm() < tolerance);

                bool inverseMatch = ((meshTriangulated.Cell0DsCoordinates.col(i_start) - meshTriangulated.Cell0DsCoordinates.col(j_end)).norm() < tolerance) &&
                                    ((meshTriangulated.Cell0DsCoordinates.col(i_end) - meshTriangulated.Cell0DsCoordinates.col(j_start)).norm() < tolerance);

                if (directMatch || inverseMatch) {
                    // Identificato un duplicato: 'i' deve essere reindirizzato a 'j'
                    candidateMaster[i] = false;

                    // Ricerca della radice (master finale) per i e j (Union-Find)
                    unsigned int root_i = i;
                    while (edgeMapping[root_i] != root_i) {
                        root_i = edgeMapping[root_i];
                    }
                    unsigned int root_j = j;
                    while (edgeMapping[root_j] != root_j) {
                        root_j = edgeMapping[root_j];
                    }

                    // Se sono diversi, uniscili, collegando root_i a root_j
                    if (root_i != root_j) {
                        edgeMapping[root_i] = root_j;
                    }
                    // Assicura che 'i' punti direttamente al master di 'j'
                    edgeMapping[i] = root_j;
                }
            }
        }
    }

    // --- FASE 3: Compressione dei percorsi per aggiornare i remap ---
    for (unsigned int k = 0; k < edgeCount; ++k) {
        unsigned int current = k;
        while (edgeMapping[current] != current) {
            current = edgeMapping[current];
            // Qui si potrebbe implementare la path compression attiva, se necessario:
            // edgeMapping[k] = current;
        }
        edgeMapping[k] = current;
    }

    // --- FASE 4: Aggiornamento della mesh con i nuovi riferimenti ---
    // Aggiorna gli ID dei lati nelle facce
    for (unsigned int faceId = 0; faceId < meshTriangulated.Cell2DsId.size(); ++faceId) {
        for (unsigned int& edgeId : meshTriangulated.Cell2DsEdges[faceId]) {
            edgeId = edgeMapping[edgeId];
        }
    }

    // Aggiorna gli ID e i flag dei lati in base al remapping
    for (unsigned int k = 0; k < edgeCount; ++k) {
        meshTriangulated.Cell1DsId[k] = edgeMapping[k];
        if (edgeMapping[k] == k) {
            // Questo lato è master o non è stato reindirizzato: assegna il flag massimo
            meshTriangulated.Cell1DsFlag[k] = MAX_FLAG;
        } else {
            // Lato duplicato: il flag rimane invariato o si può gestire diversamente se serve
        }
    }
}
	
	void NewMesh(PolyhedralMesh& meshTriangulated, PolyhedralMesh& meshFinal, const vector<int>& dimension)
	{
		unsigned int maxFlag = numeric_limits<unsigned int>::max();
		meshFinal.Cell0DsId.resize(dimension[0]);
		meshFinal.Cell0DsCoordinates = MatrixXd::Zero(3, dimension[0]);
		
		meshFinal.Cell1DsId.resize(dimension[1]);
		meshFinal.Cell1DsExtrema = MatrixXi::Zero(dimension[1], 2);
		meshFinal.Cell1DsOriginalFlag.resize(dimension[1]);
		
		meshFinal.Cell2DsId.resize(dimension[2]);
		meshFinal.Cell2DsVertices.resize(dimension[2]);
		meshFinal.Cell2DsEdges.resize(dimension[2]);
		
		// VERICI
    	map<unsigned int, unsigned int> oldToNewVertexIdMap; // Mappa per tradurre vecchi ID vertice -> nuovi ID vertice
    	vector<unsigned int> uniqueOldVertexIdsInOrder; // Tiene traccia degli id originali dei vertici unici
    	
    	unsigned int k1 = 0; // Contatore per i nuovi ID dei vertici
    	for (unsigned int i = 0; i < meshTriangulated.Cell0DsCoordinates.cols(); ++i) {
			if (meshTriangulated.Cell0DsFlag[i][0] == maxFlag) { 
				oldToNewVertexIdMap[i] = k1; // Mappa il vecchio ID unico al nuovo
				uniqueOldVertexIdsInOrder.push_back(i); // Salva il vecchio ID per recuperare le coordinate dopo
				k1++; // Incrementa il contatore del nuovo ID
			}
		}

		// Ora popola le strutture finali per i vertici
		for (unsigned int i = 0; i < k1; ++i) {
			meshFinal.Cell0DsId[i] = i; // Nuovi ID consecutivi
			meshFinal.Cell0DsCoordinates.col(i) = meshTriangulated.Cell0DsCoordinates.col(uniqueOldVertexIdsInOrder[i]);
		}
		
		// SPIGOLI
		map<unsigned int, unsigned int> oldToNewEdgeIdMap; // Mappa per tradurre vecchi ID spigolo -> nuovi ID spigolo
		vector<pair<unsigned int, unsigned int>> temp_edge_extrema; // Vettore temporaneo per gli estremi degli spigoli (con i nuovi ID vertici)
		vector<bool> temp_edge_original_max_flag; // Temp per il nuovo flag
	
		unsigned int k2 = 0; // Contatore per i nuovi ID degli spigoli
		for (unsigned int i = 0; i < meshTriangulated.Cell1DsExtrema.rows(); ++i) {
			if (meshTriangulated.Cell1DsFlag[i] == maxFlag) {
				unsigned int old_v1_id = meshTriangulated.Cell1DsExtrema(i, 0);
				unsigned int old_v2_id = meshTriangulated.Cell1DsExtrema(i, 1);
	
				// Ottieni i nuovi ID dei vertici.
				// Se uno dei vertici originali non è stato mantenuto, lo spigolo non è valido.
				// oldToNewVertexIdMap.count() restituisce 1 se la chiave esiste, 0 altrimenti.
				if (oldToNewVertexIdMap.count(old_v1_id) > 0 && oldToNewVertexIdMap.count(old_v2_id) > 0) {
					unsigned int new_v1_id = oldToNewVertexIdMap[old_v1_id];
					unsigned int new_v2_id = oldToNewVertexIdMap[old_v2_id];
	
					// Ordina per consistenza (min, max) per la chiave se usata in future mappe
					temp_edge_extrema.push_back({min(new_v1_id, new_v2_id), max(new_v1_id, new_v2_id)});
					oldToNewEdgeIdMap[i] = k2; // Mappa il vecchio ID spigolo al nuovo
					temp_edge_original_max_flag.push_back(meshTriangulated.Cell1DsOriginalFlag[i]);
					k2++;
				} 
			}
		}
		// Popola le strutture finali per gli spigoli

		for (unsigned int i = 0; i < k2; ++i) {
			meshFinal.Cell1DsId[i] = i; // Nuovi ID consecutivi
			meshFinal.Cell1DsExtrema(i, 0) = temp_edge_extrema[i].first;
			meshFinal.Cell1DsExtrema(i, 1) = temp_edge_extrema[i].second;
			meshFinal.Cell1DsOriginalFlag[i] = temp_edge_original_max_flag[i];
		}
		
		// FACCE
		
		vector<vector<unsigned int>> temp_face_vertices; // Vettore temporaneo per i vertici delle facce (con i nuovi ID vertici)
		vector<vector<unsigned int>> temp_face_edges; // Vettore temporaneo per gli spigoli delle facce (con i nuovi ID spigoli)
		
		unsigned int k3 = 0; // Contatore per i nuovi ID delle facce
		for (unsigned int i = 0; i < meshTriangulated.Cell2DsId.size(); ++i) { // Copia e aggiorna i vertici della faccia
			vector<unsigned int> current_face_new_vertices;
			bool face_valid_vertices = true;
			for (unsigned int old_vertex_id : meshTriangulated.Cell2DsVertices[i]) {
				if (oldToNewVertexIdMap.count(old_vertex_id) > 0) { // Se il vertice è stato mantenuto
					current_face_new_vertices.push_back(oldToNewVertexIdMap[old_vertex_id]);
				} else {
					face_valid_vertices = false; // Vertice non valido, la faccia è invalida
					break;
				}
			}
	
			// Copia e aggiorna gli spigoli della faccia
			vector<unsigned int> current_face_new_edges;
			bool face_valid_edges = true;
			for (unsigned int old_edge_id : meshTriangulated.Cell2DsEdges[i]) {
				if (oldToNewEdgeIdMap.count(old_edge_id) > 0) { // Se lo spigolo è stato mantenuto
					current_face_new_edges.push_back(oldToNewEdgeIdMap[old_edge_id]);
				} else {
					face_valid_edges = false; // Spigolo non valido, la faccia è invalida
					break;
				}
			}
			
			// Se la faccia è valida (tutti i suoi vertici e spigoli sono stati mantenuti e riassegnati)
			if (face_valid_vertices && face_valid_edges) {
				temp_face_vertices.push_back(current_face_new_vertices);
				temp_face_edges.push_back(current_face_new_edges);
				k3++;
			}
		}
	
		// Popola le strutture finali per le facce
		for (unsigned int i = 0; i < k3; ++i) {
			meshFinal.Cell2DsId[i] = i; // Nuovi ID consecutivi
			meshFinal.Cell2DsVertices[i] = temp_face_vertices[i];
			meshFinal.Cell2DsEdges[i] = temp_face_edges[i];
		}
		
		
	}
}