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
	
	// Calcolo i VEF dei poliedri, funzione semplice e chiara
	
	vector<int> CalcoloVEFPoliedro(const int q, const int b, const int c)
	{
		vector<int> result(3); // inizializza un vettore con 3 valori, tutti -1 
	
		int T = 0;
		T = b * b + b * c + c * c;
		int V = 0;
		int E = 0;
		int F = 0;
	
		if (q == 3) {
			V = 2 * T + 2;
			E = 6 * T;
			F = 4 * T;
		}
		else if (q == 4) {
			V = 4 * T + 2;
			E = 12 * T;
			F = 8 * T;
		}
		else {
			V = 10 * T + 2;
			E = 30 * T;
			F = 20 * T;
		}
	
		result[0] = V;  // Primo elemento: V (vertici) 
		result[1] = E;  // Secondo elemento: E (spigoli)
		result[2] = F;  // Terzo elemento: F (facce)
	
		return result;  // Restituisce il vettore con i valori di V, E, F
	}
	
	// CalcoloDimensione2 calcola una stima teorica del numero di vertici (V), spigoli (E)
    // e facce (F) di una mesh derivata da un poliedro geodetico di tipo specificato da q,
    // con un livello di suddivisione indicato dal parametro b.
    // Il calcolo si basa su formule matematiche approssimate che dipendono dal valore di q (3,4 o 5)
    // e dal livello di suddivisione b. Il risultato è un vettore contenente V, E e F, che rappresentano
    // le dimensioni teoriche della mesh prima della sua effettiva generazione.
	
	vector<int> CalcoloDimensione2(const int b, const int q)
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

		result[0] = V; result[1] = E; result[2] = F;
		
		return result;
	}
	
	// Funzione che calcola il numero aggiornato di vertici (V) e spigoli (E)
    // dopo una duplicazione della mesh, mantenendo invariato il numero di facce (F).
    // La duplicazione dipende dal livello di suddivisione `b + c` e dal tipo di poliedro `q`.
    // I valori iniziali di V, E e F vengono passati tramite il parametro `dimension`.

	vector<int> CalcoloDuplicato(const int q, const int b, const int c, const vector<int>& dimension)
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
	
	// Funzione che rimuove i vertici duplicati nella mesh triangolata.
    // Due vertici sono considerati duplicati se hanno coordinate quasi coincidenti
    // (entro una certa tolleranza) e condividono almeno una cella 1D (spigolo).
	
	void RimuoviVerticiDuplicati(PolyhedralMesh& meshTriangulated)
{
    double tol = 1e-12;
    unsigned int maxFlag = numeric_limits<unsigned int>::max();
    size_t n = meshTriangulated.Cell0DsCoordinates.cols(); 

    // Vettore per tenere traccia del reindirizzamento finale degli ID dei vertici.
    // Inizialmente, ogni vertice punta a se stesso.
    vector<unsigned int> id_remap(n);
    for (unsigned int k = 0; k < n; ++k) {
        id_remap[k] = k;
    }

    // Vettore che tiene traccia di quali vertici sono ancora candidati a essere non duplicati. 
    // Se un vertice viene identificato come duplicato di un altro, la sua flag is_master_candidate viene impostata a false
    vector<bool> is_master_candidate(n, true);

    for (size_t i = 0; i < n; ++i) {
        // Se il vertice 'i' è già stato marcato per essere tenuto (maxFlag)
        // o se è già stato reindirizzato da un vertice precedente nel ciclo
        // (il che significa che è un duplicato), salta.
        if (meshTriangulated.Cell0DsFlag[i][0] == maxFlag || !is_master_candidate[i])
            continue;

        for (size_t j = i + 1; j < n; ++j) {
            if (meshTriangulated.Cell0DsFlag[j][0] == maxFlag || !is_master_candidate[j])
                continue;

            bool commonSide = false;
            for (unsigned int fi : meshTriangulated.Cell0DsFlag[i]) {
                for (unsigned int fj : meshTriangulated.Cell0DsFlag[j]) {
                    if (fi == fj) {
                        commonSide = true;
                        break;
                    }
                }
                if (commonSide) break;
            }

            if (commonSide) {
                if ((meshTriangulated.Cell0DsCoordinates.col(i) - meshTriangulated.Cell0DsCoordinates.col(j)).norm() < tol) {
                    // Il vertice 'i' è un duplicato di 'j', vogliamo mantenere 'j' e reindirizzare 'i' a 'j'.

                    // Marchiamo 'i' come duplicato
                    is_master_candidate[i] = false;

                    // Reindirizziamo 'i' a 'j'.
                    // tutti i vertici che precedentemente puntavano a 'i'
                    // ora devono puntare a 'j'.
                    // E 'i' stesso punterà a 'j'.
                    unsigned int root_i = i;
                    while (id_remap[root_i] != root_i) {
                        root_i = id_remap[root_i];
                    }
                    unsigned int root_j = j;
                    while (id_remap[root_j] != root_j) {
                        root_j = id_remap[root_j];
                    }

                    if (root_i != root_j) {
                         id_remap[root_i] = root_j;
                    }
                    // vertice 'i' (e i suoi precedenti duplicati)
                    // puntano a 'j' o al suo master.
                    id_remap[i] = root_j;

                    // Assegna le coordinate di 'j' a 'i'
                    meshTriangulated.Cell0DsCoordinates.col(i) = meshTriangulated.Cell0DsCoordinates.col(j);
                }
            }
        }
    }

    // Applichiamo il remapping finale a tutti i vertici
    // Dobbiamo propagare le catene di reindirizzamento.
    for (unsigned int k = 0; k < n; ++k) {
        unsigned int current_id = k;
        while (id_remap[current_id] != current_id) {
            current_id = id_remap[current_id];
        }
        id_remap[k] = current_id; // Imposta il reindirizzamento finale per k
    }

    // Aggiorniamo le strutture della mesh in base al remapping finale
    for (int edgeId = 0; edgeId < meshTriangulated.Cell1DsExtrema.rows(); ++edgeId) {
		meshTriangulated.Cell1DsExtrema(edgeId, 0) = id_remap[meshTriangulated.Cell1DsExtrema(edgeId, 0)];
		meshTriangulated.Cell1DsExtrema(edgeId, 1) = id_remap[meshTriangulated.Cell1DsExtrema(edgeId, 1)];
	}
    
    for (unsigned int faceId = 0; faceId < meshTriangulated.Cell2DsId.size(); ++faceId) {
        for (unsigned int& vertexOriginalId : meshTriangulated.Cell2DsVertices[faceId]) {
            vertexOriginalId = id_remap[vertexOriginalId];
        }
    }

    // Aggiorniamo Cell0DsId e Cell0DsFlag in base al remapping.
    // I vertici che sono "master" avranno il loro id_remap[k] == k.
    // I vertici che sono duplicati avranno id_remap[k] != k.
    for (unsigned int k = 0; k < n; ++k) {
        meshTriangulated.Cell0DsId[k] = id_remap[k];
        if (id_remap[k] == k) {
            // Questo vertice è un master o non è stato reindirizzato
            meshTriangulated.Cell0DsFlag[k] = {maxFlag};
        } 
    }
}

    // Questa funzione individua e gestisce i lati duplicati all'interno di una mesh poliedrica triangolata.
    // Due lati (Cell1Ds) sono considerati duplicati se:
	//   1. Hanno la stessa "flag" (ossia appartengono allo stesso tipo o gruppo topologico).
	//   2. I loro vertici estremi coincidono geometricamente, anche in ordine inverso, entro una piccola tolleranza (tol).
	//
	// Invece di eliminare fisicamente i lati duplicati, la funzione costruisce un vettore di reindirizzamento (`id_remap`)
	// in cui, per ogni lato duplicato, si indica quale sia il lato "master" da utilizzare al suo posto.
	//
	// La funzione:
	// - Analizza ogni coppia di lati per verificare se sono duplicati.
	// - Marca come duplicati quelli ridondanti (non master) e aggiorna `id_remap` per riflettere il reindirizzamento.
	//
	// Il risultato è una mesh in cui i lati duplicati sono logicamente unificati

	void RimuoviLatiDuplicati(PolyhedralMesh& meshTriangulated)
{
    double tol = 1e-12;
    unsigned int maxFlag = numeric_limits<unsigned int>::max();
    size_t n = meshTriangulated.Cell1DsExtrema.rows();
    
    if (meshTriangulated.Cell1DsOriginalFlag.empty() || meshTriangulated.Cell1DsOriginalFlag.size() != n) {
        meshTriangulated.Cell1DsOriginalFlag.resize(n);
    }
    
    for (size_t k = 0; k < n; ++k) {
        meshTriangulated.Cell1DsOriginalFlag[k] = (meshTriangulated.Cell1DsFlag[k] == maxFlag);
    }

    // Vettore per tenere traccia del reindirizzamento finale degli ID dei lati.
    // Inizialmente, ogni lato punta a se stesso.
    vector<unsigned int> id_remap(n);
    for (unsigned int k = 0; k < n; ++k) {
        id_remap[k] = k;
    }

    // Vettore che tiene traccia di quali lati sono ancora candidati a essere non duplicati. 
    // Se un lato viene identificato come duplicato di un altro, la sua flag is_master_candidate viene impostata a false
    vector<bool> is_master_candidate(n, true);

    for (size_t i = 0; i < n; ++i) {
        // Se il lato corrente 'i' è già stato marcato per essere tenuto (maxFlag)
        // o se è già stato reindirizzato da un lato precedente nel ciclo, salta.
        if (meshTriangulated.Cell1DsFlag[i] == maxFlag || !is_master_candidate[i])
            continue;

        for (size_t j = i + 1; j < n; ++j) {
            if (meshTriangulated.Cell1DsFlag[j] == maxFlag || !is_master_candidate[j])
                continue;

            if (meshTriangulated.Cell1DsFlag[i] == meshTriangulated.Cell1DsFlag[j]) {
                int i0 = meshTriangulated.Cell1DsExtrema(i, 0); // Primo estremo del lato i
                int i1 = meshTriangulated.Cell1DsExtrema(i, 1); // Secondo estremo del lato i
                int j0 = meshTriangulated.Cell1DsExtrema(j, 0); // Primo estremo del lato j
                int j1 = meshTriangulated.Cell1DsExtrema(j, 1); // Secondo estremo del lato j

                // Controlliamo se i vertici estremi corrispondono in ordine diretto o inverso
                bool match_direct = ((meshTriangulated.Cell0DsCoordinates.col(i0) - meshTriangulated.Cell0DsCoordinates.col(j0)).norm() < tol &&
                                     (meshTriangulated.Cell0DsCoordinates.col(i1) - meshTriangulated.Cell0DsCoordinates.col(j1)).norm() < tol);

                bool match_inverse = ((meshTriangulated.Cell0DsCoordinates.col(i0) - meshTriangulated.Cell0DsCoordinates.col(j1)).norm() < tol &&
                                      (meshTriangulated.Cell0DsCoordinates.col(i1) - meshTriangulated.Cell0DsCoordinates.col(j0)).norm() < tol);

                if (match_direct || match_inverse) {
                    //  Il lato 'i' è un duplicato del lato 'j', vogliamo mantenere 'j' e reindirizzare 'i' a 'j'.

                    is_master_candidate[i] = false; // Marca 'i' duplicato.

                    unsigned int root_i = i;
                    while (id_remap[root_i] != root_i) {
                        root_i = id_remap[root_i];
                    }
                    unsigned int root_j = j;
                    while (id_remap[root_j] != root_j) {
                        root_j = id_remap[root_j];
                    }

                    if (root_i != root_j) {
                         id_remap[root_i] = root_j;
                    }
                    // Assicurati che il lato 'i' stesso (e i suoi precedenti duplicati)
                    // puntino a 'j' o al suo master.
                    id_remap[i] = root_j; // Imposta il reindirizzamento diretto per 'i'

                    //meshTriangulated.Cell1DsExtrema.row(i) = meshTriangulated.Cell1DsExtrema.row(j);
                }
            }
        }
    }

    for (unsigned int k = 0; k < n; ++k) {
        unsigned int current_id = k;
        while (id_remap[current_id] != current_id) {
            current_id = id_remap[current_id];
            id_remap[k] = current_id;
        }
        id_remap[k] = current_id;
    }

    // Aggiorniamo le strutture della mesh in base al remapping finale
    // Aggiorna i riferimenti ai lati nelle facce (Cell2DsEdges)
    for (unsigned int faceId = 0; faceId < meshTriangulated.Cell2DsId.size(); ++faceId) {
        for (unsigned int& edgeOriginalId : meshTriangulated.Cell2DsEdges[faceId]) {
            edgeOriginalId = id_remap[edgeOriginalId];
        }
    }

    // Aggiorniamo Cell1DsId e Cell1DsFlag in base al remapping.
    for (unsigned int k = 0; k < n; ++k) {
        meshTriangulated.Cell1DsId[k] = id_remap[k];
        if (id_remap[k] == k) {
            meshTriangulated.Cell1DsFlag[k] = maxFlag;
        } else {
            // Questo lato è un duplicato e punterà a un master
        }
    }
}
	
	void NewMesh(const PolyhedralMesh& meshTriangulated, PolyhedralMesh& meshFinal, const vector<int>& dimension)
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