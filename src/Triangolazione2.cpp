#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <limits>
#include <algorithm>
#include <map>

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{
	// Funzione che prende una sequenza di indici (unsigned int) e la ruota ciclicamente
// in modo tale che l'elemento minimo della sequenza diventi il primo elemento.
// Questa operazione serve a "normalizzare" la sequenza, rendendo più facile il confronto
// tra sequenze cicliche che possono essere uguali ma con punti di partenza diversi.
vector<unsigned int> get_cyclic_normalized(const vector<unsigned int>& current_edges) {
    // Caso base: se la sequenza è vuota o ha un solo elemento,
    // non serve fare alcuna rotazione perché è già normalizzata
    if (current_edges.empty() || current_edges.size() == 1) {
        return current_edges;  // restituisce la sequenza così com'è
    }

    // Creiamo una copia della sequenza di input, così da non modificarla direttamente
    vector<unsigned int> temp_edges = current_edges;

    // Troviamo l'iteratore che punta al primo elemento minimo nella sequenza
    // std::min_element restituisce un iteratore al primo elemento che ha valore minimo
    auto min_it = min_element(temp_edges.begin(), temp_edges.end());

    // Ruotiamo la sequenza in modo ciclico utilizzando std::rotate:
    // spostiamo tutti gli elementi a sinistra dell'elemento minimo verso la fine della sequenza,
    // facendo sì che l'elemento minimo diventi il primo elemento della sequenza.
    // L'ordine relativo degli elementi rimane invariato.
    rotate(temp_edges.begin(), min_it, temp_edges.end());

    // Restituiamo la sequenza ruotata, ora "normalizzata" con l'elemento minimo in testa
    return temp_edges;
}


// Funzione per normalizzare la sequenza di indici degli spigoli di una faccia
// tenendo conto che la sequenza può essere letta sia in senso orario che antiorario.
// Restituisce la sequenza "normalizzata" che è la più piccola lessicograficamente
// tra la sequenza originale e la sua inversione.
vector<unsigned int> NormalizeFaceEdges(const vector<unsigned int>& face_edges) {
    // Se la sequenza è vuota o ha un solo elemento, non serve normalizzare
    if (face_edges.empty() || face_edges.size() == 1) {
        return face_edges;
    }

    // Ottieni la forma normalizzata ciclicamente della sequenza originale
    vector<unsigned int> normalized_original = get_cyclic_normalized(face_edges);

    // Crea una copia della sequenza originale per invertirla
    vector<unsigned int> reversed_edges = face_edges;

    // Inverte la sequenza (cioè legge la sequenza al contrario)
    reverse(reversed_edges.begin(), reversed_edges.end());

    // Ottieni la forma normalizzata ciclicamente anche della sequenza invertita
    vector<unsigned int> normalized_reversed = get_cyclic_normalized(reversed_edges);

    // Confronta lessicograficamente le due sequenze normalizzate
    // Restituisce quella che risulta "minore" secondo l'ordinamento lessicografico
    // Questo permette di trattare la sequenza e la sua inversione come equivalenti
    if (normalized_original < normalized_reversed) {
        return normalized_original;
    } else {
        return normalized_reversed;
    }
}

    // Questa funzione aggiunge una nuova faccia alla mesh triangolata (meshTriangulated), a condizione che non esista già. 
	// Riceve in input i vertici e gli spigoli (edges) che compongono la nuova faccia. Per evitare duplicati, normalizza gli spigoli della nuova faccia 
	// e li confronta con quelli già presenti nella mesh (anch'essi normalizzati con la funzione NormalizeFaceEdges). 
	// Se una faccia con gli stessi spigoli è già presente, non fa nulla; altrimenti, aggiunge la faccia nei vettori Cell2DsVertices, 
	// Cell2DsEdges e Cell2DsId della mesh, utilizzando l'indice y3 come identificatore univoco, e infine incrementa y3 per prepararlo alla prossima aggiunta. 
	// Questo è utile in fase di raffinamento o generazione di mesh, dove si possono generare molte facce simili e si vuole evitare di duplicarle.

	void FindAddFace(const vector<unsigned int>& new_face_vertices,
							   const vector<unsigned int>& new_face_edges,
							   PolyhedralMesh& meshTriangulated,
							   unsigned int& y3)
	{
		// Normalizza la sequenza degli spigoli
		vector<unsigned int> new_edges_normalized = NormalizeFaceEdges(new_face_edges);
	
		// Itero attraverso le facce esistenti per trovare un duplicato
		bool found = false;
		for (unsigned int i = 0; i < meshTriangulated.Cell2DsId.size(); ++i) {
			// Normalizza gli spigoli della faccia corrente per il confronto
			vector<unsigned int> normalized_existing_edges = NormalizeFaceEdges(meshTriangulated.Cell2DsEdges[i]);
	
			// Confronta le sequenze normalizzate. Se sono uguali, la faccia esiste già.
			if (new_edges_normalized == normalized_existing_edges) {
				found = true;
				break;
			}
		}
	
		// Se il ciclo è terminato e non siamo usciti dalla funzione vuol dire che la faccia non è stata trovata e la devo aggiungere.
		if (!found){
			meshTriangulated.Cell2DsVertices[y3] = new_face_vertices;
			meshTriangulated.Cell2DsEdges[y3] = new_face_edges; // Salva la versione non normalizzata
			meshTriangulated.Cell2DsId[y3] = y3; // Assegna il nuovo ID
		
			y3++; // Restituisci il vecchio valore di y3, poi incrementalo per la prossima faccia
		}
	}

	
	unsigned int FindAddVertice(const Vector3d& coord, PolyhedralMesh& meshTriangulated, unsigned int& y1)
	{
		double tol = 1e-12;
		for (unsigned int i = 0; i < meshTriangulated.Cell0DsCoordinates.cols(); i++) {
			if ((meshTriangulated.Cell0DsCoordinates.col(i) - coord).norm() < tol) {
				return i;
			}
		}

		// se il ciclo è terminato e non siamo usciti dalla funzione vuole dire che il vertice non è stato trovato
		meshTriangulated.Cell0DsCoordinates.col(y1) = coord;
		meshTriangulated.Cell0DsId[y1] = y1;
		return y1++;
	}
	
	
	unsigned int FindAddEdge2(const unsigned int a, const unsigned int b, PolyhedralMesh& meshTriangulated, unsigned int& y2)
	{
		// Controllo se l'edge (a,b) o (b,a) esiste già
		for (unsigned int i = 0; i < meshTriangulated.Cell1DsExtrema.rows(); ++i) {
			if ((meshTriangulated.Cell1DsExtrema(i, 0) == static_cast<int>(a) && meshTriangulated.Cell1DsExtrema(i, 1) == static_cast<int>(b)) ||
    			(meshTriangulated.Cell1DsExtrema(i, 0) == static_cast<int>(b) && meshTriangulated.Cell1DsExtrema(i, 1) == static_cast<int>(a))) {

				// Edge già presente ⇒ associarlo alla faccia
				return i; // restituisci l'ID dell'edge
			}
		}

		// Edge non presente ⇒ lo aggiungo
		meshTriangulated.Cell1DsExtrema(y2, 0) = a;
		meshTriangulated.Cell1DsExtrema(y2, 1) = b;
		meshTriangulated.Cell1DsId[y2] = y2;
		return y2++; // restituisci il vecchio valore di y2, poi incrementalo
	}
	
	
	// Questa funzione, dato un edge (edgeId) e una faccia corrente (currentFaceId), restituisce il baricentro della faccia adiacente 
	// che condivide lo stesso edge. Utilizza una mappa edgeToFacesMap che associa a ogni coppia ordinata di vertici (cioè un lato) 
	// le facce che lo contengono. Se l’edge è condiviso da due facce, viene identificata quella diversa da quella corrente, e viene 
	// calcolato e restituito il suo baricentro. Se invece l’edge è usato da una sola faccia (cioè è un bordo), viene restituito un vettore nullo.
	
	Vector3d FindNearbaricentroycenter(const PolyhedralMesh& meshTriangulated, const unsigned int edgeId, const unsigned int currentFaceId, map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap) {

		// Ottieni i vertici che compongono il lato
		unsigned int v1_id = meshTriangulated.Cell1DsExtrema(edgeId, 0);
		unsigned int v2_id = meshTriangulated.Cell1DsExtrema(edgeId, 1);
		pair<unsigned int, unsigned int> edgeKey = {min(v1_id, v2_id), max(v1_id, v2_id)};

		// Trova le facce associate a questo lato
		auto it = edgeToFacesMap.find(edgeKey);
		if (it == edgeToFacesMap.end()) {
			cerr << "Errore: Edge " << edgeId << " non trovato nella mappa delle facce adiacenti." << endl;
			return Vector3d::Zero(); // Nessuna faccia associata a questo edge
		}

		const vector<unsigned int>& facesSharingEdge = it->second;

		if (facesSharingEdge.size() == 2) {
			// Ci sono due facce. Una è currentFaceId, l'altra è quella che cerchiamo.
			unsigned int faceId1 = facesSharingEdge[0];
			unsigned int faceId2 = facesSharingEdge[1];

			unsigned int targetFaceId;
			if (faceId1 == currentFaceId) {
				targetFaceId = faceId2;
			} else if (faceId2 == currentFaceId) {
				targetFaceId = faceId1;
			} else {
				// Questo caso significa che currentFaceId non è una delle facce che condividono l'edge dato.
				cerr << "Errore: Faccia corrente " << currentFaceId << " non condivide l'edge " << edgeId << "." << endl;
				return Vector3d::Zero();
			}

			// Calcola e restituisci il baricentroicentro della faccia adiacente.
			return getFaceBarycenter(meshTriangulated, targetFaceId);

		} else if (facesSharingEdge.size() == 1) {
			cerr << "Warning: Edge " << edgeId << " è un bordo della mesh. Nessuna faccia adiacente trovata." << endl;
			return Vector3d::Zero(); 
		} else {
			cerr << "Errore: Edge " << edgeId << " è condiviso da " << facesSharingEdge.size() << " facce" << endl;
			return Vector3d::Zero();
		}
}

    void triangulateAndStore2(PolyhedralMesh& mesh, PolyhedralMesh& meshTriangulated, const vector<int>& dimension, map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap) {

    // Alloca spazio per i vertici (0D) nella mesh triangolata
    // - dimension[0] è il numero totale di vertici previsti dopo la triangolazione
    meshTriangulated.Cell0DsId.resize(dimension[0]);                             // ID dei vertici
    meshTriangulated.Cell0DsCoordinates = MatrixXd::Zero(3, dimension[0]);      // Coordinate dei vertici (3D)
    meshTriangulated.Cell0DsFlag.resize(dimension[0]);                          // Flag associati ai vertici (es. edge di provenienza)

    // Alloca spazio per gli spigoli (1D)
    // - dimension[1] è il numero totale di spigoli previsti
    meshTriangulated.Cell1DsId.resize(dimension[1]);                            // ID degli spigoli
    meshTriangulated.Cell1DsExtrema = MatrixXi::Zero(dimension[1], 2);         // Estremi di ogni spigolo (coppie di vertici)
    meshTriangulated.Cell1DsFlag.resize(dimension[1]);                          // Flag associati agli spigoli

    // Alloca spazio per le facce (2D)
    // - dimension[2] è il numero totale di facce triangolari previste
    meshTriangulated.Cell2DsId.resize(dimension[2]);                            // ID delle facce
    meshTriangulated.Cell2DsVertices.resize(dimension[2]);                      // Elenco dei vertici di ciascuna faccia (3 per triangoli)
    meshTriangulated.Cell2DsEdges.resize(dimension[2]);                         // Elenco degli spigoli associati a ciascuna faccia

        unsigned int y1 = 0; // Contatore nodi globali (0D)
        unsigned int y2 = 0; // Contatore edge globali (1D)
        unsigned int y3 = 0; // Contatore facce/triangoli globali (2D)
        
        for (unsigned faceId = 0; faceId < mesh.Cell2DsId.size() ; ++faceId){
			const auto& faceVertices = mesh.Cell2DsVertices[faceId]; // id dei vertici della faccia originale
				
			Vector3d V0 = mesh.Cell0DsCoordinates.col(faceVertices[0]);
			Vector3d V1 = mesh.Cell0DsCoordinates.col(faceVertices[1]);
			Vector3d V2 = mesh.Cell0DsCoordinates.col(faceVertices[2]);
			
			Vector3d baricentroycenter = ((V0 + V1 + V2)/3.0);
			
			const auto& faceEdges = mesh.Cell2DsEdges[faceId]; // id dei lati della faccia originale

			for (unsigned int e = 0; e <3; e++){
				
				if (!mesh.Cell1DsOriginalFlag[faceEdges[e]]) { // lato di bordo
					Vector3d medioPoint = (mesh.Cell0DsCoordinates.col(faceVertices[e])+mesh.Cell0DsCoordinates.col(faceVertices[(e+1)%3]))/2.0;

					// TRIANGOLO SOPRA IL PUNTO MEDIO

					unsigned int vertice = FindAddVertice(mesh.Cell0DsCoordinates.col(faceVertices[e]), meshTriangulated, y1);
					unsigned int medio = FindAddVertice(medioPoint, meshTriangulated, y1);
					unsigned int baricentro = FindAddVertice(baricentroycenter, meshTriangulated, y1);
					
					vector<unsigned int> new_face_vertices = {vertice, medio, baricentro};
					
					// lato vertice - punto medio
					unsigned int latoA = FindAddEdge2(vertice, medio, meshTriangulated, y2);
					
					// lato punto medio- baricentroicentro
					unsigned int latoB = FindAddEdge2(medio, baricentro, meshTriangulated, y2);
					
					// lato baricentroicentro - vertice
					unsigned int latoC = FindAddEdge2(baricentro, vertice, meshTriangulated, y2);
					
					vector<unsigned int> new_face_edges = {latoA, latoB, latoC};

					FindAddFace(new_face_vertices, new_face_edges, meshTriangulated, y3);


					// TRIANGOLO SOTTO IL PUNTO MEDIO
					
					unsigned int vertice1 = FindAddVertice(mesh.Cell0DsCoordinates.col(faceVertices[(e+1)%3]), meshTriangulated, y1);
					vector<unsigned int> new_face_vertices1 = {vertice1, medio, baricentro};
					
					// lato vertice - punto medio
					unsigned int edgeD = FindAddEdge2(vertice1, medio, meshTriangulated, y2);
					
					// lato baricentroicentro - vertice
					unsigned int edgeE = FindAddEdge2(baricentro, vertice1, meshTriangulated, y2);
					
					vector<unsigned int> new_face_edges1 = {edgeD, latoB, edgeE};
					
					FindAddFace(new_face_vertices1, new_face_edges1, meshTriangulated, y3);
					
				
				} else {
					
				// TRIANGOLO A SINISTRA
				// Trova il punto baricentrico della faccia adiacente (vicino al lato corrente)
				Vector3d baricentroicentro2_coord = FindNearbaricentroycenter(mesh, faceEdges[e], faceId, edgeToFacesMap);
				
				// Aggiunge (o recupera) il vertice corrente della faccia originale nella mesh triangolata
				unsigned int vertice2 = FindAddVertice(mesh.Cell0DsCoordinates.col(faceVertices[e]), meshTriangulated, y1);
				
				// Aggiunge il baricentro della faccia originale
				unsigned int baricentroicentro1 = FindAddVertice(baricentroycenter, meshTriangulated, y1);
				
				// Aggiunge il punto baricentrico dell'altra faccia (vicina lungo l'edge)
				unsigned int baricentroicentro2 = FindAddVertice(baricentroicentro2_coord, meshTriangulated, y1);
				
				// Definisce i vertici del nuovo triangolo sinistro
				vector<unsigned int> new_face_vertices2 = {vertice2, baricentroicentro1, baricentroicentro2};
				
				// Aggiunge gli spigoli tra i vertici del triangolo
				unsigned int edgeF = FindAddEdge2(vertice2, baricentroicentro1, meshTriangulated, y2);
				unsigned int edgeG = FindAddEdge2(baricentroicentro1, baricentroicentro2, meshTriangulated, y2);
				unsigned int edgeH = FindAddEdge2(baricentroicentro2, vertice2, meshTriangulated, y2);
					
					// Definisce gli spigoli del triangolo
				vector<unsigned int> new_face_edges2 = {edgeF, edgeG, edgeH};
				
				// Aggiunge la nuova faccia triangolare alla mesh
				FindAddFace(new_face_vertices2, new_face_edges2, meshTriangulated, y3);
					
				
				// TRIANGOLO A DESTRA
				
				// Aggiunge (o recupera) il vertice successivo della faccia originale
				unsigned int vertice3 = FindAddVertice(mesh.Cell0DsCoordinates.col(faceVertices[(e+1)%3]), meshTriangulated, y1);
				
				// Definisce i vertici del triangolo destro
				vector<unsigned int> new_face_vertices3 = {vertice3, baricentroicentro1, baricentroicentro2};
				
				// Aggiunge i bordi tra i vertici
				unsigned int edgeI = FindAddEdge2(vertice3, baricentroicentro1, meshTriangulated, y2);
				unsigned int edgeL = FindAddEdge2(baricentroicentro2, vertice3, meshTriangulated, y2);
					
				// I bordi sono gli stessi del triangolo sinistro, riutilizzati dove possibile
				vector<unsigned int> new_face_edges3 = {edgeI, edgeG, edgeL};
				
				// Aggiunge anche questa faccia triangolare alla mesh
				FindAddFace(new_face_vertices3, new_face_edges3, meshTriangulated, y3);	
				
				}
			}
		}
	}

}
