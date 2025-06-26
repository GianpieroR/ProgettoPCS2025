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
	// Funzione per trovare la forma normalizzata ciclicamente
	// cioè, ruotare la sequenza in modo che inizi con l'elemento più piccolo
	vector<unsigned int> get_cyclic_normalized(const vector<unsigned int>& current_edges) {
		if (current_edges.empty() || current_edges.size() == 1) {
			return current_edges; // Non c'è nulla da normalizzare per 0 o 1 elemento
		}
		vector<unsigned int> temp_edges = current_edges; // Lavora su una copia
		auto min_it = min_element(temp_edges.begin(), temp_edges.end());
		rotate(temp_edges.begin(), min_it, temp_edges.end());
		return temp_edges;
	}

	vector<unsigned int> NormalizeFaceEdges(const vector<unsigned int>& face_edges) {
		if (face_edges.empty() || face_edges.size() == 1) {
			return face_edges; // Non c'è nulla da normalizzare per 0 o 1 elemento
		}
	
		// Ottieni la forma normalizzata ciclicamente della sequenza originale
		vector<unsigned int> normalized_original = get_cyclic_normalized(face_edges);
	
		// Ottieni la sequenza inversa dell'originale
		vector<unsigned int> reversed_edges = face_edges;
		reverse(reversed_edges.begin(), reversed_edges.end()); // Inverte la copia
	
		// Ottieni la forma normalizzata ciclicamente anche della sequenza inversa
		vector<unsigned int> normalized_reversed = get_cyclic_normalized(reversed_edges);
	
		// Confronta le due forme normalizzate e restituisci la minore
		// confronta gli elementi uno a uno
		if (normalized_original < normalized_reversed) {
			return normalized_original;
		} else {
			return normalized_reversed;
		}
	}

	void FindAddFace(const vector<unsigned int>& new_face_vertices,
							   const vector<unsigned int>& new_face_edges,
							   PolyhedralMesh& meshTriangulated,
							   unsigned int& k3)
	{
		// Normalizza la sequenza degli spigoli
		vector<unsigned int> normalized_new_edges = NormalizeFaceEdges(new_face_edges);
	
		// Itero attraverso le facce esistenti per trovare un duplicato
		bool found = false;
		// mi serve perchè se no incrementerei k3 tutte le volte anche quando la faccia esiste già
		for (unsigned int i = 0; i < meshTriangulated.Cell2DsId.size(); ++i) {
			// Normalizza gli spigoli della faccia corrente per il confronto
			vector<unsigned int> normalized_existing_edges = NormalizeFaceEdges(meshTriangulated.Cell2DsEdges[i]);
	
			// Confronta le sequenze normalizzate. Se sono uguali, la faccia esiste già.
			if (normalized_new_edges == normalized_existing_edges) {
				found = true;
				break;
			}
		}
	
		// Se il ciclo è terminato e non siamo usciti dalla funzione,
		// significa che la faccia non è stata trovata. La aggiungo.
		if (!found){
			meshTriangulated.Cell2DsVertices[k3] = new_face_vertices;
			meshTriangulated.Cell2DsEdges[k3] = new_face_edges; // Salva la versione non normalizzata
			meshTriangulated.Cell2DsId[k3] = k3; // Assegna il nuovo ID
		
			k3++; // Restituisci il vecchio valore di k3, poi incrementalo per la prossima faccia
		}
	}

	
	unsigned int FindAddVertice(const Vector3d& coord, PolyhedralMesh& meshTriangulated, unsigned int& k1)
	{
		double tol = 1e-12;
		for (unsigned int i = 0; i < meshTriangulated.Cell0DsCoordinates.cols(); i++) {
			if ((meshTriangulated.Cell0DsCoordinates.col(i) - coord).norm() < tol) {
				return i;
			}
		}

		// se il ciclo è terminato e non siamo usciti dalla funzione vuole dire che il vertice non è stato trovato
		meshTriangulated.Cell0DsCoordinates.col(k1) = coord;
		meshTriangulated.Cell0DsId[k1] = k1;
		return k1++;
	}
	
	
	unsigned int FindAddEdge2(const unsigned int a, const unsigned int b, PolyhedralMesh& meshTriangulated, unsigned int& k2)
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
		meshTriangulated.Cell1DsExtrema(k2, 0) = a;
		meshTriangulated.Cell1DsExtrema(k2, 1) = b;
		meshTriangulated.Cell1DsId[k2] = k2;
		return k2++; // restituisci il vecchio valore di k2, poi incrementalo
	}
	
	
	// trova la faccia adiacente a face tramite edge, e poi restituire il baricentro di questa faccia adiacente
	Vector3d FindNearBarycenter(const PolyhedralMesh& meshTriangulated, const unsigned int edgeId, const unsigned int currentFaceId, map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap) {

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

			// Calcola e restituisci il baricentro della faccia adiacente.
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

        meshTriangulated.Cell0DsId.resize(dimension[0]);
        meshTriangulated.Cell0DsCoordinates = MatrixXd::Zero(3, dimension[0]);
        meshTriangulated.Cell0DsFlag.resize(dimension[0]);

        meshTriangulated.Cell1DsId.resize(dimension[1]);
        meshTriangulated.Cell1DsExtrema = MatrixXi::Zero(dimension[1], 2);
        meshTriangulated.Cell1DsFlag.resize(dimension[1]);

        meshTriangulated.Cell2DsId.resize(dimension[2]);
        meshTriangulated.Cell2DsVertices.resize(dimension[2]);
        meshTriangulated.Cell2DsEdges.resize(dimension[2]);

        unsigned int k1 = 0; // Contatore nodi globali (0D)
        unsigned int k2 = 0; // Contatore edge globali (1D)
        unsigned int k3 = 0; // Contatore facce/triangoli globali (2D)
        
        for (unsigned faceId = 0; faceId < mesh.Cell2DsId.size() ; ++faceId){
			const auto& faceVertices = mesh.Cell2DsVertices[faceId]; // id dei vertici della faccia originale
				
			Vector3d V0 = mesh.Cell0DsCoordinates.col(faceVertices[0]);
			Vector3d V1 = mesh.Cell0DsCoordinates.col(faceVertices[1]);
			Vector3d V2 = mesh.Cell0DsCoordinates.col(faceVertices[2]);
			
			Vector3d barycenter = ((V0 + V1 + V2)/3.0);
			
			const auto& faceEdges = mesh.Cell2DsEdges[faceId]; // id dei lati della faccia originale

			for (unsigned int e = 0; e <3; e++){
				
				if (!mesh.Cell1DsOriginalFlag[faceEdges[e]]) { // lato di bordo
					Vector3d mediumPoint = (mesh.Cell0DsCoordinates.col(faceVertices[e])+mesh.Cell0DsCoordinates.col(faceVertices[(e+1)%3]))/2.0;

					// TRIANGOLO SOPRA IL PUNTO MEDIO

					unsigned int vertex = FindAddVertice(mesh.Cell0DsCoordinates.col(faceVertices[e]), meshTriangulated, k1);
					unsigned int medium = FindAddVertice(mediumPoint, meshTriangulated, k1);
					unsigned int bar = FindAddVertice(barycenter, meshTriangulated, k1);
					
					vector<unsigned int> new_face_vertices = {vertex, medium, bar};
					
					// lato vertice - punto medio
					unsigned int edgeA = FindAddEdge2(vertex, medium, meshTriangulated, k2);
					
					// lato punto medio- baricentro
					unsigned int edgeB = FindAddEdge2(medium, bar, meshTriangulated, k2);
					
					// lato baricentro - vertice
					unsigned int edgeC = FindAddEdge2(bar, vertex, meshTriangulated, k2);
					
					vector<unsigned int> new_face_edges = {edgeA, edgeB, edgeC};

					FindAddFace(new_face_vertices, new_face_edges, meshTriangulated, k3);


					// TRIANGOLO SOTTO IL PUNTO MEDIO
					
					unsigned int vertex1 = FindAddVertice(mesh.Cell0DsCoordinates.col(faceVertices[(e+1)%3]), meshTriangulated, k1);
					vector<unsigned int> new_face_vertices1 = {vertex1, medium, bar};
					
					// lato vertice - punto medio
					unsigned int edgeD = FindAddEdge2(vertex1, medium, meshTriangulated, k2);
					
					// lato baricentro - vertice
					unsigned int edgeE = FindAddEdge2(bar, vertex1, meshTriangulated, k2);
					
					vector<unsigned int> new_face_edges1 = {edgeD, edgeB, edgeE};
					
					FindAddFace(new_face_vertices1, new_face_edges1, meshTriangulated, k3);
					
				
				} else {
					
				// TRIANGOLO A SINISTRA
				Vector3d bar2_coord = FindNearBarycenter(mesh, faceEdges[e], faceId, edgeToFacesMap);
				unsigned int vertex2 = FindAddVertice(mesh.Cell0DsCoordinates.col(faceVertices[e]), meshTriangulated, k1);
				unsigned int bar1 = FindAddVertice(barycenter, meshTriangulated, k1);
				unsigned int bar2 = FindAddVertice(bar2_coord, meshTriangulated, k1);
				
				vector<unsigned int> new_face_vertices2 = {vertex2, bar1, bar2};
				
				// lato vertice - bar1
				unsigned int edgeF = FindAddEdge2(vertex2, bar1, meshTriangulated, k2);
					
				// lato bar1- bar2
				unsigned int edgeG = FindAddEdge2(bar1, bar2, meshTriangulated, k2);
					
				// lato bar2 - vertice
				unsigned int edgeH = FindAddEdge2(bar2, vertex2, meshTriangulated, k2);
					
				vector<unsigned int> new_face_edges2 = {edgeF, edgeG, edgeH};
				
				FindAddFace(new_face_vertices2, new_face_edges2, meshTriangulated, k3);
					
				
				// TRIANGOLO A DESTRA
				unsigned int vertex3 = FindAddVertice(mesh.Cell0DsCoordinates.col(faceVertices[(e+1)%3]), meshTriangulated, k1);
				vector<unsigned int> new_face_vertices3 = {vertex3, bar1, bar2};
				
				// lato vertice - bar1
				unsigned int edgeI = FindAddEdge2(vertex3, bar1, meshTriangulated, k2);
					
				// lato bar2 - vertice
				unsigned int edgeL = FindAddEdge2(bar2, vertex3, meshTriangulated, k2);
					
				vector<unsigned int> new_face_edges3 = {edgeI, edgeG, edgeL};
				
				FindAddFace(new_face_vertices3, new_face_edges3, meshTriangulated, k3);	
				
				}
			}
		}
	}

}
