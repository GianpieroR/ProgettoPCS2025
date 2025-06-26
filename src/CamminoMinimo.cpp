#include "Utils.hpp"          
#include "PolyhedralMesh.hpp" 
#include <iostream>           // Per std::cout, std::cerr
#include <vector>             // Per std::vector
#include <queue>              // Per std::priority_queue
#include <map>                // Per std::map
#include <limits>             // Per std::numeric_limits
#include <cmath>              // Per std::sqrt (usato in calculateDistanceById)
#include <algorithm>          // Per std::min, std::max
#include <tuple>              // Per std::tie o destructuring assignment (C++17)
#include <utility>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen; 
namespace PolyhedralLibrary{

// IMPLEMENTAZIONE DEL COSTRUTTORE DI SHORTESTPATHRESULT
ShortestPathResult::ShortestPathResult(unsigned int nEdges, double len)
    : numEdges(nEdges), totalLength(len)
{}

// Funzione per calcolare la distanza euclidea tra due punti.
double calculateDistanceById(const PolyhedralMesh& mesh, const unsigned int id1, const unsigned int id2) {
    
    VectorXd p1 = mesh.Cell0DsCoordinates.col(id1);
    VectorXd p2 = mesh.Cell0DsCoordinates.col(id2);

    return (p1 - p2).norm(); // Calcola la norma (distanza euclidea)
}


// Funzione per calcolare la matrice di adiacenza
MatrixXi calculateAdjacencyMatrix(const PolyhedralMesh& mesh) {
    const unsigned int numVertices = mesh.Cell0DsCoordinates.cols();
	MatrixXi adjMatrix = MatrixXi::Zero(numVertices, numVertices);

    // Itera su tutti i lati (Cell1Ds) della mesh
    for (unsigned int i = 0; i < mesh.Cell1DsId.size(); ++i) {

        unsigned int v1 = mesh.Cell1DsExtrema(i, 0);
        unsigned int v2 = mesh.Cell1DsExtrema(i, 1);
		
        adjMatrix(v1, v2) = 1;
        adjMatrix(v2, v1) = 1;
    }

    return adjMatrix;
}

ShortestPathResult findShortestPathDijkstra(
    PolyhedralMesh& mesh,
    const MatrixXi& adjMatrix,
    const unsigned int startVertexId,
    const unsigned int endVertexId
) {
	const unsigned int numVertices = mesh.Cell0DsCoordinates.cols(); // Numero totale di vertici
    const unsigned int numEdgesInMesh = mesh.Cell1DsId.size();     // Numero totale di lati nella mesh
    
    if (startVertexId >= numVertices) {
        cerr << "Errore: startVertexId (" << startVertexId << ") è fuori dal range valido di vertici [0, " << numVertices - 1 << "]." << endl;
        // Restituisci un risultato vuoto per indicare un errore
        return ShortestPathResult(0, 0.0);
    }
    if (endVertexId >= numVertices) {
        cerr << "Errore: endVertexId (" << endVertexId << ") è fuori dal range valido di vertici [0, " << numVertices - 1 << "]." << endl;
        // Restituisci un risultato vuoto per indicare un errore
        return ShortestPathResult(0, 0.0);
    }

	// Inizializza il risultato
    ShortestPathResult result(0, 0.0);

    // Caso banale: partenza e arrivo sono lo stesso vertice
    if (startVertexId == endVertexId) {
        cout << "Partenza e arrivo coincidono. Cammino nullo." << endl;
        mesh.Cell0DsMarker.assign(numVertices, 0); 
        mesh.Cell1DsMarker.assign(numEdgesInMesh, 0);
        mesh.Cell0DsMarker[startVertexId] = 1; // Marca il vertice sulla mesh
        return result;
    }

    // Mappa per collegare una coppia di vertici (tramite INDICI)
    // all'ID del lato e alla sua lunghezza.
    map<pair<unsigned int, unsigned int>, pair<unsigned int, double>> edgeInfoMap;
    
	for (unsigned int i = 0; i < numEdgesInMesh; ++i) {
        unsigned int v1 = mesh.Cell1DsExtrema(i, 0);
        unsigned int v2 = mesh.Cell1DsExtrema(i, 1);

        double length = calculateDistanceById(mesh, v1, v2);

        // Memorizziamo l'informazione usando gli INDICI dei vertici, garantendo ordine per la chiave della mappa
        edgeInfoMap[{min(v1, v2), max(v1, v2)}] = {mesh.Cell1DsId[i], length};
    }

    // Variabili Dijkstra
	
    // dist[i]: distanza minima conosciuta dal vertice di partenza a i
    vector<double> dist(numVertices, numeric_limits<double>::infinity());
	
    // predVertex[i]: indice del vertice precedente nel cammino minimo a i
    vector<unsigned int> predVertex(numVertices, -1); // Usiamo -1 per indicare nessun predecessore
	
    // predEdge[i]: ID del Cell1D usato per raggiungere i dal suo predecessore
    vector<unsigned int> predEdge(numVertices, -1);

    // visited[i]: true se il cammino più breve a 'i' è stato finalizzato
    vector<bool> visited(numVertices, false); 

    // Coda di priorità: memorizza coppie {distanza, indice_vertice}
    // std::greater per creare un min-heap (estrae l'elemento con la distanza minore)
    using QueueElem = pair<double, unsigned int>;
    priority_queue<QueueElem, vector<QueueElem>, greater<>> pq;
	
	// priority_queue<QueueElem, vector<QueueElem>, greater<QueueElem>> pq;

    // Imposta la distanza del vertice di partenza a 0 e aggiungilo alla coda
    dist[startVertexId] = 0.0;
    pq.push({0.0, startVertexId});

    // algoritmo Dijkstra
    while (!pq.empty()) {
        // Estrai il vertice 'u' con la distanza minima corrente dalla coda
		auto [current_distance, u] = pq.top(); // accedo all'elemento
        pq.pop(); // rimuovo l'elemento

        // Se il vertice è già stato visitato, significa che abbiamo già trovato
        // il cammino più breve per esso, quindi ignoriamo questa istanza.
        if (visited[u]) {
            continue;
        }
        visited[u] = true; // Marca il vertice come visitato/finalizzato

        // Se abbiamo raggiunto il vertice di arrivo, possiamo terminare l'algoritmo
        if (u == endVertexId) {
            break;
        }

        // Itera su tutti i possibili vicini 'v' di 'u' usando la matrice di adiacenza
        for (unsigned int v = 0; v < numVertices; ++v) {
            // Se c'è un lato tra u e v (adjMatrix(u, v) == 1)
            if (adjMatrix(u, v) == 1) {
                // Recupera le informazioni sul lato (ID reale e lunghezza/peso)
                auto it_edge = edgeInfoMap.find({min(u, v), max(u, v)});
                if (it_edge == edgeInfoMap.end()) {
                    cerr << "Avviso: Lato tra indici (" << u << "," << v << ") presente in AdjacencyMatrix ma non trovato in edgeInfoMap.\n";
                    continue;
                }
                double weight = it_edge->second.second; // Peso del lato (lunghezza euclidea)

                // Operazione di "rilassamento":
                // Se la distanza calcolata a 'v' passando per 'u' è minore della distanza attuale di 'v'
                if (dist[u] + weight < dist[v]) {
                    dist[v] = dist[u] + weight; // Aggiorna la distanza minima per 'v'
                    predVertex[v] = u;          // Imposta 'u' come predecessore di 'v'
                    predEdge[v] = it_edge->second.first; // Memorizza l'ID reale del lato usato
                    pq.push({dist[v], v});      // Inserisci 'v' nella coda di priorità con la nuova distanza
                }
            }
        }
    }
	
    // Se la distanza al vertice di arrivo è ancora infinito, significa che non è stato trovato alcun cammino
    if (dist[endVertexId] == numeric_limits<double>::infinity()) {
        cout << "Nessun cammino trovato tra il vertice " << startVertexId
                  << " e il vertice " << endVertexId << endl;
        mesh.Cell0DsMarker.assign(numVertices, 0);
        mesh.Cell1DsMarker.assign(numEdgesInMesh, 0);
		return result;
    }
	
	
    // Ricostruzione del cammino e calcolo delle statistiche
    // Ricostruisci il cammino a ritroso dal vertice di arrivo al vertice di partenza
	
	// Inizializzo i marker della mesh a 0
    mesh.Cell0DsMarker.assign(numVertices, 0); // Inizializza i marker a 0
    mesh.Cell1DsMarker.assign(numEdgesInMesh, 0);
	
	result.totalLength = dist[endVertexId]; // distanza totale

    unsigned int current_idx = endVertexId; // Partiamo dall'indice del vertice di arrivo
	
	while (current_idx != startVertexId) {
        // Marca il vertice corrente come parte del cammino minimo
        mesh.Cell0DsMarker[current_idx] = 1;

        unsigned int prev_vertex_idx = predVertex[current_idx]; // Indice del vertice precedente nel cammino
        unsigned int edge_used_id = predEdge[current_idx];      // ID del lato usato per raggiungere current_idx

        mesh.Cell1DsMarker[edge_used_id] = 1;
        result.numEdges++;
        
        current_idx = prev_vertex_idx; // Spostati al vertice precedente e continua la ricostruzione
    }
	
    // Marca anche il vertice di partenza come parte del cammino
    mesh.Cell0DsMarker[startVertexId] = 1;	
	
	return result;
}

}