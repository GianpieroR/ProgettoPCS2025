#include "Utils.hpp"          
#include "PolyhedralMesh.hpp" 
#include <iostream>           
#include <vector>             
#include <queue>              
#include <map>                
#include <limits>             
#include <cmath>              
#include <algorithm>          
#include <tuple>              
#include <utility>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen; 

namespace PolyhedralLibrary{

// COSTRUTTORE DEL RISULTATO DEL CAMMINO MINIMO
ShortestPathResult::ShortestPathResult(unsigned int edgesCount, double pathLength)
    : numEdges(edgesCount), totalLength(pathLength)
{}

// Calcola la distanza euclidea tra due vertici dati i loro ID
double calculateDistanceById(const PolyhedralMesh& mesh, const unsigned int vertexId1, const unsigned int vertexId2) {
    VectorXd point1 = mesh.Cell0DsCoordinates.col(vertexId1);
    VectorXd point2 = mesh.Cell0DsCoordinates.col(vertexId2);
    return (point1 - point2).norm();  // Norma vettoriale = distanza euclidea
}

// Costruisce la matrice di adiacenza della mesh
MatrixXi calculateAdjacencyMatrix(const PolyhedralMesh& mesh) {
    const unsigned int totalVertices = mesh.Cell0DsCoordinates.cols();
    MatrixXi adjacencyMatrix = MatrixXi::Zero(totalVertices, totalVertices);

    // Per ogni lato nella mesh, indica connessione bidirezionale tra i due vertici
    for (unsigned int edgeIdx = 0; edgeIdx < mesh.Cell1DsId.size(); ++edgeIdx) {
        unsigned int vertexA = mesh.Cell1DsExtrema(edgeIdx, 0);
        unsigned int vertexB = mesh.Cell1DsExtrema(edgeIdx, 1);
        adjacencyMatrix(vertexA, vertexB) = 1;
        adjacencyMatrix(vertexB, vertexA) = 1;
    }

    return adjacencyMatrix;
}

// Implementazione dell'algoritmo di Dijkstra per il cammino minimo tra due vertici
ShortestPathResult findShortestPathDijkstra(
    PolyhedralMesh& mesh,
    const MatrixXi& adjacencyMatrix,
    const unsigned int startVertex,
    const unsigned int targetVertex
) {
    const unsigned int totalVertices = mesh.Cell0DsCoordinates.cols();
    const unsigned int totalEdges = mesh.Cell1DsId.size();

    // Controllo validità degli ID di partenza e arrivo
    if (startVertex >= totalVertices) {
        cerr << "Errore: startVertex (" << startVertex << ") non valido." << endl;
        return ShortestPathResult(0, 0.0);
    }
    if (targetVertex >= totalVertices) {
        cerr << "Errore: targetVertex (" << targetVertex << ") non valido." << endl;
        return ShortestPathResult(0, 0.0);
    }

    // Risultato iniziale (cammino vuoto)
    ShortestPathResult pathResult(0, 0.0);

    // Caso semplice: partenza e arrivo coincidono
    if (startVertex == targetVertex) {
        cout << "Vertice di partenza e arrivo coincidono. Cammino nullo." << endl;
        mesh.Cell0DsMarker.assign(totalVertices, 0);
        mesh.Cell1DsMarker.assign(totalEdges, 0);
        mesh.Cell0DsMarker[startVertex] = 1; // Marca solo il vertice di partenza
        return pathResult;
    }

    // Mappa che associa ogni coppia ordinata di vertici all'ID del lato e alla sua lunghezza
    map<pair<unsigned int, unsigned int>, pair<unsigned int, double>> edgeDataMap;
    for (unsigned int edgeIdx = 0; edgeIdx < totalEdges; ++edgeIdx) {
        unsigned int vertex1 = mesh.Cell1DsExtrema(edgeIdx, 0);
        unsigned int vertex2 = mesh.Cell1DsExtrema(edgeIdx, 1);
        double edgeLength = calculateDistanceById(mesh, vertex1, vertex2);

        // Salva le informazioni usando coppie ordinate (min, max) come chiave
        edgeDataMap[{min(vertex1, vertex2), max(vertex1, vertex2)}] = {mesh.Cell1DsId[edgeIdx], edgeLength};
    }

    // Vettori per l'algoritmo di Dijkstra
    vector<double> minDistance(totalVertices, numeric_limits<double>::infinity()); // Distanza minima da startVertex
    vector<unsigned int> previousVertex(totalVertices, -1);                        // Predecessore del vertice nel cammino minimo
    vector<unsigned int> previousEdge(totalVertices, -1);                          // ID del lato usato per arrivare al vertice
    vector<bool> isProcessed(totalVertices, false);                                // Indica se il vertice è stato processato definitivamente

    // Coda di priorità per vertici da esplorare: (distanza, vertice)
    using QueueElement = pair<double, unsigned int>;
    priority_queue<QueueElement, vector<QueueElement>, greater<>> vertexQueue;

    // Inizializza vertice di partenza
    minDistance[startVertex] = 0.0;
    vertexQueue.push({0.0, startVertex});

    // Ciclo principale di Dijkstra
    while (!vertexQueue.empty()) {
        auto [currentDist, currentVertex] = vertexQueue.top();
        vertexQueue.pop();

        // Se il vertice è già stato processato, saltalo
        if (isProcessed[currentVertex]) continue;

        isProcessed[currentVertex] = true;

        // Se siamo arrivati al vertice di destinazione, usciamo
        if (currentVertex == targetVertex) break;

        // Esamina tutti i vicini del vertice corrente
        for (unsigned int neighbor = 0; neighbor < totalVertices; ++neighbor) {
            if (adjacencyMatrix(currentVertex, neighbor) == 1) {
                // Recupera info del lato tra currentVertex e neighbor
                auto edgeIt = edgeDataMap.find({min(currentVertex, neighbor), max(currentVertex, neighbor)});
                if (edgeIt == edgeDataMap.end()) {
                    cerr << "Avviso: lato (" << currentVertex << "," << neighbor << ") non trovato in edgeDataMap.\n";
                    continue;
                }

                double edgeWeight = edgeIt->second.second;
                // Se passando per currentVertex si migliora la distanza a neighbor, aggiorna
                if (minDistance[currentVertex] + edgeWeight < minDistance[neighbor]) {
                    minDistance[neighbor] = minDistance[currentVertex] + edgeWeight;
                    previousVertex[neighbor] = currentVertex;
                    previousEdge[neighbor] = edgeIt->second.first;
                    vertexQueue.push({minDistance[neighbor], neighbor});
                }
            }
        }
    }

    // Se non è stato trovato un cammino valido
    if (minDistance[targetVertex] == numeric_limits<double>::infinity()) {
        cout << "Nessun cammino trovato tra i vertici " << startVertex << " e " << targetVertex << "." << endl;
        mesh.Cell0DsMarker.assign(totalVertices, 0);
        mesh.Cell1DsMarker.assign(totalEdges, 0);
        return pathResult;
    }

    // Ricostruisci il cammino minimo a ritroso
    mesh.Cell0DsMarker.assign(totalVertices, 0); // Reset marker vertici
    mesh.Cell1DsMarker.assign(totalEdges, 0);    // Reset marker lati

    pathResult.totalLength = minDistance[targetVertex];
    unsigned int currentIndex = targetVertex;

    while (currentIndex != startVertex) {
        mesh.Cell0DsMarker[currentIndex] = 1; // Marca il vertice come parte del cammino

        unsigned int prevIndex = previousVertex[currentIndex];
        unsigned int usedEdgeId = previousEdge[currentIndex];

        mesh.Cell1DsMarker[usedEdgeId] = 1;   // Marca il lato usato
        pathResult.numEdges++;

        currentIndex = prevIndex;
    }

    // Marca anche il vertice di partenza
    mesh.Cell0DsMarker[startVertex] = 1;

    return pathResult;
}

}