#include <iostream>
#include <cstdlib> 
#include "PolyhedralMesh.hpp"  
#include "Utils.hpp"           

using namespace std;
using namespace Eigen;
using namespace PolyhedralLibrary;

int main(int argc, char *argv[]) {
    
    // === [1] Inizializzazione variabili ===
    // ID dei vertici per il calcolo del cammino minimo
    unsigned int startVertexId = 0;
    unsigned int endVertexId = 0;
    // Flag che attiva il calcolo del cammino minimo
    bool calculatePath = false;

    // === [2] Parsing argomenti da linea di comando ===
    // L'utente può inserire 4 o 6 parametri:
    // - 4 per la generazione della mesh: p q b c
    // - 6 per includere anche il calcolo del cammino minimo: p q b c start end
    if (argc == 5) {
        cout << "Modalità: Generazione mesh.\n";
    } 
    else if (argc == 7) {
        cout << "Modalità: Generazione mesh e calcolo cammino minimo.\n";
        startVertexId = stoul(argv[5]);
        endVertexId = stoul(argv[6]);
        calculatePath = true;
    } 
    else {
        // Formato non valido
        cerr << "Uso corretto:\n";
        cerr << "  " << argv[0] << " p q b c\n";
        cerr << "  " << argv[0] << " p q b c startVertexId endVertexId\n";
        return 1;
    }

    // === [3] Conversione degli argomenti in interi ===
    int p = stoi(argv[1]);  // numero di facce che si incontrano in un vertice
    int q = stoi(argv[2]);  // numero di lati per faccia
    int b = stoi(argv[3]);  // numero di strati base
    int c = stoi(argv[4]);  // numero di strati superiori

    cout << "Parametri inseriti: p=" << p << ", q=" << q << ", b=" << b << ", c=" << c << "\n";
    if (calculatePath)
        cout << "Cammino minimo da vertice " << startVertexId << " a vertice " << endVertexId << "\n";

    // === [4] Controllo validità dei parametri p e q ===
    // Solo i valori regolari sono ammessi (3, 4 o 5)
    if ((p < 3 || p > 5) || (q < 3 || q > 5)) {
        cerr << "Errore: p e q devono essere 3, 4 o 5.\n";
        return 1;
    }

    // === [5] Dichiarazione delle strutture mesh ===
    PolyhedralMesh mesh;       // Mesh di partenza (grezza)
    PolyhedralMesh meshFinal;  // Mesh finale (dopo triangolazione)

    // === [6] Costruzione e triangolazione della mesh ===
    if (b != c) {
        // Caso: mesh "asimmetrica" (strati base ≠ strati superiori)

        if (p == 3 && q == 3) {
            generaTetraedro(mesh);
            Triangolazione(q, b, c, mesh, meshFinal);

        } else if (p == 3) {
            // Se p = 3, il poliedro base è un solido platonico con 3 lati per faccia
            (q == 4) ? generaOttaedro(mesh) : generaIcosaedro(mesh);
            Triangolazione(q, b, c, mesh, meshFinal);

        } else if (q == 3) {
            // Invertiamo p e q se q = 3, perché servono per il duale
            (p == 4) ? generaOttaedro(mesh) : generaIcosaedro(mesh);
            invertiValori(p, q);
            TriangolazioneDuale(q, b, c, mesh, meshFinal);
        }

    } else {
        // Caso: mesh "simmetrica" (stessi strati sopra e sotto)

        if (p == 3 && q == 3) {
            generaTetraedro(mesh);
            Triangolazione2(q, b, mesh, meshFinal);

        } else if (p == 3) {
            (q == 4) ? generaOttaedro(mesh) : generaIcosaedro(mesh);
            Triangolazione2(q, b, mesh, meshFinal);

        } else if (q == 3) {
            (p == 4) ? generaOttaedro(mesh) : generaIcosaedro(mesh);
            invertiValori(p, q);
            Triangolazione2(q, b, mesh, meshFinal);
        }
    }

    // === [7] Calcolo del cammino minimo (se richiesto) ===
    if (calculatePath) {
        // Costruisce la matrice di adiacenza della mesh finale
        MatrixXi adjMatrix = calculateAdjacencyMatrix(meshFinal);

        // Esegue l'algoritmo di Dijkstra per trovare il cammino più breve
        ShortestPathResult pathResult = findShortestPathDijkstra(
            meshFinal,
            adjMatrix,
            startVertexId,
            endVertexId
        );

        if (pathResult.numEdges > 0 || startVertexId == endVertexId) {
            // Cammino trovato o i vertici coincidono
            cout << "\n--- Cammino Minimo ---\n";
            cout << "Numero di lati nel cammino: " << pathResult.numEdges << endl;
            cout << "Lunghezza totale del cammino: " << pathResult.totalLength << endl;
        } else {
            // Nessun percorso valido
            cout << "\nNessun cammino trovato tra " << startVertexId 
                 << " e " << endVertexId << ".\n";
        }
    }

    // === [8] Proiezione su sfera unitaria e visualizzazione ===
    ProjectMeshToUnitSphere(meshFinal);  // Normalizza i vertici sulla sfera
    ExportParaview(meshFinal);           // Esporta per visualizzazione con Paraview

    // === [9] Esportazione su file .txt per analisi (celle di ogni dimensione) ===
    WriteCell0Ds(meshFinal);  // Vertici
    WriteCell1Ds(meshFinal);  // Spigoli
    WriteCell2Ds(meshFinal);  // Facce
    WriteCell3Ds(meshFinal);  // Celle (poliedri)

    return 0;
	
}