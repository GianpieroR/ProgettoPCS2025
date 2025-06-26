#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <fstream>
#include <sstream>

using namespace std;
namespace PolyhedralLibrary {

    // Funzione per scrivere le coordinate dei vertici (Cell0Ds) del mesh su un file di testo
void WriteCell0Ds(const PolyhedralMesh& mesh) {
    // Apri il file "Cell0Ds.txt" in modalità scrittura (sovrascrive il file esistente)
    ofstream file("Cell0Ds.txt");

    // Controlla se il file è stato aperto correttamente
    if (!file.is_open()) {
        // Se non è stato possibile aprire il file, stampa un messaggio di errore e termina la funzione
        cerr << "Errore: impossibile aprire il file Cell0Ds.txt per la scrittura." << endl;
        return;
    }

    // Scrivi l'intestazione delle colonne nel file, separata da punto e virgola
    // Indica che le colonne saranno: Id, x, y, z
    file << "Id;x;y;z\n";

    // Itera su tutti gli indici dei vertici presenti nel mesh
    for (size_t i = 0; i < mesh.Cell0DsId.size(); ++i) {
        // Scrivi per ogni vertice:
        // - il suo Id
        // - la coordinata x (prima riga della matrice Cell0DsCoordinates, colonna i)
        // - la coordinata y (seconda riga, colonna i)
        // - la coordinata z (terza riga, colonna i)
        // Ogni valore è separato da un punto e virgola e ogni vertice su una nuova riga
        file << mesh.Cell0DsId[i] << ";"
             << mesh.Cell0DsCoordinates(0, i) << ";"
             << mesh.Cell0DsCoordinates(1, i) << ";"
             << mesh.Cell0DsCoordinates(2, i) << "\n";
    }

    // Chiudi il file dopo aver scritto tutti i dati
    file.close();
}

    // Funzione per scrivere gli archi (Cell1Ds) del mesh su un file di testo
void WriteCell1Ds(const PolyhedralMesh& mesh) {
    // Apro il file "Cell1Ds.txt" in modalità scrittura
    ofstream file("Cell1Ds.txt");

    // Controllo se il file è stato aperto correttamente
    if (!file.is_open()) {
        cerr << "Errore: impossibile aprire il file Cell1Ds.txt per la scrittura." << endl;
        return; // Esco dalla funzione in caso di errore
    }

    // Scrivo l'intestazione delle colonne: Id;Start;End
    file << "Id;Start;End\n";

    // Ciclo su tutti gli archi del mesh
    for (size_t i = 0; i < mesh.Cell1DsId.size(); ++i) {
        // Per ogni arco scrivo:
        // - il suo Id
        // - l'indice del vertice di inizio (estremo 0)
        // - l'indice del vertice di fine (estremo 1)
        file << mesh.Cell1DsId[i] << ";"
             << mesh.Cell1DsExtrema(i, 0) << ";"
             << mesh.Cell1DsExtrema(i, 1) << "\n";
    }

    // Chiudo il file dopo aver scritto tutti i dati
    file.close();
}


// Funzione per scrivere le facce (Cell2Ds) del mesh su un file di testo
void WriteCell2Ds(const PolyhedralMesh& mesh){
    // Apro il file "Cell2Ds.txt" in modalità scrittura
    std::ofstream file("Cell2Ds.txt");

    // Controllo se il file è stato aperto correttamente
    if (!file.is_open()) {
        cerr << "Errore: impossibile aprire il file Cell2Ds.txt per la scrittura." << endl;
        return; // Esco dalla funzione in caso di errore
    }

    // Scrivo l'intestazione delle colonne:
    // Id;NumVertices;Vertices;NumEdges;Edges
    // Dove NumVertices e NumEdges indicano il numero degli elementi nelle liste
    file << "Id;NumVertices;Vertices;NumEdges;Edges\n";

    // Ciclo su tutte le facce del mesh
    for (size_t i = 0; i < mesh.Cell2DsId.size(); ++i) {
        // Scrivo l'Id della faccia
        file << mesh.Cell2DsId[i] << ";";

        // Scrivo il numero di vertici (qui sempre 3 per triangolo)
        file << "3;";

        // Scrivo la lista degli indici dei vertici separati da ';'
        for (size_t j = 0; j < mesh.Cell2DsVertices[i].size(); ++j) {
            file << mesh.Cell2DsVertices[i][j] << ";";
        }

        // Scrivo il numero di archi (qui sempre 3 per triangolo)
        file << "3;";

        // Scrivo la lista degli indici degli archi separati da ';'
        for (size_t j = 0; j < mesh.Cell2DsEdges[i].size(); ++j) {
            file << mesh.Cell2DsEdges[i][j] << ";";
        }

        // Vado a capo per la prossima faccia
        file << "\n";
    }

    // Chiudo il file dopo aver scritto tutti i dati
    file.close();
}
    // stesso procedimento delle precedenti WriteCell

    void WriteCell3Ds(const PolyhedralMesh& mesh) {
        ofstream file("Cell3Ds.txt");
		if (!file.is_open()) {
        cerr << "Errore: impossibile aprire il file Cell3Ds.txt per la scrittura." << endl;
        return;
    	}
        file << "Id;NumVertices;Vertices;NumEdges;Edges;NumFaces;Faces\n";
        file << mesh.Cell3DsId << ";";
        
        file << mesh.NumCells0Ds << ";";
        for (size_t i = 0; i < mesh.Cell3DsVertices.size(); ++i) {
            file << mesh.Cell3DsVertices[i] << ";";
        }
        file << mesh.NumCells1Ds << ";";
        for (size_t i = 0; i < mesh.Cell3DsEdges.size(); ++i) {
            file << mesh.Cell3DsEdges[i] << ";";
        }
        file << mesh.NumCells2Ds << ";";
        for (size_t i = 0; i < mesh.Cell3DsFaces.size(); ++i) {
            file << mesh.Cell3DsFaces[i] << ";";
        }

        file.close();
    }

}