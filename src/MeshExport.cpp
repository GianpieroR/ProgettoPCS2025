#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <fstream>
#include <sstream>

using namespace std;
namespace PolyhedralLibrary {

    void WriteCell0Ds(const PolyhedralMesh& mesh) {
        ofstream file("Cell0Ds.txt");
		if (!file.is_open()) {
        cerr << "Errore: impossibile aprire il file Cell0Ds.txt per la scrittura." << endl;
        return;
    	}
        file << "Id;x;y;z\n";
        for (size_t i = 0; i < mesh.Cell0DsId.size(); ++i) {
            file << mesh.Cell0DsId[i] << ";"
                 << mesh.Cell0DsCoordinates(0, i) << ";"
                 << mesh.Cell0DsCoordinates(1, i) << ";"
                 << mesh.Cell0DsCoordinates(2, i) << "\n";
    	}
    file.close();
    }

    void WriteCell1Ds(const PolyhedralMesh& mesh) {
        ofstream file("Cell1Ds.txt");
		if (!file.is_open()) {
        cerr << "Errore: impossibile aprire il file Cell1Ds.txt per la scrittura." << endl;
        return;
    	}
        file << "Id;Start;End";
        for (size_t i = 0; i < mesh.Cell1DsId.size(); ++i) {
            file << mesh.Cell1DsId[i] << ";"
                 << mesh.Cell1DsExtrema(i, 0) << ";"
                 << mesh.Cell1DsExtrema(i, 1) << "\n";
        }
        file.close();
    }

    void WriteCell2Ds(const PolyhedralMesh& mesh){
        std::ofstream file("Cell2Ds.txt");
		if (!file.is_open()) {
        cerr << "Errore: impossibile aprire il file Cell2Ds.txt per la scrittura." << endl;
        return;
    }
        file << "Id;NumVertices;Vertices;NumEdges;Edges\n";
        for (size_t i = 0; i < mesh.Cell2DsId.size(); ++i) {
            file << mesh.Cell2DsId[i] << ";" << "3;";
            for (size_t j = 0; j < mesh.Cell2DsVertices[i].size(); ++j) {
                file << mesh.Cell2DsVertices[i][j] << ";";
            }
            file << "3;";
            for (size_t j = 0; j < mesh.Cell2DsEdges[i].size(); ++j) {
                file << mesh.Cell2DsEdges[i][j] << ";";
            }
            file << "\n";
        }
        file.close();
    }

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