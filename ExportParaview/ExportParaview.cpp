#include "UCDUtilities.hpp"
#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <Eigen/Dense>

using namespace std;
namespace PolyhedralLibrary
{
	// Funzione per esportare la mesh in formato Paraview tramite file .inp
	void ExportParaview(const PolyhedralMesh& meshData){

        Gedim::UCDUtilities ucdUtil;

        // Vettore per memorizzare gli ID dei vertici in formato Eigen::VectorXi
        Eigen::VectorXi vertexIndices(meshData.Cell0DsId.size());
        for (size_t idx = 0; idx < meshData.Cell0DsId.size(); ++idx)
            vertexIndices[idx] = static_cast<int>(meshData.Cell0DsId[idx]);

        // Vettore per memorizzare gli ID dei lati (edges) in formato Eigen::VectorXi
        Eigen::VectorXi edgeIndices(meshData.Cell1DsId.size());
        for (size_t idx = 0; idx < meshData.Cell1DsId.size(); ++idx)
            edgeIndices[idx] = static_cast<int>(meshData.Cell1DsId[idx]);

        // Controlla se i marker sui vertici e sugli spigoli sono presenti e validi
        bool validVertexMarkers = !meshData.Cell0DsMarker.empty();
        bool validEdgeMarkers = !meshData.Cell1DsMarker.empty();

        if (validVertexMarkers && validEdgeMarkers) {
            // Creo vettore per le proprietà scalari sui vertici (es. marcatori di percorso)
            vector<double> scalarDataVertices;
            scalarDataVertices.reserve(meshData.Cell0DsMarker.size());
            for (auto mark : meshData.Cell0DsMarker)
                scalarDataVertices.push_back(static_cast<double>(mark));

            // Definisco la proprietà scalare associata ai vertici
            Gedim::UCDProperty<double> scalarPropertyVertices;
            scalarPropertyVertices.Label = "PathVertex";   // Nome proprietà
            scalarPropertyVertices.UnitLabel = "-";        // Unità (nessuna)
            scalarPropertyVertices.NumComponents = 1;      // Componente singola
            scalarPropertyVertices.Data = scalarDataVertices.data();
            scalarPropertyVertices.Size = scalarDataVertices.size();

            // Creo vettore per le proprietà scalari sugli spigoli
            vector<double> scalarDataEdges;
            scalarDataEdges.reserve(meshData.Cell1DsMarker.size());
            for (auto mark : meshData.Cell1DsMarker)
                scalarDataEdges.push_back(static_cast<double>(mark));

            // Definisco la proprietà scalare associata agli spigoli
            Gedim::UCDProperty<double> scalarPropertyEdges;
            scalarPropertyEdges.Label = "PathEdge";
            scalarPropertyEdges.UnitLabel = "-";
            scalarPropertyEdges.NumComponents = 1;
            scalarPropertyEdges.Data = scalarDataEdges.data();
            scalarPropertyEdges.Size = scalarDataEdges.size();

            // Esporto i punti (vertici) con le proprietà scalari
            ucdUtil.ExportPoints("./Cell0Ds.inp",
                                 meshData.Cell0DsCoordinates,
                                 { scalarPropertyVertices },
                                 vertexIndices);

            // Esporto i segmenti (lati) con le proprietà scalari su vertici e spigoli
            ucdUtil.ExportSegments("./Cell1Ds.inp",
                                   meshData.Cell0DsCoordinates,
                                   meshData.Cell1DsExtrema.transpose(),
                                   { scalarPropertyVertices },
                                   { scalarPropertyEdges },
                                   edgeIndices);

        } else {
            // Esporto senza proprietà scalari, solo la geometria
            ucdUtil.ExportPoints("./Cell0Ds.inp",
                                 meshData.Cell0DsCoordinates,
                                 {},  // nessuna proprietà sui vertici
                                 Eigen::VectorXi());

            ucdUtil.ExportSegments("./Cell1Ds.inp",
                                   meshData.Cell0DsCoordinates,
                                   meshData.Cell1DsExtrema.transpose(),
                                   {}, // nessuna proprietà sui vertici
                                   {}, // nessuna proprietà sugli spigoli
                                   Eigen::VectorXi());
        }
	}

	// Funzione di debug per stampare a console i dettagli della mesh triangolata
	void printMeshTriangulated(const PolyhedralMesh& meshData) {

		// Stampa gli ID dei vertici
		cout << "Vertici (Cell0DsId): ";
		for (auto vertexId : meshData.Cell0DsId) cout << vertexId << " ";
		cout << "\nCoordinate dei vertici (per colonne):" << endl;
		for (int col = 0; col < meshData.Cell0DsCoordinates.cols(); ++col) {
			cout << "Colonna " << col << ": ";
			for (int row = 0; row < meshData.Cell0DsCoordinates.rows(); ++row) {
				cout << meshData.Cell0DsCoordinates(row, col) << " ";
			}
			cout << endl;
		}

		// Stampa i flag associati ai vertici
		cout << "Flag vertici (Cell0DsFlag):" << endl;
		for (const auto& flagRow : meshData.Cell0DsFlag) {
			for (auto flagValue : flagRow) cout << flagValue << " ";
			cout << endl;
		}

		// Stampa gli ID dei lati (edges)
		cout << "Lati (Cell1DsId): ";
		for (auto edgeId : meshData.Cell1DsId) cout << edgeId << " ";
		cout << "\nEstremi dei lati (Cell1DsExtrema, per righe):" << endl;
		for (int row = 0; row < meshData.Cell1DsExtrema.rows(); ++row) {
			for (int col = 0; col < meshData.Cell1DsExtrema.cols(); ++col) {
				cout << meshData.Cell1DsExtrema(row, col) << " ";
			}
			cout << endl;
		}

		// Stampa i flag associati ai lati
		cout << "Flag lati (Cell1DsFlag):" << endl;
		for (const auto& flagRow : meshData.Cell1DsFlag) {
			cout << flagRow << " ";
			cout << endl;
		}

		// Stampa i flag originali dei lati
		cout << "Flag originali lati (Cell1DsOriginalFlag):" << endl;
		for (const auto& flagRow : meshData.Cell1DsOriginalFlag) {
			cout << flagRow << " ";
			cout << endl;
		}

		// Stampa gli ID delle facce (Cell2Ds)
		cout << "Facce (Cell2DsId): ";
		for (auto faceId : meshData.Cell2DsId) cout << faceId << " ";

		cout << "\nVertici delle facce (Cell2DsVertices):" << endl;
		for (const auto& faceVertices : meshData.Cell2DsVertices) {
			for (auto vertex : faceVertices) cout << vertex << " ";
			cout << endl;
		}

		cout << "Lati delle facce (Cell2DsEdges):" << endl;
		for (const auto& faceEdges : meshData.Cell2DsEdges) {
			for (auto edge : faceEdges) cout << edge << " ";
			cout << endl;
		}

		// Informazioni sul poliedro (Cell3Ds)
		cout << "ID poliedro (Cell3DsId): " << meshData.Cell3DsId << endl;
		cout << "Numero di vertici: " << meshData.NumCells0Ds << endl;
		cout << "Numero di lati: " << meshData.NumCells1Ds << endl;
		cout << "Numero di facce: " << meshData.NumCells2Ds << endl;

		cout << "Vertici del poliedro (Cell3DsVertices): ";
		for (auto vertId : meshData.Cell3DsVertices) cout << vertId << " ";
		cout << endl;

		cout << "Lati del poliedro (Cell3DsEdges): ";
		for (auto edgeId : meshData.Cell3DsEdges) cout << edgeId << " ";
		cout << endl;

		cout << "Facce del poliedro (Cell3DsFaces): ";
		for (auto faceId : meshData.Cell3DsFaces) cout << faceId << " ";
		cout << endl;

		cout << "\n--- Fine struttura mesh ---" << endl;
	}
}