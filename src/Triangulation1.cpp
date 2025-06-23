#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <limits>

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary {

// Trova o crea un edge tra i vertici a e b per il triangolo corrente
void FindAddEdge(
    const int vertA,
    const int vertB,
    PolyhedralMesh& meshOut,
    unsigned int& nextEdgeIndex,
    const unsigned int faceIndex)
{
    bool edgeFound = false;

    for (unsigned int e = 0; e <= nextEdgeIndex; ++e) {
        // Se l'edge esiste (in ordine qualsiasi)
        const auto& endpoints = meshOut.Cell1DsExtrema.row(e);
        if ((endpoints(0) == vertA && endpoints(1) == vertB) ||
            (endpoints(0) == vertB && endpoints(1) == vertA)) {

            meshOut.Cell2DsEdges[faceIndex].push_back(e);
            edgeFound = true;
            break;
        }
    }

    if (!edgeFound) {
        // Nuovo edge: aggiunta e configurazione flag
        meshOut.Cell1DsExtrema.row(nextEdgeIndex) << vertA, vertB;
        meshOut.Cell1DsId[nextEdgeIndex] = nextEdgeIndex;

        const auto& flagsA = meshOut.Cell0DsFlag[vertA];
        const auto& flagsB = meshOut.Cell0DsFlag[vertB];
        bool hasCommonFlag = false;

        // Determina se i vertici condividono uno stesso flag
        for (auto fA : flagsA) {
            for (auto fB : flagsB) {
                if (fA == fB) {
                    meshOut.Cell1DsFlag[nextEdgeIndex] = fA;
                    hasCommonFlag = true;
                    break;
                }
            }
            if (hasCommonFlag) break;
        }

        if (!hasCommonFlag) {
            meshOut.Cell1DsFlag[nextEdgeIndex] = numeric_limits<unsigned int>::max();
        }

        meshOut.Cell2DsEdges[faceIndex].push_back(nextEdgeIndex);
        ++nextEdgeIndex;
    }
}

// Triangolazione della mesh: suddivisione di ogni faccia in piccoli triangoli
void triangulateAndStore(
    PolyhedralMesh& meshIn,
    PolyhedralMesh& meshOut,
    const unsigned int b,
    const unsigned int c,
    const vector<int>& dimDup)
{
    unsigned int subdivision = b + c;

    // --- Inizializzazione storage ---
    meshOut.Cell0DsId.resize(dimDup[0]);
    meshOut.Cell0DsCoordinates = MatrixXd::Zero(3, dimDup[0]);
    meshOut.Cell0DsFlag.resize(dimDup[0]);

    meshOut.Cell1DsId.resize(dimDup[1]);
    meshOut.Cell1DsExtrema = MatrixXi::Zero(dimDup[1], 2);
    meshOut.Cell1DsFlag.resize(dimDup[1]);

    meshOut.Cell2DsId.resize(dimDup[2]);
    meshOut.Cell2DsVertices.resize(dimDup[2]);
    meshOut.Cell2DsEdges.resize(dimDup[2]);

    unsigned int vertexCounter = 0;
    unsigned int edgeCounter = 0;
    unsigned int faceCounter = 0;

    // Loop su ogni faccia della mesh originale
    for (unsigned int fId = 0; fId < meshIn.Cell2DsId.size(); ++fId) {
        const auto& faceVerts = meshIn.Cell2DsVertices[fId];
        Vector3d vA = meshIn.Cell0DsCoordinates.col(faceVerts[0]);
        Vector3d vB = meshIn.Cell0DsCoordinates.col(faceVerts[1]);
        Vector3d vC = meshIn.Cell0DsCoordinates.col(faceVerts[2]);

        // Genera griglia di vertici suddivisi
        vector<vector<int>> vertexGrid;
        vertexGrid.reserve(subdivision + 1);

        for (unsigned int i = 0; i <= subdivision; ++i) {
            vector<int> rowIndices;
            rowIndices.reserve(i + 1);

            Vector3d lineStart = (i/subdivision) * vB + ((subdivision-i)/static_cast<double>(subdivision)) * vA;
            Vector3d lineEnd   = (i/subdivision) * vC + ((subdivision-i)/static_cast<double>(subdivision)) * vA;

            for (unsigned int j = 0; j <= i; ++j) {
                Vector3d newPoint = (i == 0)
                    ? vA
                    : (j/static_cast<double>(i)) * lineEnd + ((i-j)/static_cast<double>(i)) * lineStart;

                meshOut.Cell0DsCoordinates.col(vertexCounter) = newPoint;
                meshOut.Cell0DsId[vertexCounter] = vertexCounter;

                // Imposta i flag di posizione sui bordi
                if (i == 0) {
                    meshOut.Cell0DsFlag[vertexCounter] = { meshIn.Cell2DsEdges[fId][0],
                                                           meshIn.Cell2DsEdges[fId][2] };
                } else if (i == subdivision) {
                    if (j == 0)
                        meshOut.Cell0DsFlag[vertexCounter] = { meshIn.Cell2DsEdges[fId][0],
                                                               meshIn.Cell2DsEdges[fId][1] };
                    else if (j == subdivision)
                        meshOut.Cell0DsFlag[vertexCounter] = { meshIn.Cell2DsEdges[fId][1],
                                                               meshIn.Cell2DsEdges[fId][2] };
                    else
                        meshOut.Cell0DsFlag[vertexCounter] = { meshIn.Cell2DsEdges[fId][1] };
                } else if (j == 0) {
                    meshOut.Cell0DsFlag[vertexCounter] = { meshIn.Cell2DsEdges[fId][0] };
                } else if (j == i) {
                    meshOut.Cell0DsFlag[vertexCounter] = { meshIn.Cell2DsEdges[fId][2] };
                } else {
                    meshOut.Cell0DsFlag[vertexCounter] = { numeric_limits<unsigned int>::max() };
                }

                rowIndices.push_back(vertexCounter++);
            }
            vertexGrid.emplace_back(move(rowIndices));
        }

        // Creazione triangoli nella griglia
        for (unsigned int i = 0; i < subdivision; ++i) {
            for (unsigned int j = 0; j < i; ++j) {
                // Triangolo basso
                auto triA = vector<unsigned int>{
                    vertexGrid[i][j],
                    vertexGrid[i+1][j],
                    vertexGrid[i+1][j+1]
                };
                meshOut.Cell2DsVertices[faceCounter] = triA;
                meshOut.Cell2DsId[faceCounter] = faceCounter;
                for (unsigned int k = 0; k < 3; ++k)
                    FindAddEdge(triA[k], triA[(k+1)%3], meshOut, edgeCounter, faceCounter);
                faceCounter++;

                // Triangolo alto
                auto triB = vector<unsigned int>{
                    vertexGrid[i][j],
                    vertexGrid[i+1][j+1],
                    vertexGrid[i][j+1]
                };
                meshOut.Cell2DsVertices[faceCounter] = triB;
                meshOut.Cell2DsId[faceCounter] = faceCounter;
                for (unsigned int k = 0; k < 3; ++k)
                    FindAddEdge(triB[k], triB[(k+1)%3], meshOut, edgeCounter, faceCounter);
                faceCounter++;
            }

            // Triangolo di fine riga
            auto triC = vector<unsigned int>{
                vertexGrid[i][i],
                vertexGrid[i+1][i],
                vertexGrid[i+1][i+1]
            };
            meshOut.Cell2DsVertices[faceCounter] = triC;
            meshOut.Cell2DsId[faceCounter] = faceCounter;
            for (unsigned int k = 0; k < 3; ++k)
                FindAddEdge(triC[k], triC[(k+1)%3], meshOut, edgeCounter, faceCounter);
            faceCounter++;
        }
    }
}