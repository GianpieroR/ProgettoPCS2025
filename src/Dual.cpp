#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <fstream>
#include <sstream>
#include <Eigen/Dense> 
#include <cmath>

using namespace std;
using namespace Eigen;
namespace PolyhedralLibrary{

Vector3d getFaceBarycenter(const PolyhedralMesh& meshTriangulated, const unsigned int faceId) {
    Vector3d barycenter = Vector3d::Zero();
    const auto& verticesOfFace = meshTriangulated.Cell2DsVertices[faceId]; // Vertici della faccia
    for (unsigned int vertexId : verticesOfFace) {
        barycenter += meshTriangulated.Cell0DsCoordinates.col(vertexId);
    }
    return barycenter /= 3.0;
}

map<pair<unsigned int, unsigned int>, vector<unsigned int>> buildEdgeToFacesMap(const PolyhedralMesh& meshTriangulated) {
    map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap;

    for (unsigned int faceId = 0; faceId < meshTriangulated.Cell2DsId.size(); ++faceId) {
        const vector<unsigned int>& edgesOfFace = meshTriangulated.Cell2DsEdges[faceId];
        for (unsigned int edgeId : edgesOfFace) {
            unsigned int vertex1 = meshTriangulated.Cell1DsExtrema(edgeId, 0);
            unsigned int vertex2 = meshTriangulated.Cell1DsExtrema(edgeId, 1);
            pair<unsigned int, unsigned int> sortedEdge = {min(vertex1, vertex2), max(vertex1, vertex2)};
            edgeToFacesMap[sortedEdge].push_back(faceId);
        }
    }
    return edgeToFacesMap;
}

map<unsigned int, vector<unsigned int>> buildVertexToFacesMap(const PolyhedralMesh& meshTriangulated) {
    map<unsigned int, vector<unsigned int>> vertexToFacesMap;
    for (unsigned int faceId = 0; faceId < meshTriangulated.Cell2DsId.size(); ++faceId) {
        for (unsigned int vertexId : meshTriangulated.Cell2DsVertices[faceId]) {
            vertexToFacesMap[vertexId].push_back(faceId);
        }
    }
    return vertexToFacesMap;
}

map<unsigned int, vector<unsigned int>> buildVertexToEdgesMap(const PolyhedralMesh& meshTriangulated) {
    map<unsigned int, vector<unsigned int>> vertexToEdgesMap;
    for (unsigned int edgeId = 0; edgeId < meshTriangulated.Cell1DsExtrema.rows(); ++edgeId) {
        unsigned int vertex1 = meshTriangulated.Cell1DsExtrema(edgeId, 0);
        unsigned int vertex2 = meshTriangulated.Cell1DsExtrema(edgeId, 1);
        vertexToEdgesMap[vertex1].push_back(edgeId);
        vertexToEdgesMap[vertex2].push_back(edgeId);
    }
    return vertexToEdgesMap;
}

void CalculateDual(PolyhedralMesh& meshTriangulated, PolyhedralMesh& meshDual, map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap) {
    // VERTICI DEL DUALE
    meshDual.Cell0DsId.resize(meshTriangulated.Cell2DsId.size());
    meshDual.Cell0DsCoordinates = MatrixXd::Zero(3, meshTriangulated.Cell2DsId.size());
    
    for (unsigned int faceId = 0; faceId < meshTriangulated.Cell2DsId.size(); ++faceId){
        meshDual.Cell0DsCoordinates.col(faceId) = getFaceBarycenter(meshTriangulated, faceId);
        meshDual.Cell0DsId[faceId] = faceId;
    }

    // SPIGOLI DEL DUALE
    vector<pair<unsigned int, unsigned int>> dualEdges;

    for (const auto& edgeFacesEntry : edgeToFacesMap) {
        const vector<unsigned int>& facesSharingEdge = edgeFacesEntry.second;
        if (facesSharingEdge.size() == 2) {
            unsigned int face1 = facesSharingEdge[0];
            unsigned int face2 = facesSharingEdge[1];
            pair<unsigned int, unsigned int> dualEdge = {min(face1, face2), max(face1, face2)};
            dualEdges.push_back(dualEdge);
        }
    }

    meshDual.Cell1DsId.resize(dualEdges.size());
    meshDual.Cell1DsExtrema = MatrixXi::Zero(dualEdges.size(), 2);

    for (unsigned int i = 0; i < dualEdges.size(); ++i) {
        meshDual.Cell1DsId[i] = i;
        meshDual.Cell1DsExtrema(i, 0) = dualEdges[i].first;
        meshDual.Cell1DsExtrema(i, 1) = dualEdges[i].second;
    }

    // FACCE DEL DUALE
    map<pair<unsigned int, unsigned int>, unsigned int> dualEdgeToIdMap;
    for (unsigned int i = 0; i < meshDual.Cell1DsId.size(); ++i) {
        unsigned int v1 = meshDual.Cell1DsExtrema(i, 0);
        unsigned int v2 = meshDual.Cell1DsExtrema(i, 1);
        dualEdgeToIdMap[{min(v1, v2), max(v1, v2)}] = i;
    }

    map<unsigned int, vector<unsigned int>> vertexToFacesMap = buildVertexToFacesMap(meshTriangulated);
    map<unsigned int, vector<unsigned int>> vertexToEdgesMap = buildVertexToEdgesMap(meshTriangulated);

    meshDual.Cell2DsVertices.resize(vertexToFacesMap.size());
    meshDual.Cell2DsId.resize(vertexToFacesMap.size());
    meshDual.Cell2DsEdges.resize(vertexToFacesMap.size());

    unsigned int dualFaceId = 0;

    for (const auto& vertexFacesEntry : vertexToFacesMap) {
        unsigned int originalVertexId = vertexFacesEntry.first;
        vector<unsigned int> incidentFaces = vertexFacesEntry.second;
        vector<unsigned int> incidentEdges = vertexToEdgesMap[originalVertexId];

        if (incidentFaces.size() < 3 || incidentEdges.size() < 3) {
            cerr << "Warning: Original vertex " << originalVertexId << " has less than 3 incident faces or edges. Dual face cannot be formed." << endl;
            continue;
        }

        vector<unsigned int> orderedDualVertices;
        vector<unsigned int> dualFaceEdges;

        unsigned int currentOriginalEdgeId = incidentEdges[0];

        pair<unsigned int, unsigned int> currentEdgeKey = {
            min(meshTriangulated.Cell1DsExtrema(currentOriginalEdgeId, 0), meshTriangulated.Cell1DsExtrema(currentOriginalEdgeId, 1)),
            max(meshTriangulated.Cell1DsExtrema(currentOriginalEdgeId, 0), meshTriangulated.Cell1DsExtrema(currentOriginalEdgeId, 1))
        };

        vector<unsigned int> facesOnCurrentEdge = edgeToFacesMap.at(currentEdgeKey);

        if (facesOnCurrentEdge.size() != 2) {
            cerr << "Error: Original edge " << currentOriginalEdgeId << " not shared by exactly two faces." << endl;
            continue;
        }

        unsigned int startFaceId = facesOnCurrentEdge[0];
        unsigned int previousDualVertexId = startFaceId;

        orderedDualVertices.push_back(previousDualVertexId);

        for (size_t i = 0; i < incidentEdges.size(); ++i) {
            vector<unsigned int> currentEdgeFaces = edgeToFacesMap.at(currentEdgeKey);

            unsigned int faceA = currentEdgeFaces[0];
            unsigned int faceB = currentEdgeFaces[1];

            unsigned int nextDualVertexId;
            if (faceA == previousDualVertexId) nextDualVertexId = faceB;
            else if (faceB == previousDualVertexId) nextDualVertexId = faceA;
            else {
                cerr << "Logical error: current face not found among faces sharing edge." << endl;
                break;
            }

            pair<unsigned int, unsigned int> dualEdgeKey = {min(previousDualVertexId, nextDualVertexId), max(previousDualVertexId, nextDualVertexId)};
            auto dualEdgeIt = dualEdgeToIdMap.find(dualEdgeKey);

            if (dualEdgeIt != dualEdgeToIdMap.end()) {
                dualFaceEdges.push_back(dualEdgeIt->second);
            } else {
                cerr << "Error: Dual edge not found for faces " << previousDualVertexId << " and " << nextDualVertexId
                     << " around original vertex " << originalVertexId << endl;
                break;
            }

            if (nextDualVertexId == startFaceId) {
                if (orderedDualVertices.size() == incidentFaces.size()) break;
                else {
                    cerr << "Topology error: incomplete cycle for vertex " << originalVertexId << endl;
                    break;
                }
            }

            bool alreadyAdded = false;
            for (unsigned int v : orderedDualVertices) {
                if (v == nextDualVertexId) {
                    alreadyAdded = true;
                    break;
                }
            }
            if (!alreadyAdded) orderedDualVertices.push_back(nextDualVertexId);

            previousDualVertexId = nextDualVertexId;

            bool foundNextEdge = false;
            const auto& edgesOfNextFace = meshTriangulated.Cell2DsEdges[nextDualVertexId];
            for (unsigned int edgeNextFace : edgesOfNextFace) {
                if (edgeNextFace == currentOriginalEdgeId) continue;

                unsigned int eV1 = meshTriangulated.Cell1DsExtrema(edgeNextFace, 0);
                unsigned int eV2 = meshTriangulated.Cell1DsExtrema(edgeNextFace, 1);

                if (eV1 == originalVertexId || eV2 == originalVertexId) {
                    currentOriginalEdgeId = edgeNextFace;
                    currentEdgeKey = {min(eV1, eV2), max(eV1, eV2)};
                    foundNextEdge = true;
                    break;
                }
            }

            if (!foundNextEdge) {
                cerr << "Error: Could not find next edge in cycle around vertex " << originalVertexId << endl;
                break;
            }
        }

        if (!orderedDualVertices.empty()) {
            meshDual.Cell2DsId[dualFaceId] = dualFaceId;
            meshDual.Cell2DsVertices[dualFaceId] = orderedDualVertices;
            meshDual.Cell2DsEdges[dualFaceId] = dualFaceEdges;
            dualFaceId++;
        }
    }
}

void ProjectMeshToUnitSphere(PolyhedralMesh& meshTriangulated) {
    for (int vertexIdx = 0; vertexIdx < meshTriangulated.Cell0DsCoordinates.cols(); ++vertexIdx) {
        Vector3d vertexCoords = meshTriangulated.Cell0DsCoordinates.col(vertexIdx);
        double length = vertexCoords.norm();
        if (length < 1e-12) {
            cerr << "Warning: Vertex " << vertexIdx << " too close to origin; projection skipped." << endl;
            continue;
        }
        meshTriangulated.Cell0DsCoordinates.col(vertexIdx) = vertexCoords / length;
    }
}

}