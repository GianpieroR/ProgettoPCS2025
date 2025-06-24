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
	/ Restituisce una versione ciclicamente normalizzata della sequenza di edge,
// cioè ruotata in modo da partire dal valore minimo.
vector<unsigned int> get_cyclic_normalized(const vector<unsigned int>& edges) {
	if (edges.size() <= 1)
		return edges;

	vector<unsigned int> rotated = edges;
	auto min_it = min_element(rotated.begin(), rotated.end());
	rotate(rotated.begin(), min_it, rotated.end());
	return rotated;
}

// Normalizza una faccia in modo ciclico e lessicografico (considera anche il verso opposto)
vector<unsigned int> NormalizeFaceEdges(const vector<unsigned int>& face_edges) {
	if (face_edges.size() <= 1)
		return face_edges;

	vector<unsigned int> norm_forward = get_cyclic_normalized(face_edges);
	vector<unsigned int> reversed = face_edges;
	reverse(reversed.begin(), reversed.end());
	vector<unsigned int> norm_reversed = get_cyclic_normalized(reversed);

	return (norm_forward < norm_reversed) ? norm_forward : norm_reversed;
}

// Aggiunge una nuova faccia solo se non è già presente (evita duplicati)
void FindAddFace(const vector<unsigned int>& vertex_ids,
                 const vector<unsigned int>& edge_ids,
                 PolyhedralMesh& meshOut,
                 unsigned int& faceCounter)
{
	vector<unsigned int> normalizedNew = NormalizeFaceEdges(edge_ids);
	bool alreadyExists = false;

	for (unsigned int i = 0; i < meshOut.Cell2DsId.size(); ++i) {
		vector<unsigned int> normalizedExisting = NormalizeFaceEdges(meshOut.Cell2DsEdges[i]);
		if (normalizedNew == normalizedExisting) {
			alreadyExists = true;
			break;
		}
	}

	if (!alreadyExists) {
		meshOut.Cell2DsVertices[faceCounter] = vertex_ids;
		meshOut.Cell2DsEdges[faceCounter] = edge_ids;
		meshOut.Cell2DsId[faceCounter] = faceCounter;
		faceCounter++;
	}
}

// Cerca un vertice con coordinate quasi uguali. Se esiste lo riusa, altrimenti lo crea.
unsigned int FindAddVertice(const Vector3d& coord, PolyhedralMesh& meshOut, unsigned int& vertexCounter) {
	double tol = 1e-12;
	for (unsigned int i = 0; i < meshOut.Cell0DsCoordinates.cols(); i++) {
		if ((meshOut.Cell0DsCoordinates.col(i) - coord).norm() < tol)
			return i;
	}
	meshOut.Cell0DsCoordinates.col(vertexCounter) = coord;
	meshOut.Cell0DsId[vertexCounter] = vertexCounter;
	return vertexCounter++;
}

// Cerca un edge esistente (in entrambi i versi) o ne crea uno nuovo
unsigned int FindAddEdge2(unsigned int a, unsigned int b, PolyhedralMesh& meshOut, unsigned int& edgeCounter) {
	for (unsigned int i = 0; i < meshOut.Cell1DsExtrema.rows(); ++i) {
		if ((meshOut.Cell1DsExtrema(i, 0) == static_cast<int>(a) && meshOut.Cell1DsExtrema(i, 1) == static_cast<int>(b)) ||
			(meshOut.Cell1DsExtrema(i, 0) == static_cast<int>(b) && meshOut.Cell1DsExtrema(i, 1) == static_cast<int>(a)))
			return i;
	}

	meshOut.Cell1DsExtrema(edgeCounter, 0) = a;
	meshOut.Cell1DsExtrema(edgeCounter, 1) = b;
	meshOut.Cell1DsId[edgeCounter] = edgeCounter;
	return edgeCounter++;
}

// Trova il baricentro della faccia adiacente a quella data tramite un edge
Vector3d FindNearBarycenter(const PolyhedralMesh& mesh, unsigned int edgeId, unsigned int currentFaceId,
							map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap) {

	unsigned int vA = mesh.Cell1DsExtrema(edgeId, 0);
	unsigned int vB = mesh.Cell1DsExtrema(edgeId, 1);
	auto edgeKey = make_pair(min(vA, vB), max(vA, vB));

	auto it = edgeToFacesMap.find(edgeKey);
	if (it == edgeToFacesMap.end())
		return Vector3d::Zero();

	const auto& adjacentFaces = it->second;

	if (adjacentFaces.size() == 2) {
		unsigned int otherFace = (adjacentFaces[0] == currentFaceId) ? adjacentFaces[1] : adjacentFaces[0];
		return getFaceBarycenter(mesh, otherFace);
	} else {
		return Vector3d::Zero(); // bordo o errore
	}
}

// Triangola ogni faccia e salva triangoli in meshTriangulated
void triangulateAndStore2(PolyhedralMesh& meshIn, PolyhedralMesh& meshOut, const vector<int>& dims,
						  map<pair<unsigned int, unsigned int>, vector<unsigned int>> edgeToFacesMap) {

	meshOut.Cell0DsId.resize(dims[0]);
	meshOut.Cell0DsCoordinates = MatrixXd::Zero(3, dims[0]);
	meshOut.Cell0DsFlag.resize(dims[0]);

	meshOut.Cell1DsId.resize(dims[1]);
	meshOut.Cell1DsExtrema = MatrixXi::Zero(dims[1], 2);
	meshOut.Cell1DsFlag.resize(dims[1]);

	meshOut.Cell2DsId.resize(dims[2]);
	meshOut.Cell2DsVertices.resize(dims[2]);
	meshOut.Cell2DsEdges.resize(dims[2]);

	unsigned int vertexCounter = 0;
	unsigned int edgeCounter = 0;
	unsigned int faceCounter = 0;

	for (unsigned int faceId = 0; faceId < meshIn.Cell2DsId.size(); ++faceId) {
		const auto& verts = meshIn.Cell2DsVertices[faceId];
		Vector3d V0 = meshIn.Cell0DsCoordinates.col(verts[0]);
		Vector3d V1 = meshIn.Cell0DsCoordinates.col(verts[1]);
		Vector3d V2 = meshIn.Cell0DsCoordinates.col(verts[2]);
		Vector3d bary = (V0 + V1 + V2) / 3.0;

		const auto& edges = meshIn.Cell2DsEdges[faceId];

		for (unsigned int e = 0; e < 3; e++) {
			unsigned int vA = verts[e];
			unsigned int vB = verts[(e + 1) % 3];

			if (!meshIn.Cell1DsOriginalFlag[edges[e]]) {
				Vector3d midpoint = (meshIn.Cell0DsCoordinates.col(vA) + meshIn.Cell0DsCoordinates.col(vB)) / 2.0;

				unsigned int id_vA = FindAddVertice(meshIn.Cell0DsCoordinates.col(vA), meshOut, vertexCounter);
				unsigned int id_vB = FindAddVertice(meshIn.Cell0DsCoordinates.col(vB), meshOut, vertexCounter);
				unsigned int id_mid = FindAddVertice(midpoint, meshOut, vertexCounter);
				unsigned int id_bary = FindAddVertice(bary, meshOut, vertexCounter);

				// Primo triangolo
				FindAddFace({id_vA, id_mid, id_bary},
							{FindAddEdge2(id_vA, id_mid, meshOut, edgeCounter),
							 FindAddEdge2(id_mid, id_bary, meshOut, edgeCounter),
							 FindAddEdge2(id_bary, id_vA, meshOut, edgeCounter)},
							meshOut, faceCounter);

				// Secondo triangolo
				FindAddFace({id_vB, id_mid, id_bary},
							{FindAddEdge2(id_vB, id_mid, meshOut, edgeCounter),
							 FindAddEdge2(id_mid, id_bary, meshOut, edgeCounter),
							 FindAddEdge2(id_bary, id_vB, meshOut, edgeCounter)},
							meshOut, faceCounter);
			} else {
				Vector3d otherBary = FindNearBarycenter(meshIn, edges[e], faceId, edgeToFacesMap);

				unsigned int id_vA = FindAddVertice(meshIn.Cell0DsCoordinates.col(vA), meshOut, vertexCounter);
				unsigned int id_vB = FindAddVertice(meshIn.Cell0DsCoordinates.col(vB), meshOut, vertexCounter);
				unsigned int id_bary = FindAddVertice(bary, meshOut, vertexCounter);
				unsigned int id_otherBary = FindAddVertice(otherBary, meshOut, vertexCounter);

				unsigned int edge_b1_b2 = FindAddEdge2(id_bary, id_otherBary, meshOut, edgeCounter);

				// Primo triangolo
				FindAddFace({id_vA, id_bary, id_otherBary},
							{FindAddEdge2(id_vA, id_bary, meshOut, edgeCounter),
							 edge_b1_b2,
							 FindAddEdge2(id_otherBary, id_vA, meshOut, edgeCounter)},
							meshOut, faceCounter);

				// Secondo triangolo
				FindAddFace({id_vB, id_bary, id_otherBary},
							{FindAddEdge2(id_vB, id_bary, meshOut, edgeCounter),
							 edge_b1_b2,
							 FindAddEdge2(id_otherBary, id_vB, meshOut, edgeCounter)},
							meshOut, faceCounter);
			}
		}
	}
}