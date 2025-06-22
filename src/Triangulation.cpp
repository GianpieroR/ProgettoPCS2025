#include "Utils.hpp"
#include "PolyhedralMesh.hpp"
#include <vector>
#include <array>
#include <Eigen/Dense>
#include <limits>

using namespace std;
using namespace Eigen;

namespace PolyhedralLibrary
{

	//============================
	// Triangulation
	//============================
	void Triangulation(int q, int b, int c, PolyhedralMesh& mesh, PolyhedralMesh& meshFinal)
	{
		PolyhedralMesh meshTriangulated;

		vector<int> dimension = ComputePolyhedronVEF(q, b, c);
		vector<int> dimensionDuplicated = CalculateDuplicated(q, b, c, dimension);

		triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
		RemoveDuplicatedVertices(meshTriangulated);
		RemoveDuplicatedEdges(meshTriangulated);
		
		NewMesh(meshTriangulated, meshFinal, dimension);
		PopulateCell3D(meshFinal, dimension);
	}


	//============================
	// Triangulation + Dual
	//============================
	void TriangulationDual(int q, int b, int c, PolyhedralMesh& mesh, PolyhedralMesh& meshDual)
	{
		PolyhedralMesh meshTriangulated;
		PolyhedralMesh meshFinal;

		vector<int> dimension = ComputePolyhedronVEF(q, b, c);
		vector<int> dimensionDuplicated = CalculateDuplicated(q, b, c, dimension);

		triangulateAndStore(mesh, meshTriangulated, b, c, dimensionDuplicated);
		RemoveDuplicatedVertices(meshTriangulated);
		RemoveDuplicatedEdges(meshTriangulated);
		
		NewMesh(meshTriangulated, meshFinal, dimension);

		auto edgeToFacesMap = buildEdgeToFacesMap(meshFinal);

		CalculateDual(meshFinal, meshDual, edgeToFacesMap);
		PopulateCell3D(meshDual, dimension);
	}


	//============================
	// Triangulation2
	//============================
	void Triangulation2(int q, int b, PolyhedralMesh& mesh, PolyhedralMesh& meshTriangulated2)
	{
		PolyhedralMesh meshTriangulated;
		PolyhedralMesh meshFinal;

		vector<int> dimension = ComputePolyhedronVEF(q, b, 0);
		vector<int> dimensionDuplicated = CalculateDuplicated(q, b, 0, dimension);

		triangulateAndStore(mesh, meshTriangulated, b, 0, dimensionDuplicated);
		RemoveDuplicatedVertices(meshTriangulated);
		RemoveDuplicatedEdges(meshTriangulated);
		
		NewMesh(meshTriangulated, meshFinal, dimension);

		vector<int> dimension2 = CalculateDimension2(b, q);
		auto edgeToFacesMap = buildEdgeToFacesMap(meshFinal);

		triangulateAndStore2(meshFinal, meshTriangulated2, dimension2, edgeToFacesMap);
		PopulateCell3D(meshTriangulated2, dimension);
	}


	//============================
	// Triangulation2 + Dual
	//============================
	void Triangulation2Dual(int q, int b, PolyhedralMesh& mesh, PolyhedralMesh& meshDual)
	{
		PolyhedralMesh meshTriangulated;
		PolyhedralMesh meshTriangulated2;
		PolyhedralMesh meshFinal;

		vector<int> dimension = ComputePolyhedronVEF(q, b, 0);
		vector<int> dimensionDuplicated = CalculateDuplicated(q, b, 0, dimension);

		triangulateAndStore(mesh, meshTriangulated, b, 0, dimensionDuplicated);
		RemoveDuplicatedVertices(meshTriangulated);
		RemoveDuplicatedEdges(meshTriangulated);
		
		NewMesh(meshTriangulated, meshFinal, dimension);

		vector<int> dimension2 = CalculateDimension2(b, q);
		auto edgeToFacesMap = buildEdgeToFacesMap(meshFinal);

		triangulateAndStore2(meshFinal, meshTriangulated2, dimension2, edgeToFacesMap);
		CalculateDual(meshTriangulated2, meshDual, edgeToFacesMap);
		PopulateCell3D(meshDual, dimension);
	}


	//============================
	// Helper: Populate Cell3D
	//============================
	void PopulateCell3D(PolyhedralMesh& meshTriangulated, const vector<int>& dimension)
	{
		meshTriangulated.Cell3DsId = {0};
		meshTriangulated.NumCells0Ds = dimension[0];
		meshTriangulated.NumCells1Ds = dimension[1];
		meshTriangulated.NumCells2Ds = dimension[2];

		for (unsigned int i = 0; i < meshTriangulated.Cell0DsId.size(); i++)
			meshTriangulated.Cell3DsVertices.push_back(meshTriangulated.Cell0DsId[i]);

		for (unsigned int i = 0; i < meshTriangulated.Cell1DsId.size(); i++)
			meshTriangulated.Cell3DsEdges.push_back(meshTriangulated.Cell1DsId[i]);

		for (unsigned int i = 0; i < meshTriangulated.Cell2DsId.size(); i++)
			meshTriangulated.Cell3DsFaces.push_back(i);
	}

}