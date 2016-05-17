#include <algorithm>
#include "Graph.h"
#include "DetourNavMesh.h"
#include "Recast.h"
#include "DetourStatus.h"
#include "DetourNavMeshQuery.h"
#include "DetourNode.h"
#include "DetourAssert.h"
#include "RecastAssert.h"
#include <fstream>

extern "C" {
 #include "metis.h"
}

typedef int idxtype;

class Navigation
{

public:

	void createHierarchicalGraph(int p_levels,int p_level,int p_mergedPolys,rcContext* ctx, const dtMeshTile* ptile, const dtNavMesh* pm_navMesh, const dtNavMeshQuery* pm_navQuery, std::map<dtPolyRef, int> &nodesInCluster);

	dtStatus findPathNav(rcContext* ctx, dtPolyRef startRef, dtPolyRef endRef, const float* startPos, const float* endPos, const dtQueryFilter* filter, dtPolyRef* path, int &pathCount, const int maxPath);

private:
	static const int MAX_POLYS = 256;
	const dtMeshTile* tile;
	const dtNavMesh* m_navMesh;
	const dtNavMeshQuery* m_navQuery;
	int numParts;
	dtPolyRef refBase;
	
	int levels;
	int level;
    int maxPolyInNode;

	class dtNodePool* m_nodePool;		///< Pointer to node pool.
	class dtNodeQueue* m_openList;		///< Pointer to open list queue.

	Graph mainGraph;
	Graph currentGraph;
	Graph parentGraph;

	int numTotalEdges;
	int numLevel;

	std::map<dtPolyRef, dtPolyRef> nodeCluster;
	std::multimap<dtPolyRef, dtPolyRef> clusterNode;

	Graph *graphs;
	int numGraphs;

	void mergeNodes();
	void buildHierarchy();
	void buildNodes();
	void buildEdges();

	void init();

	dtStatus findHierarchicalPath(dtPolyRef startRef, dtPolyRef endRef, int  startIdPos, int endIdPos, const float* startPos, const float* endPos, dtPolyRef* tempPathNodes,int* tempPathPolys, int *nTempPath, const int maxPath);
	float findPathLocal(dtPolyRef startRef, dtPolyRef endRef,int  startIdPos, int endIdPos, const float* startPos, const float* endPos, dtPolyRef* path, int &pathCount, const int maxPath,  std::vector<dtPolyRef> subGraphNodes);
	float findPath(dtPolyRef startRef, dtPolyRef endRef, const float* startPos, const float* endPos, dtPolyRef* path, int &pathCount, const int maxPath);

	Graph::Node * getNode(dtPolyRef ref, int l);
	void getPath(int fromPosId, int toPosId, Graph::Node *node, int l, dtPolyRef* tempPath, int &nTempPath);
	void setGraph();
	float getCost(const Graph::Node *node, int startPosId, int endPosId, const float * startPos, const float * endPos);

	void linkStartToGraph(Graph::Node *node, dtPolyRef ref, const float *pos, int startIdPos);
	void linkEndToGraph(Graph::Node *node, dtPolyRef ref, const float *pos, int startIdPos);

	void checkPartition(int* part, const int numNodes, const int numParts);
	void explorePartition(int idNode, int* newPart, int* part);
};

