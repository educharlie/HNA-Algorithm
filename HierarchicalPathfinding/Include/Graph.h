#ifndef	__GRAPH_H__
#define	__GRAPH_H__

#include "DetourCommon.h"
#include "DetourAlloc.h"
#include <vector>
#include <map>


const int maxInternalPath = 256;
const int maxNodes = 2048;

static const float H_SCALE = 0.999f;

class Graph{

typedef unsigned int dtPolyRef;

public:
	
	struct IntraEdge
	{
		int startPosId;
		int endPosId;
		int nPath;
		float cost;
		dtPolyRef path[maxInternalPath];
	};

	struct Edge
	{
		dtPolyRef targetNodeId;
		int idPos;
		int idPoly;
		float pos[3];
	};

	struct Node
	{
		dtPolyRef idNode;
		dtPolyRef idParent;

		Edge *edges;
		unsigned int numEdges;

		IntraEdge *intraEdges;
		unsigned int numIntraEdges;

		void InitEdge(int maxEdgesPerNode)
		{
			numEdges = 0;
			edges = 0;
			edges = (Edge*)dtAlloc(sizeof(Edge)*maxEdgesPerNode, DT_ALLOC_PERM);
			memset(edges, 0, sizeof(Edge)*maxEdgesPerNode);
		}

		void DestroyEdge()
		{
			dtFree(edges);
			numEdges = 0;
		}

		void InitIntraEdge()
		{
			int maxIntraEdgesPerNode = numEdges * numEdges;

			if(maxIntraEdgesPerNode > 0)
			{
				numIntraEdges = 0;
				intraEdges = 0;
				intraEdges = (IntraEdge*)dtAlloc(sizeof(IntraEdge)*maxIntraEdgesPerNode, DT_ALLOC_PERM);
				memset(intraEdges, 0, sizeof(IntraEdge)*maxIntraEdgesPerNode);
			}
		}

		void DestroyIntraEdge()
		{
			dtFree(intraEdges);
			numIntraEdges = 0;
		}
	};

	Node* nodes;
	int numNodes;

	void AddEdge(dtPolyRef sourceIdNode, dtPolyRef targetIdNode, float* pos, int idPos, int idPoly)
	{
		Node *node = &nodes[sourceIdNode];
		Edge *edge = 0;

		edge = &node->edges[node->numEdges++];
		edge->targetNodeId = targetIdNode;
		edge->idPos = idPos;
		edge->idPoly = idPoly;
		dtVcopy(edge->pos, pos);
	}

	void AddNode(int nodeId)
	{
		Node *node = 0;
		node = &nodes[numNodes++];
		node->idNode =  nodeId;
	}

	void AddIntraEdge(dtPolyRef clusterId, int startPosId, int endPosId, float cost, dtPolyRef* path, int nPath)
	{
		Node *node = &nodes[clusterId];
		IntraEdge *intraEdge = 0;

		intraEdge = &node->intraEdges[node->numIntraEdges++];

		intraEdge->startPosId = startPosId;
		intraEdge->endPosId = endPosId;
		intraEdge->cost = cost;
		intraEdge->nPath = nPath;
		memcpy(intraEdge->path, path, sizeof(dtPolyRef)*nPath); 
	}

	void AddParent(dtPolyRef nodeId, dtPolyRef idParent)
	{
		Node *node = &nodes[nodeId];
		node->idParent = idParent;
	}

	void Init(int numMaxNodes)
	{
		numNodes = 0;
		nodes = 0;
		nodes = (Node*)dtAlloc(sizeof(Node)*numMaxNodes, DT_ALLOC_PERM);
		memset(nodes, 0, sizeof(Node)*numMaxNodes);
	}

	void InitEdge(int nodeId, int maxEdges)
	{
		Node *node = &nodes[nodeId];
		node->InitEdge(maxEdges);
	}

	void Destroy()
	{
		dtFree(nodes);
		numNodes = 0;
	}
};

#endif	// __GRAPH_H__