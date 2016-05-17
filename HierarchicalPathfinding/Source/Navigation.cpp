#include "Navigation.h"

void Navigation::init()
{
	//liberate memory
	for(int i = 0; i < numGraphs; i++)
	{
		graphs[i].nodes->DestroyIntraEdge();
		graphs[i].nodes->DestroyEdge();
		graphs[i].Destroy();			
	}
	
	numGraphs = 0;
	//reserve memory por graphs
	graphs = (Graph*)dtAlloc(sizeof(Graph)*levels, DT_ALLOC_PERM);
	memset(graphs, 0, sizeof(Graph)*levels);
	numLevel = 0;

	//reserve memory por pool nodes
	m_nodePool = 0;
	m_openList = 0;
	if (!m_nodePool || m_nodePool->getMaxNodes() < maxNodes)
	{
		if (m_nodePool)
		{
			m_nodePool->~dtNodePool();
			dtFree(m_nodePool);
			m_nodePool = 0;
		}
		m_nodePool = new (dtAlloc(sizeof(dtNodePool), DT_ALLOC_PERM)) dtNodePool(maxNodes, dtNextPow2(maxNodes/4));
	}
	else
	{
		m_nodePool->clear();
	}

	//reserve memory for open list
	if (!m_openList || m_openList->getCapacity() < maxNodes)
	{
		if (m_openList)
		{
			m_openList->~dtNodeQueue();
			dtFree(m_openList);
			m_openList = 0;
		}
		m_openList = new (dtAlloc(sizeof(dtNodeQueue), DT_ALLOC_PERM)) dtNodeQueue(maxNodes);
	}
	else
	{
		m_openList->clear();
	}

	//get the base poly ref
	refBase = m_navMesh->getPolyRefBase(tile);
}

void Navigation::createHierarchicalGraph(int p_levels,int p_level,int p_mergedPolys, rcContext* ctx,const dtMeshTile* ptile,  const dtNavMesh* pm_navMesh, const dtNavMeshQuery* pm_navQuery, std::map<dtPolyRef, int> &nodesInCluster)
{
	rcAssert(ctx);

	if(maxPolyInNode == p_mergedPolys && levels == p_levels)
	{
		level = p_level;
		setGraph();
	}
	else
	{

		ctx->startTimer(RC_TIMER_BUILD_HIERARCHICAL_GRAPH);

		levels = p_levels;
		level = p_level;
		maxPolyInNode = p_mergedPolys;

		tile = ptile;
		m_navMesh = pm_navMesh;
		m_navQuery = pm_navQuery;

		init();

		const dtPoly* poly = 0;
		const dtPoly* neighbourPoly = 0;
		dtPolyRef ref;
		dtPolyRef neighbourRef;
		int numEdgesNode = 0;

		//Main graph
		parentGraph.Init(tile->header->polyCount);

		//create nodes
		for (int i = 0; i < tile->header->polyCount; ++i)
		{
			parentGraph.AddNode(i);
		}

		//reserve memory for edges
		for(int i = 0; i < parentGraph.numNodes; i++)
		{
			ref = parentGraph.nodes[i].idNode + refBase;
			m_navMesh->getTileAndPolyByRefUnsafe(ref, &tile, &poly);

			for (unsigned int j = poly->firstLink; j != DT_NULL_LINK; j = tile->links[j].next)
			{
				numEdgesNode++;
			}

			if(numEdgesNode > 0)
			{
				parentGraph.InitEdge(ref - refBase, numEdgesNode);
				numEdgesNode = 0;
			}
		}

		//create edges
		int idPos = 1;

		std::map<std::pair<dtPolyRef, dtPolyRef>, float *> positions;

		for(int i = 0; i < parentGraph.numNodes; i++)
		{
			//current node
			ref = parentGraph.nodes[i].idNode + refBase;

			m_navMesh->getTileAndPolyByRefUnsafe(ref, &tile, &poly);

			for (unsigned int j = poly->firstLink; j != DT_NULL_LINK; j = tile->links[j].next)
			{
				neighbourRef = tile->links[j].ref;
				m_navMesh->getTileAndPolyByRefUnsafe(neighbourRef, &tile, &neighbourPoly);

				//get position
				float mid[3];
				m_navQuery->getEdgeMidPoint(ref, poly, tile, neighbourRef, neighbourPoly, tile, mid);

				std::pair<dtPolyRef,dtPolyRef> p1;
				std::pair<dtPolyRef,dtPolyRef> p2;

				p1 = std::make_pair(ref,neighbourRef);
				auto it1 = positions.find(p1);

				if(it1 == positions.end())
				{
					p2 = std::make_pair(neighbourRef,ref);

					//add edge
					parentGraph.AddEdge(ref - refBase,neighbourRef - refBase, mid, idPos,ref - refBase);
					parentGraph.AddEdge(neighbourRef- refBase,ref - refBase, mid, idPos,neighbourRef- refBase);

					positions.insert(std::make_pair<std::pair<dtPolyRef,dtPolyRef>,float*>(p1,mid));
					positions.insert(std::make_pair<std::pair<dtPolyRef,dtPolyRef>,float*>(p2,mid));

					idPos++;
				}
			}
		}
		numTotalEdges = idPos;
		positions.clear();

		//hierarchical subdivisions
		mainGraph = parentGraph;
		graphs[numGraphs++] = parentGraph;
		buildHierarchy();

		setGraph();

		ctx->stopTimer(RC_TIMER_BUILD_HIERARCHICAL_GRAPH);

		for(int i = 0; i < numGraphs; i++)
		{
			ctx->log(RC_LOG_PROGRESS, ">>Level %d = %d nodes.", i, graphs[i].numNodes);
		}

		std::ofstream fout("results.txt", std::ios_base::app | std::ios_base::out);
		fout << "\n\nMerged Polys: " << p_mergedPolys << "\n";
		for(int i = 0; i < numGraphs; i++)
		{
			fout << "level " << i << ": " <<  graphs[i].numNodes << "\n";
		}
		fout.close();
	}

	//get painted nodes
	if(level !=0)
	{
		for(int i = 0; i < graphs[0].numNodes; i++)
		{
			dtPolyRef nodeId = graphs[0].nodes[i].idNode;
			Graph::Node* node = getNode(nodeId,0);
			if(node == NULL)
				nodesInCluster.insert(std::make_pair(nodeId + refBase, -1));
			else
				nodesInCluster.insert(std::make_pair(nodeId + refBase, node->idNode));
		}
	}
}

dtStatus Navigation::findPathNav(rcContext* ctx, dtPolyRef startRef, dtPolyRef endRef,const float* startPos, const float* endPos, const dtQueryFilter* filter, dtPolyRef* path, int &pathCount, const int maxPath)
{
	rcAssert(ctx);

	std::ofstream fout("results.txt", std::ios_base::app | std::ios_base::out);
	fout << "\n" << maxPolyInNode << " " << level << " " << currentGraph.numNodes << " ";

	dtPolyRef tempPathNodes[128];
	int tempPathPolys[128];
	int nTempPath = 0;

	ctx->startTimer(RC_TIMER_GET_PARENT_NODES);

	Graph::Node *sNode =  getNode(startRef-refBase,0);
	Graph::Node *eNode =  getNode(endRef-refBase,0);

	if(sNode == NULL || eNode==NULL)
		return DT_FAILURE;

	fout << sNode->numEdges << " " << eNode->numEdges << " ";

	int nSNodeIntraedges = sNode->numIntraEdges;
	int nENodeIntraedges = eNode->numIntraEdges;

	ctx->stopTimer(RC_TIMER_GET_PARENT_NODES);

	fout << ctx->getAccumulatedTime(RC_TIMER_GET_PARENT_NODES)/1000.0f << " ";

	if((level == 0) || (sNode->idNode == eNode->idNode))
	{
		ctx->startTimer(RC_TIMER_FIND_HIERACHICAL_PATH);
		m_navQuery->findPath(startRef, endRef, startPos, endPos, filter, path, &pathCount, MAX_POLYS);
		ctx->stopTimer(RC_TIMER_FIND_HIERACHICAL_PATH);

		fout << ctx->getAccumulatedTime(RC_TIMER_FIND_HIERACHICAL_PATH)/1000.0f << " ";

		return DT_SUCCESS;
	}

	//new start and end position id
	int startIdPos = numTotalEdges+1;
	int endIdPos  = startIdPos+1;

	//add node start to graph
	ctx->startTimer(RC_TIMER_LINK_START_END_NODES);

	linkStartToGraph(sNode,startRef,startPos,startIdPos);
	linkEndToGraph(eNode,endRef,endPos, endIdPos);

	ctx->stopTimer(RC_TIMER_LINK_START_END_NODES);

	fout << ctx->getAccumulatedTime(RC_TIMER_LINK_START_END_NODES)/1000.0f << " ";

	//find path
	ctx->startTimer(RC_TIMER_FIND_HIERACHICAL_PATH);

	startRef = sNode->idNode;
	endRef = eNode->idNode;
	findHierarchicalPath(startRef, endRef, startIdPos, endIdPos, startPos, endPos, tempPathNodes, tempPathPolys, &nTempPath, maxPath);

	ctx->stopTimer(RC_TIMER_FIND_HIERACHICAL_PATH);

	fout << ctx->getAccumulatedTime(RC_TIMER_FIND_HIERACHICAL_PATH)/1000.0f << " ";

	int posId;
	int n = 0;

	//start
	ctx->startTimer(RC_TIMER_GET_SUBPATH);

	if(nTempPath > 1)
	{
		posId = tempPathNodes[1];

		for(unsigned int i = nSNodeIntraedges; i < sNode->numIntraEdges; i++)
		{
			if((sNode->intraEdges[i].startPosId ==  startIdPos) && (sNode->intraEdges[i].endPosId ==  posId))
			{
				memcpy(path, sNode->intraEdges[i].path, sizeof(dtPolyRef)*sNode->intraEdges[i].nPath); 
				n = sNode->intraEdges[i].nPath;
				break;
			}
		}
	}

	//middle
	for(int i = 1; i < nTempPath-1; i++)
	{
		Graph::Node *node = &currentGraph.nodes[tempPathPolys[i]];

		getPath(tempPathNodes[i], tempPathNodes[i+1], node, level, path, n);
	}

	//end
	posId = tempPathNodes[nTempPath-1];

	for(unsigned int i = nENodeIntraedges; i < eNode->numIntraEdges; i++)
	{
		if((eNode->intraEdges[i].startPosId ==  posId) && (eNode->intraEdges[i].endPosId ==  endIdPos))
		{
			memcpy(path + n, eNode->intraEdges[i].path, sizeof(dtPolyRef)*eNode->intraEdges[i].nPath);
			n += eNode->intraEdges[i].nPath;
			break;
		}
	}

	ctx->stopTimer(RC_TIMER_GET_SUBPATH);

	fout << ctx->getAccumulatedTime(RC_TIMER_GET_SUBPATH)/1000.0f << " ";

	pathCount = n;

	//delete start node and links
	ctx->startTimer(RC_TIMER_DELETE_START_END_NODES);

	sNode->numIntraEdges = nSNodeIntraedges;
	eNode->numIntraEdges = nENodeIntraedges;
		
	ctx->stopTimer(RC_TIMER_DELETE_START_END_NODES);

	fout << ctx->getAccumulatedTime(RC_TIMER_DELETE_START_END_NODES)/1000.0f;

	fout.close();
	return DT_SUCCESS;
}

void Navigation::getPath(int fromPosId, int toPosId, Graph::Node *node, int l, dtPolyRef* tempPath, int &nTempPath)
{
	if(l==0)
		return;

	for(unsigned int i = 0; i < node->numIntraEdges; i++)
	{
		if((node->intraEdges[i].startPosId ==  fromPosId) && (node->intraEdges[i].endPosId ==  toPosId))
		{
			Graph::IntraEdge intraEdge =node->intraEdges[i];

			l--;
			if(l==0)
			{
				memcpy(tempPath + nTempPath, intraEdge.path, sizeof(dtPolyRef)*intraEdge.nPath);
				nTempPath += intraEdge.nPath;
			
				return;
			}
			else
			{
				dtPolyRef nextRef;
				int startIdPos;
				int endIdPos;

				startIdPos = fromPosId;
		
				//for each subNode
				int n = intraEdge.nPath;
				for(int j = 0; j < n; j++ )  
				{
					Graph::Node * subNode = &graphs[l].nodes[intraEdge.path[j]];
					nextRef = intraEdge.path[(j+1)%n];

					if(j==n-1)
					{
						endIdPos = toPosId;
					}
					else
					{
						unsigned int numSubEdges =  subNode->numEdges;

						for(unsigned int k = 0; k < numSubEdges; k++ )  
						{
							if(subNode->edges[k].targetNodeId == nextRef)
							{
								endIdPos = subNode->edges[k].idPos;
								break;
							}
						}
					}

					getPath(startIdPos, endIdPos, subNode, l,tempPath, nTempPath);

					startIdPos = endIdPos;
				}

			}
			break;
		}
	}
}

dtStatus Navigation::findHierarchicalPath(dtPolyRef startRef, dtPolyRef endRef, int  startIdPos, int endIdPos, const float* startPos, const float* endPos,  dtPolyRef* tempPathNodes, int* tempPathPolys, int *nTempPath, const int maxPath)
{

	dtAssert(m_nodePool);
	dtAssert(m_openList);

	m_nodePool->clear();
	m_openList->clear();

	dtNode* startNode = m_nodePool->getNode(startIdPos);

	dtVcopy(startNode->pos, startPos);
	startNode->pidx = 0;
	startNode->cost = 0;
	startNode->total = dtVdist(startPos, endPos) * H_SCALE;
	startNode->id = startIdPos;
	startNode->idPos = startRef;
	startNode->flags = DT_NODE_OPEN;
	m_openList->push(startNode);
	
	dtNode* lastBestNode = startNode;
	float lastBestNodeCost = startNode->total;

	dtStatus status = DT_SUCCESS;

	while (!m_openList->empty())
	{
		// Remove node from open list and put it in closed list.
		dtNode* bestNode = m_openList->pop();
		bestNode->flags &= ~DT_NODE_OPEN;
		bestNode->flags |= DT_NODE_CLOSED;
		
		// Reached the goal, stop searching.
		if (bestNode->id == endIdPos)
		{
			lastBestNode = bestNode;
			break;
		}

		const Graph::Node* node = &currentGraph.nodes[bestNode->idPos];

		// Get parent
		dtPolyRef parentRef = -1;

		if (bestNode->pidx)
			parentRef = m_nodePool->getNodeAtIdx(bestNode->pidx)->id;

		for(unsigned int i=0; i < node->numEdges; i++)
		{
			dtPolyRef neighbourRef = node->edges[i].idPos;

			if(bestNode->id ==neighbourRef)
				continue;

			if (neighbourRef == parentRef)
				continue;

			dtNode* neighbourNode = m_nodePool->getNode(neighbourRef);

			if(!neighbourNode)
				continue;

			// If the node is visited the first time, calculate node position.
			if (neighbourNode->flags == 0)
			{
				dtVcopy(neighbourNode->pos,  node->edges[i].pos);
			}

			// Calculate cost and heuristic.
			float cost = 0.0f;
			float heuristic = 0.0f;

			// Special case for last node.
			if (neighbourRef == endIdPos)
			{
				const float curCost = getCost(node, bestNode->id, neighbourNode->id, bestNode->pos, neighbourNode->pos);
				//const Graph::Node* eNode = currentGraph.nodes[endRef];
				//const float endCost = getCost(eNode, neighbourNode->id, endIdPos, neighbourNode->pos, endPos);
				const float endCost = getCost(node, neighbourNode->id, endIdPos, neighbourNode->pos, endPos);

				cost = bestNode->cost + curCost + endCost;
				heuristic = 0.0f;
			}
			else
			{

				const float curCost = getCost(node, bestNode->id, neighbourNode->id,bestNode->pos, neighbourNode->pos);
				
				cost = bestNode->cost + curCost;
				heuristic = dtVdist(neighbourNode->pos, endPos)*H_SCALE;
			}

			const float total = cost + heuristic;
			
			// The node is already in open list and the new result is worse, skip.
			if ((neighbourNode->flags & DT_NODE_OPEN) && total >= neighbourNode->total)
				continue;
			// The node is already visited and process, and the new result is worse, skip.
			if ((neighbourNode->flags & DT_NODE_CLOSED) && total >= neighbourNode->total)
				continue;
			
			// Add or update the node.
			neighbourNode->pidx = m_nodePool->getNodeIdx(bestNode);
			neighbourNode->idPos = node->edges[i].targetNodeId;
			neighbourNode->id = neighbourRef;
			neighbourNode->flags = (neighbourNode->flags & ~DT_NODE_CLOSED);
			neighbourNode->cost = cost;
			neighbourNode->total = total;
			
			if (neighbourNode->flags & DT_NODE_OPEN)
			{
				// Already in open, update node location.
				m_openList->modify(neighbourNode);
			}
			else
			{
				// Put the node in open list.
				neighbourNode->flags |= DT_NODE_OPEN;
				m_openList->push(neighbourNode);
			}
			
			// Update nearest node to target so far.
			if (heuristic < lastBestNodeCost)
			{
				lastBestNodeCost = heuristic;
				lastBestNode = neighbourNode;
			}		
		}
	}
	
	// Reverse the path.
	dtNode* prev = 0;
	dtNode* node = lastBestNode;
	do
	{
		dtNode* next = m_nodePool->getNodeAtIdx(node->pidx);
		node->pidx = m_nodePool->getNodeIdx(prev);
		prev = node;
		node = next;
	}
	while (node);
	
	// Store path

	node = prev;
	int n = 0;
	do
	{
	
		tempPathNodes[n] = node->id;
		tempPathPolys[n] = node->idPos;
		n++;

		if (n >= maxPath)
		{
			break;
		}
		node = m_nodePool->getNodeAtIdx(node->pidx);
	}
	while (node);
	
	*nTempPath = n;
	return status;
}

void Navigation::linkStartToGraph(Graph::Node *node, dtPolyRef startRef, const float *startPos, int startIdPos)
{
	for(unsigned int i=0; i < node->numEdges; i++)
	{
		dtPolyRef m_polys[MAX_POLYS];
		int m_npolys;	

		float cost = findPath(startRef - refBase, node->edges[i].idPoly, startPos, node->edges[i].pos,m_polys,m_npolys,MAX_POLYS);

		currentGraph.AddIntraEdge(node->idNode, startIdPos, node->edges[i].idPos, cost, m_polys,m_npolys);
	}
}

void Navigation::linkEndToGraph(Graph::Node *node, dtPolyRef startRef, const float *startPos, int startIdPos)
{
	for(unsigned int i=0; i < node->numEdges; i++)
	{
		dtPolyRef m_polys[MAX_POLYS];
		int m_npolys;	

		float cost = findPath(node->edges[i].idPoly, startRef - refBase, node->edges[i].pos, startPos,m_polys,m_npolys,MAX_POLYS);

		currentGraph.AddIntraEdge(node->idNode, node->edges[i].idPos, startIdPos, cost, m_polys,m_npolys);
	}
}

float Navigation::getCost(const Graph::Node *node, int startPosId, int endPosId, const float * startPos, const float * endPos)
{
	if(startPosId == endPosId )
		return 0.0f;
	
	/*for(unsigned int i = 0; i < node->numIntraEdges; i++)
	{
		if((node->intraEdges[i].startPosId ==  startPosId) && (node->intraEdges[i].endPosId ==  endPosId))
		{
			return node->intraEdges[i].cost;
		}
	}*/
		
	return dtVdist(startPos, endPos);
}

void Navigation::buildHierarchy()
{
	do{
		//merge nodes

		//numParts =	floor((float(parentGraph.numNodes) / float(maxPolyInNode)) +0.5f);
		numParts = 	parentGraph.numNodes / maxPolyInNode;

		if(numParts < 2 || numParts >= parentGraph.numNodes)
			break;

		//merge nodes and build cluster
		mergeNodes();

		//create subgraph
		buildNodes();
		buildEdges();

		graphs[numGraphs++] = currentGraph;
		
		parentGraph.nodes = currentGraph.nodes;
		parentGraph.numNodes = currentGraph.numNodes;

		nodeCluster.clear();
		clusterNode.clear();
		numLevel++;

	}while(numLevel < levels);
}

void Navigation::mergeNodes()
{
	int maxLinkCount = 0;
	
	//number edges
	for(int j = 0; j < parentGraph.numNodes; j++)
	{
		Graph::Node * node = &parentGraph.nodes[j];
		for(unsigned int i=0; i < node->numEdges; i++)
		{
			maxLinkCount++;
		}
	}

	idxtype* xadj = (idxtype*)dtAlloc( (parentGraph.numNodes+1)*sizeof(idxtype), DT_ALLOC_TEMP );
	idxtype* adjncy = (idxtype*)dtAlloc( maxLinkCount*2*sizeof(idxtype), DT_ALLOC_TEMP );
	idxtype* part = (idxtype*)dtAlloc( parentGraph.numNodes * sizeof(idxtype), DT_ALLOC_TEMP );

	//get parameters
	int xID = 0;
	int adjID = 0;
	unsigned int ip;

	for(int j = 0; j < parentGraph.numNodes; j++)
	{
		xadj[xID++] = adjID;
		Graph::Node * node = &parentGraph.nodes[j];
		for(unsigned int i=0; i < node->numEdges; i++)
		{
			ip = node->edges[i].targetNodeId;
			adjncy[adjID++] = ip;
		}
	}	
	xadj[xID] = adjID;

	idxtype wgtflag = 0;
	idxtype numflag = 0;
	int options[5] = {0,0,0,0,0};
	idxtype edgecut = -1;

	//Metis
	METIS_PartGraphKway(&parentGraph.numNodes, xadj, adjncy, NULL, NULL, &wgtflag, &numflag, &numParts, options, &edgecut, part);

	//get partitions (check merged nodes are neighbours
	checkPartition(part, parentGraph.numNodes, numParts);

	//set the partition number for each node
	int partCount = 0;

	for(int j = 0; j < parentGraph.numNodes; ++j)
	{
		int idC = part[partCount++];

		Graph::Node *node = &parentGraph.nodes[j];

		if(!node)
			continue;

		if(node->numEdges <= 0)
			continue;

		clusterNode.insert(std::make_pair<dtPolyRef, dtPolyRef>(idC, j));
	}

	//add parent for each child
	int nodeId = 0;

	for(auto it = clusterNode.begin(), end = clusterNode.end(); it != end; it = clusterNode.upper_bound(it->first))
	{
		auto ret = clusterNode.equal_range(it->first);
		for (auto it1=ret.first; it1!=ret.second; ++it1)
		{
			nodeCluster.insert(std::make_pair<dtPolyRef, dtPolyRef>(it1->second, nodeId));
			parentGraph.AddParent(it1->second, nodeId);
		}
		nodeId++;
	}

	//create memory for nodes in graph
	clusterNode.clear();
	currentGraph.Init(nodeId);

	for(auto it = nodeCluster.begin(); it != nodeCluster.end(); ++it ) 
	{
		clusterNode.insert(std::make_pair<dtPolyRef, dtPolyRef>(it->second, it->first));	
	}

	dtFree(xadj);
	dtFree(adjncy);
	dtFree(part);
}

void Navigation::buildNodes()
{
	//nodes
	for(auto it = clusterNode.begin(), end = clusterNode.end(); it != end; it = clusterNode.upper_bound(it->first))
	{
		currentGraph.AddNode(it->first);
	}
}

void Navigation::buildEdges()
{
	std::multimap<std::pair<dtPolyRef,dtPolyRef>, std::pair<int ,float*>> portals;
	std::vector<dtPolyRef> subgraphNodes;

	dtPolyRef clusterId = 0;
	dtPolyRef neighbourClusterId = 0;

	
	// create memory for edges
	int numEdgesNode = 0;
	for(auto it = clusterNode.begin(), end = clusterNode.end(); it != end; it = clusterNode.upper_bound(it->first))
	{
		clusterId = it->first;
		auto ret = clusterNode.equal_range(clusterId);
		for (auto it1=ret.first; it1!=ret.second; ++it1)
		{
			Graph::Node * node = &parentGraph.nodes[it1->second];
			for(unsigned int i=0; i < node->numEdges; i++)
			{
				Graph::Node *neighbour =  &parentGraph.nodes[node->edges[i].targetNodeId];
				neighbourClusterId = nodeCluster.find(neighbour->idNode)->second;
			
				if(neighbourClusterId == clusterId)
					continue;

				Graph::Node *tempNode = &graphs[numLevel].nodes[node->idNode];

				for(unsigned int j=0; j < tempNode->numEdges; j++)
				{
					if(tempNode->edges[j].idPos == node->edges[i].idPos)
					{
						numEdgesNode++;
					}
				}
			}			
		}

		if(numEdgesNode > 0)
		{
			currentGraph.InitEdge(clusterId, numEdgesNode);
			numEdgesNode = 0;
		}
	}

	
	//create edges
	for(auto it = clusterNode.begin(), end = clusterNode.end(); it != end; it = clusterNode.upper_bound(it->first))
	{
		clusterId = it->first;

		//group of nodes
		auto ret = clusterNode.equal_range(clusterId);

		//for each node in the cluster
		for (auto it1=ret.first; it1!=ret.second; ++it1)
		{
			Graph::Node * node = &parentGraph.nodes[it1->second];

			//detect portal nodes between two clusters
			for(unsigned int i=0; i < node->numEdges; i++)
			{
				Graph::Node *neighbour =  &parentGraph.nodes[node->edges[i].targetNodeId];
				neighbourClusterId = nodeCluster.find(neighbour->idNode)->second;

				//add subgraph nodes of thecluster
				if(std::find(subgraphNodes.begin(), subgraphNodes.end(), node->edges[i].idPos)==subgraphNodes.end())
				{
					subgraphNodes.push_back(node->edges[i].idPos);
				}
			
				if(neighbourClusterId == clusterId)
					continue;

				//add edge
				Graph::Node *tempNode = &graphs[numLevel].nodes[node->idNode];

				for(unsigned int j=0; j < tempNode->numEdges; j++)
				{
					if(tempNode->edges[j].idPos == node->edges[i].idPos)
					{
						currentGraph.AddEdge(clusterId,neighbourClusterId, node->edges[i].pos, node->edges[i].idPos,tempNode->edges[i].idPoly);
						break;
					}
				}

				//portals
				std::pair<dtPolyRef,dtPolyRef> p = std::make_pair(node->idNode, neighbour->idNode);
				std::pair<int,float*> p1 = std::make_pair(node->edges[i].idPos, node->edges[i].pos);
				portals.insert(std::make_pair(p, p1));		
			}			
		}

		Graph::Node *currentNode = &currentGraph.nodes[clusterId];
		currentNode->InitIntraEdge();

		//add intraedges
		for(auto it1 = portals.begin(); it1 != portals.end(); ++it1 ) 
		{
			std::pair<dtPolyRef,dtPolyRef> p = it1->first;

			dtPolyRef idNode = p.first;

			std::pair<int,float *> p1 = it1->second;

			float nodePos[3];
			int idPos = p1.first;
			dtVcopy(nodePos, p1.second);

			//intraedges
			for(auto it2 = portals.begin(); it2 != portals.end(); ++it2 ) 
			{
				if(it1==it2)
					continue;

				std::pair<dtPolyRef,dtPolyRef> p2 = it2->first;
				dtPolyRef idNeighbour = p2.first;
			
				float neighbourPos[3];
				std::pair<int,float*> p3 = it2->second;
				int neiughbourIdPos = p3.first;
				dtVcopy(neighbourPos, p3.second);

				dtPolyRef m_polys[MAX_POLYS];
				int m_npolys;

				float cost = findPathLocal(idNode, idNeighbour, idPos, neiughbourIdPos, nodePos,neighbourPos,m_polys,m_npolys,MAX_POLYS,subgraphNodes);
				
				currentGraph.AddIntraEdge(clusterId, idPos, neiughbourIdPos, cost, m_polys, m_npolys);
			}
		}

		subgraphNodes.clear();
		portals.clear();
	}
}

float Navigation::findPath(dtPolyRef startRef, dtPolyRef endRef,const float* startPos, const float* endPos, dtPolyRef* path, int &pathCount, const int maxPath)
{
	float cost = 0.0f;
	dtAssert(m_nodePool);
	dtAssert(m_openList);
	
	pathCount = 0;
	
	m_nodePool->clear();
	m_openList->clear();
	
	dtNode* startNode = m_nodePool->getNode(startRef);
	dtVcopy(startNode->pos, startPos);
	startNode->pidx = 0;
	startNode->cost = 0;
	startNode->total = dtVdist(startPos, endPos) * H_SCALE;
	startNode->id = startRef;
	startNode->flags = DT_NODE_OPEN;
	m_openList->push(startNode);
	
	dtNode* lastBestNode = startNode;
	float lastBestNodeCost = startNode->total;
	
	while (!m_openList->empty())
	{
		// Remove node from open list and put it in closed list.
		dtNode* bestNode = m_openList->pop();
		bestNode->flags &= ~DT_NODE_OPEN;
		bestNode->flags |= DT_NODE_CLOSED;
		
		// Reached the goal, stop searching.
		if (bestNode->id == endRef)
		{
			lastBestNode = bestNode;
			break;
		}
		
		Graph::Node* node = &mainGraph.nodes[bestNode->id];
		
		// Get parent poly and tile.
		dtPolyRef parentRef = 0;
	
		if (bestNode->pidx)
			parentRef = m_nodePool->getNodeAtIdx(bestNode->pidx)->id;
		
		for (unsigned int i=0; i < node->numEdges; i++)
		{
			dtPolyRef neighbourRef = node->edges[i].targetNodeId;

			// Skip invalid ids and do not expand back to where we came from.
			if (!neighbourRef || neighbourRef == parentRef)
				continue;

			dtNode* neighbourNode = m_nodePool->getNode(neighbourRef);
			if (!neighbourNode)
			{
				continue;
			}
			
			// If the node is visited the first time, calculate node position.
			if (neighbourNode->flags == 0)
			{
				dtVcopy(neighbourNode->pos,  node->edges[i].pos);
			}

			// Calculate cost and heuristic.
			float cost = 0.0f;
			float heuristic = 0.0f;
			
			if (neighbourRef == endRef)
			{
				const float curCost = dtVdist(bestNode->pos, neighbourNode->pos);
				const float endCost = dtVdist(bestNode->pos, endPos);

				cost = bestNode->cost + curCost + endCost;
				heuristic = 0.0f;
			}
			else
			{
				const float curCost = dtVdist(bestNode->pos, neighbourNode->pos);
				cost = bestNode->cost + curCost;
				heuristic = dtVdist(neighbourNode->pos, endPos)*H_SCALE;
			}
			
			const float total = cost + heuristic;
			
			// The node is already in open list and the new result is worse, skip.
			if ((neighbourNode->flags & DT_NODE_OPEN) && total >= neighbourNode->total)
				continue;
			// The node is already visited and process, and the new result is worse, skip.
			if ((neighbourNode->flags & DT_NODE_CLOSED) && total >= neighbourNode->total)
				continue;
			
			// Add or update the node.
			neighbourNode->pidx = m_nodePool->getNodeIdx(bestNode);
			neighbourNode->id = neighbourRef;
			neighbourNode->flags = (neighbourNode->flags & ~DT_NODE_CLOSED);
			neighbourNode->cost = cost;
			neighbourNode->total = total;
			
			if (neighbourNode->flags & DT_NODE_OPEN)
			{
				// Already in open, update node location.
				m_openList->modify(neighbourNode);
			}
			else
			{
				// Put the node in open list.
				neighbourNode->flags |= DT_NODE_OPEN;
				m_openList->push(neighbourNode);
			}
			
			// Update nearest node to target so far.
			if (heuristic < lastBestNodeCost)
			{
				lastBestNodeCost = heuristic;
				lastBestNode = neighbourNode;
			}
		}
	}
		
	// Reverse the path.
	dtNode* prev = 0;
	dtNode* node = lastBestNode;
	cost = lastBestNode->cost;
	do
	{
		dtNode* next = m_nodePool->getNodeAtIdx(node->pidx);
		node->pidx = m_nodePool->getNodeIdx(prev);
		prev = node;
		node = next;
	}
	while (node);
	
	// Store path
	node = prev;
	int n = 0;
	do
	{
		path[n++] = node->id + refBase;
		if (n >= maxPath)
		{
			break;
		}
		node = m_nodePool->getNodeAtIdx(node->pidx);
	}
	while (node);
	
	pathCount = n;
	
	return cost;
}

float Navigation::findPathLocal(dtPolyRef startRef, dtPolyRef endRef,int  startIdPos, int endIdPos, const float* startPos, const float* endPos, dtPolyRef* path, int &pathCount, const int maxPath, std::vector<dtPolyRef> subGraphNodes)
{
    dtAssert(m_nodePool);
	dtAssert(m_openList);

	float cost = 0.0f;
	pathCount = 0;
	
	if (!maxPath)
		return cost;

	m_nodePool->clear();
	m_openList->clear();

	//get start and goal node
	startRef = parentGraph.nodes[startRef].idNode;
	endRef = parentGraph.nodes[endRef].idNode;

	if (startRef == endRef)
	{
		if(numLevel == 0)
			path[0] = startRef + refBase;
		else 
			path[0] = startRef;
		pathCount = 1;
		cost = getCost(&parentGraph.nodes[startRef], startIdPos,endIdPos, startPos, endPos);
		return cost;
	}

	dtNode* startNode = m_nodePool->getNode(startIdPos);
	dtVcopy(startNode->pos, startPos);
	startNode->pidx = 0;
	startNode->cost = 0;
	startNode->total = dtVdist(startPos, endPos) * H_SCALE;
	startNode->id = startIdPos;
	startNode->idPos = startRef;
	startNode->flags = DT_NODE_OPEN;
	m_openList->push(startNode);
	
	dtNode* lastBestNode = startNode;
	float lastBestNodeCost = startNode->total;

	while (!m_openList->empty())
	{
		// Remove node from open list and put it in closed list.
		dtNode* bestNode = m_openList->pop();
		bestNode->flags &= ~DT_NODE_OPEN;
		bestNode->flags |= DT_NODE_CLOSED;
		
		// Reached the goal, stop searching.
		if (bestNode->id == endIdPos)
		{
			lastBestNode = bestNode;
			break;
		}

		Graph::Node* node = &parentGraph.nodes[bestNode->idPos];

		// Get parent poly and tile.
		dtPolyRef parentRef = 0;

		if (bestNode->pidx)
			parentRef = m_nodePool->getNodeAtIdx(bestNode->pidx)->id;

		for(unsigned int i=0; i < node->numEdges; i++)
		{
			dtPolyRef neighbourRef = node->edges[i].idPos;

			if(!subGraphNodes.empty())
			{
				if(std::find(subGraphNodes.begin(), subGraphNodes.end(), neighbourRef)==subGraphNodes.end())
					continue;
			}
		
			if (neighbourRef == parentRef)
				continue;

			dtNode* neighbourNode = m_nodePool->getNode(neighbourRef);

			if(bestNode->id ==neighbourNode->id)
				continue;

			// If the node is visited the first time, calculate node position.
			if (neighbourNode->flags == 0)
			{
				dtVcopy(neighbourNode->pos,  node->edges[i].pos);
			}

			// Calculate cost and heuristic.
			float heuristic = 0.0f;

			// Special case for last node.
			if (neighbourRef == endIdPos)
			{

				const float curCost = getCost(node, bestNode->id, neighbourNode->id, bestNode->pos, neighbourNode->pos);
				Graph::Node* eNode = &parentGraph.nodes[endRef];
				const float endCost = getCost(eNode, neighbourNode->id, endIdPos, neighbourNode->pos, endPos);

				cost = bestNode->cost + curCost + endCost;
				heuristic = 0.0f;
			}
			else
			{
				//const float curCost = dtVdist(bestNode->pos, neighbourNode->pos);

				const float curCost = getCost(node, bestNode->id, neighbourNode->id,bestNode->pos, neighbourNode->pos);

				cost = bestNode->cost + curCost;
				heuristic = dtVdist(neighbourNode->pos, endPos)*H_SCALE;
			}

			const float total = cost + heuristic;
			
			// The node is already in open list and the new result is worse, skip.
			if ((neighbourNode->flags & DT_NODE_OPEN) && total >= neighbourNode->total)
				continue;
			// The node is already visited and process, and the new result is worse, skip.
			if ((neighbourNode->flags & DT_NODE_CLOSED) && total >= neighbourNode->total)
				continue;
			
			// Add or update the node.
			neighbourNode->pidx = m_nodePool->getNodeIdx(bestNode);
			neighbourNode->idPos = node->edges[i].targetNodeId;
			neighbourNode->id = neighbourRef;
			neighbourNode->flags = (neighbourNode->flags & ~DT_NODE_CLOSED);
			neighbourNode->cost = cost;
			neighbourNode->total = total;
			
			if (neighbourNode->flags & DT_NODE_OPEN)
			{
				// Already in open, update node location.
				m_openList->modify(neighbourNode);
			}
			else
			{
				// Put the node in open list.
				neighbourNode->flags |= DT_NODE_OPEN;
				m_openList->push(neighbourNode);
			}
			
			// Update nearest node to target so far.
			if (heuristic < lastBestNodeCost)
			{
				lastBestNodeCost = heuristic;
				lastBestNode = neighbourNode;
			}
		}

	}
	
	// Reverse the path.
	dtNode* prev = 0;
	dtNode* node = lastBestNode;
	cost = lastBestNode->cost;
	do
	{
		dtNode* next = m_nodePool->getNodeAtIdx(node->pidx);
		node->pidx = m_nodePool->getNodeIdx(prev);
		prev = node;
		node = next;
	}
	while (node);
	
	// Store path
	node = prev;
	int n = 0;
	do
	{
		if(numLevel == 0)
			path[n++] = node->idPos + refBase;
		else
			path[n++] = node->idPos;
		node = m_nodePool->getNodeAtIdx(node->pidx);
	}
	while (node);
	
	pathCount = n-1;

	return cost;
}

void Navigation::setGraph()
{
	if( level < numGraphs)
	{
		currentGraph =graphs[level];
	}
	else
	{
		currentGraph =graphs[0];
	}
}

Graph::Node *Navigation::getNode(dtPolyRef ref, int l)
{
	if(l == level)
		return &graphs[l].nodes[ref];

	Graph::Node *node = &graphs[l].nodes[ref];
	
	if(node->edges <= 0)
		return NULL;

	l++;
	getNode(node->idParent,l);
}

void Navigation::checkPartition(int* part, const int numNodes, const int numParts)
{
	int* newPart = (int*) dtAlloc( sizeof(int)*numNodes, DT_ALLOC_TEMP );
	bool* used = (bool*) dtAlloc( sizeof(bool)*numParts, DT_ALLOC_TEMP );
	int nextPart = numParts;

	memset(used, 0, sizeof(bool)*numParts );
	for(int i=0; i < numNodes; i++)
		newPart[i] = -1;

	for(int i=0; i < numNodes; i++)
	{
		if( newPart[i] != -1 )
			continue;

		if( !used[ part[i] ] )
		{
			newPart[i] = part[i];
			used[ part[i] ] = true;			
		}
		else
		{
			newPart[i] = nextPart++;
		}
		explorePartition(i, newPart, part);
	}

	memcpy(part, newPart, sizeof(int)*numNodes );
	dtFree(newPart );
	dtFree( used );
}

void Navigation::explorePartition(int idNode, int* newPart, int* part)
{
	Graph::Node *node = &parentGraph.nodes[idNode];
	
	for(unsigned int i=0; i < node->numEdges; i++)
	{
		dtPolyRef neighbourRef = node->edges[i].targetNodeId;

		Graph::Node *neighbour = &parentGraph.nodes[neighbourRef];

		int neighbourID = neighbour->idNode;

		if( neighbourID == -1 )
			continue;

		if( part[idNode] != part[ neighbourID ] )
			continue;

		if( newPart[neighbourID] != -1 )
			continue;

		newPart[neighbourID] = newPart[idNode];
		explorePartition(neighbourID, newPart, part);
	}
}
