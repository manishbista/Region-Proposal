
#include "../header/disjointSet.h"
#include <iostream>
#include <algorithm>
#include <fstream>

bool operator<(const elements &a, const elements &b) 
{
	bool success = false;
	if(a.itsParentRank > b.itsParentRank) success = true;
	else if(a.itsParentRank < b.itsParentRank) success = false;
	else if((a.itsParentRank == b.itsParentRank) && (a.itsParent > b.itsParent)) success = true;
	else if((a.itsParentRank == b.itsParentRank) && (a.itsParent < b.itsParent)) success = false;

  return success;
}

universe::universe(int numElements)
{
 
//initialize elements
 universeElements = new elements[numElements];
 universe::totalNumOfElements = numElements;

 //initialize elements of the universe
 for(int i = 0; i < numElements; i++)
 {
  	universeElements[i].ID = i;
  	universeElements[i].rank = 0;
  	universeElements[i].itsParent = i;		//each node points to itself i.e. is its own parent
  	universeElements[i].size = 1;
	universeElements[i].itsParentRank = 0;
 }
 
}


universe::~universe()
{
delete universeElements;
}

//for an inputNode, which is the ID of element, find the node's topmost parent node ID
int universe::findParentOf(int inputNode) 
{
  int varNode = inputNode;

//repeat until parent of node is itself, traverse up the tree
  while (varNode != universeElements[varNode].itsParent)
    varNode = universeElements[varNode].itsParent;
//after the parent node is found for inputNode, assign values => varNode
  universeElements[inputNode].itsParent = varNode;
  return varNode;
}

//join two nodes, x and y
//the node having higher rank will be the parent of one with lower rank

void universe::joinNodes(int x, int y) 
{

  if (universeElements[x].rank > universeElements[y].rank) {
    universeElements[y].itsParent = x;				//x now points to y
    universeElements[x].size += universeElements[y].size;	//increase size of x
    universeElements[x].rank++;
  } else {							//if x.rank <= y.rank, make y parent of x
    universeElements[x].itsParent = y;
    universeElements[y].size += universeElements[x].size;
    if (universeElements[x].rank == universeElements[y].rank)
      universeElements[y].rank++;
  }
  totalNumOfElements--;						//reduce the number of elements in the universe
}

void universe::assignParentRank(int numVertices)
{

 for(int i = 0; i < numVertices; i++)
	universeElements[i].itsParentRank = universeElements[universeElements[i].itsParent].size;
 
 std::sort(universeElements, universeElements + numVertices);
 
/*
 std::cout<<" node display SORTED"<<std::endl;
 int totalCount = 0, parentID, index, parentRank;

 int prevParent, curParent, prevParentRank;

 prevParent = universeElements[0].itsParent;
 prevParentRank = universeElements[0].itsParentRank;
 int sizeCount = 0;

 for(int i = 0; i < numVertices; i++)
 {
	index = universeElements[i].ID;
	parentID = universeElements[i].itsParent;
	curParent = parentID;
	parentRank = universeElements[i].itsParentRank;

	if(prevParent != curParent)
	{
		if(totalCount < 25)
		{
			std::cout<<" pParent "<<prevParent<<" cParent "<<curParent<<" ppR "<<prevParentRank<<" cpR "<<parentRank<<" s "<<sizeCount<<std::endl;
			sizeCount = 0;
			totalCount++;
			prevParent = curParent;
			prevParentRank = parentRank;
		}
	}
	sizeCount++;	

 } 
*/


}






