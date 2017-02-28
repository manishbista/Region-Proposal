
#include "../header/segment.h"
#include <algorithm>
#include <cstdlib>

//comparision operation for sorting edges based on weight
bool operator<(const edge &a, const edge &b) {
  return a.weight < b.weight;
}

bool operator<(const similarity &a, const similarity &b)
{
 return a.similarityMeasure < b.similarityMeasure;
}

void edge::DisplayEdge()
{
 std::cout<<"  "<<this->nodeA<<"---"<<this->weight<<"---"<<this->nodeB<<std::endl;
}

segment::segment(int SCREEN_WIDTH, int SCREEN_HEIGHT, int thresConst, int minSize)
{
 thresholdConstant = thresConst;			//scaling factor
 segment::SCREEN_WIDTH = SCREEN_WIDTH;
 segment::SCREEN_HEIGHT = SCREEN_HEIGHT;
 segment::num_edges = 0;				//total number of edges, initially set to zero
 segment::minimumSize = minSize;			//minimum permissible size of a segment

//create a universe with nodes equal to the number of pixels in window
 uni = new universe(SCREEN_WIDTH * SCREEN_HEIGHT);

//initialize threshold values for all nodes
//threshold value measures internal difference in intensity
//weight of an edge measurese intensity differece between boundaries

  threshold = new float[SCREEN_WIDTH * SCREEN_HEIGHT];
  for (int i = 0; i < SCREEN_WIDTH * SCREEN_HEIGHT; i++)
    threshold[i] = calculateThreshold(1,thresholdConstant); //size of node for each pixel is initially 1


 segmentSpace = new unsigned char[SCREEN_WIDTH * SCREEN_HEIGHT * 4];
  for(int i = 0; i < SCREEN_WIDTH * SCREEN_HEIGHT; i++)
 {
 	segmentSpace[4*i + 0] = rand()%255;
 	segmentSpace[4*i + 1] = rand()%255;
 	segmentSpace[4*i + 2] = rand()%255;
 }

//for color histogram memory allocation, you have the following conditions
//there are minimumSize number of segments
//each segment has three channels H, S and V
//each channel requires a histogram of 16 bins

segmentColorHistogram = new float** [minimumSize];
for(int i = 0; i < minimumSize; i++)
	segmentColorHistogram[i] = new float* [3];			//for HSV
for(int i = 0; i < minimumSize; i++)
 {
	for(int j = 0; j < 3; j++)
		segmentColorHistogram[i][j] = new float [16];	//for 16 bins
 }

//for texture histogram memory allocation, you have the following conditions
//there are minimumSize number of segments
//each segment has three channels R, G and B
//each channel requires a histogram of 8 bins

segmentTextureHistogram = new float** [minimumSize];
for(int i = 0; i < minimumSize; i++)
	segmentTextureHistogram[i] = new float* [3];			//for RGB
for(int i = 0; i < minimumSize; i++)
 {
	for(int j = 0; j < 3; j++)
		segmentTextureHistogram[i][j] = new float [8];		//for 8 bins
 }

Regions = new region* [minimumSize];
bBox = new rect* [minimumSize];

//for similarity Measure reserve n(n-1)/2

//an edge is defined between a center node and eight nodes surrounding it
//so nodes corresponding to vertices that lie on the boundary of display are not considered to minimize nested conditional statements
//removing top and bottom boundary, there are (SCREEN_HEIGHT -2) rows of vertices

 pEdge = new edge [SCREEN_WIDTH * SCREEN_HEIGHT * 8];


 for(int i = SCREEN_WIDTH * 1; i < SCREEN_WIDTH * (SCREEN_HEIGHT-1); i++)
 {
  if(i%SCREEN_WIDTH != 0 && (i + 1)%SCREEN_WIDTH != 0)
  {
	pEdge[8 *i +0]. nodeA = i;
	pEdge[8 *i +0]. nodeB = i + SCREEN_WIDTH - 1;

	pEdge[8 *i +1]. nodeA = i;
	pEdge[8 *i +1]. nodeB = i + SCREEN_WIDTH;

	pEdge[8 *i +2]. nodeA = i;
	pEdge[8 *i +2]. nodeB = i + SCREEN_WIDTH + 1;

	pEdge[8 *i +3]. nodeA = i;
	pEdge[8 *i +3]. nodeB = i + 1;

	pEdge[8 *i +4]. nodeA = i;
	pEdge[8 *i +4]. nodeB = i - SCREEN_WIDTH + 1;

	pEdge[8 *i +5]. nodeA = i;
	pEdge[8 *i +5]. nodeB = i - SCREEN_WIDTH;

	pEdge[8 *i +6]. nodeA = i;
	pEdge[8 *i +6]. nodeB = i - SCREEN_WIDTH - 1;

	pEdge[8 *i +7]. nodeA = i;
	pEdge[8 *i +7]. nodeB = i -1;

	num_edges += 8;
  }
 }

}

void segment::fillEdges(float *squareImage, float *diamondImage)
{
 for(int i = 0; i < SCREEN_WIDTH * SCREEN_HEIGHT; i++)
 {
	pEdge[8 *i +0].weight = *(squareImage + 4 * i + 0);

	pEdge[8 *i +1].weight = *(diamondImage + 4 * i + 0);

	pEdge[8 *i +2].weight = *(squareImage + 4 * i + 1);

	pEdge[8 *i +3].weight = *(diamondImage + 4 * i + 1);

	pEdge[8 *i +4].weight = *(squareImage + 4 * i + 2);

	pEdge[8 *i +5].weight = *(diamondImage + 4 * i + 2);

	pEdge[8 *i +6].weight = *(squareImage + 4 * i + 3);

	pEdge[8 *i +7].weight = *(diamondImage + 4 * i + 3);
 }
}


void segment::segmentGraph()
{
 //sort all edges in decreasing order of weights
 std::sort(pEdge, pEdge + num_edges);

 int parentA, parentB;
  for (int i = 0; i < SCREEN_WIDTH * SCREEN_HEIGHT * 8; i++) 
  {
    edge *mEdge = &pEdge[i];
    if(mEdge->nodeA != -1 && mEdge->nodeB != -1)
    {
    // components conected by this edge
    parentA = uni->findParentOf(mEdge->nodeA);
    parentB = uni->findParentOf(mEdge->nodeB);

    if (parentA != parentB) 				//if both of them are not the same node or pixel
      {
	//if external boundary has less variation than within the boundary, join the nodes
       if ((mEdge->weight <= threshold[parentA]) && (mEdge->weight <= threshold[parentB])) 
		{
		uni->joinNodes(parentA, parentB);		//nodes joined
		parentA = uni->findParentOf(parentA);	//find the topmost parent and change its internal weight
		threshold[parentA] = mEdge->weight + calculateThreshold(uni->size(parentA), thresholdConstant);	//weight of new edge and 
        }
      }
    }
  }

}

void segment::fillRegions(unsigned char* hsvImage, unsigned char* magnitudeImage, unsigned char* angleImage)
{

	//var declarations
	int segmentCount = 0, prevParent, parentNode, nodeID, segmentSize;
    int x_min, x_max, y_min, y_max, x_pos, y_pos;

    int numBins_col = 16;
    int histCounter_h[numBins_col], histCounter_s[numBins_col], histCounter_v[numBins_col];    
    int col_h, col_s, col_v;

	int numBins_tex = 8;
    int histCounter_r[numBins_tex], histCounter_g[numBins_tex], histCounter_b[numBins_tex];
	int col_r, col_g, col_b;    

    //initialize histogram counters
    for(int i = 0; i < numBins_col; i++)
    {
        histCounter_h[i] = 0;
        histCounter_s[i] = 0;
        histCounter_v[i] = 0;
    }

	for(int i = 0; i < numBins_tex; i++)
	{
		histCounter_r[i] = 0;
		histCounter_g[i] = 0;
		histCounter_b[i] = 0;
	}


	//assign parentID to all nodes before sorting or any processing functions
	for(int i = 0; i < SCREEN_WIDTH * SCREEN_HEIGHT; i++)
	 uni->findParentOf(i);

	//sort nodes or vertices based on their parent's rank i.e the vertices forming larger
	//groups are ordered first
	uni->assignParentRank(SCREEN_WIDTH * SCREEN_HEIGHT);

	//return parentID of first vertex
	prevParent = uni->returnParentOf(0);
    //initialize extremities of bounding box
    x_min = SCREEN_WIDTH -1;
    x_max = 0;
    y_min = SCREEN_HEIGHT -1;
    y_max = 0;   


 for(int i = 0; i < SCREEN_WIDTH * SCREEN_HEIGHT; i++)
 {
	parentNode = uni->returnParentOf(i);				//parentID of vertices
	nodeID = uni->returnNodeID(i);						//nodeID of vertices

	if(parentNode != prevParent && segmentCount < minimumSize)			
	{

		segmentSize = uni->returnParentRankOf(i-1);
		bBox[segmentCount] = new rect(x_min, y_min, x_max, y_max);

		for(int k = 0; k < numBins_col; k++)
		{
		 segmentColorHistogram[segmentCount][0][k] = (float)histCounter_h[k]/segmentSize;
		 segmentColorHistogram[segmentCount][1][k] = (float)histCounter_s[k]/segmentSize;
		 segmentColorHistogram[segmentCount][2][k] = (float)histCounter_v[k]/segmentSize;
		}

		//to normalize texture histogram, find the overall area of histogram i.e. the sum of discrete values for each channel
        col_r = 0, col_g = 0, col_b = 0;
        for(int k = 0; k < numBins_tex; k++)
        {
         col_r += histCounter_r[k];
		 col_g += histCounter_g[k];
		 col_b += histCounter_b[k];
        }

		for(int k = 0; k < numBins_tex; k++)
		{
		 segmentTextureHistogram[segmentCount][0][k] = (float)histCounter_r[k]/col_r;
		 segmentTextureHistogram[segmentCount][1][k] = (float)histCounter_g[k]/col_g;
		 segmentTextureHistogram[segmentCount][2][k] = (float)histCounter_b[k]/col_b;
		}
		
		Regions[segmentCount] = new region(segmentCount, segmentCount, segmentCount, segmentSize, bBox[segmentCount], segmentColorHistogram[segmentCount], segmentTextureHistogram[segmentCount]);
		//Regions[segmentCount]-s>Display();
	
            //save
   /*      std::cout<<std::endl<<" ID "<<segmentCount;
          std::cout<<" xmin "<<x_min<<" x_max "<<x_max<<" y_min "<<y_min<<" y_max "<<y_max<<std::endl;
            std::cout<<" BBox Area "<<(y_max - y_min)*(x_max - x_min)<<" pixsize "<<uni->returnParentRankOf(i-1);
            std::cout<<std::endl<<" histogram Info "<<std::endl;
            col_h = 0, col_s = 0, col_v = 0;
            for(int k = 0; k < numBins_col; k++)
            {
             std::cout<<"  h  "<<histCounter_h[k]<<"  s  "<<histCounter_s[k]<<"  v  "<<histCounter_v[k]<<std::endl;
             col_h += histCounter_h[k];
			 col_s += histCounter_s[k];
			 col_v += histCounter_v[k];
            }
            std::cout<<" total Hist h: "<<col_h<<" s: "<<col_s<<" v: "<<col_v<<std::endl;
    

            std::cout<<std::endl<<" histogram Info "<<std::endl;
            col_r = 0, col_g = 0, col_b = 0;
            for(int k = 0; k < numBins_tex; k++)
            {
             std::cout<<"  r  "<<histCounter_r[k]<<"  g  "<<histCounter_g[k]<<"  b  "<<histCounter_b[k]<<std::endl;
             col_r += histCounter_r[k];
			 col_g += histCounter_g[k];
			 col_b += histCounter_b[k];
            }
            std::cout<<" total Hist r: "<<col_r<<" s: "<<col_g<<" v: "<<col_b<<std::endl;

        float colr = 0, colg = 0, colb = 0;
		for(int k = 0; k < numBins_tex; k++)
		 {
		  colr += segmentTextureHistogram[segmentCount][0][k];
		  colg += segmentTextureHistogram[segmentCount][1][k];
		  colb += segmentTextureHistogram[segmentCount][2][k];

			std::cout<<std::endl<<" h "<<segmentTextureHistogram[segmentCount][0][k]<<" s "<<segmentTextureHistogram[segmentCount][1][k]<<" v "<<segmentTextureHistogram[segmentCount][2][k];
		 }
			std::cout<<std::endl<<" ID "<<segmentCount<<" ho "<<colr<<" ss "<<colg<<" vo "<<colb;
*/

            //re-initialize or update
			segmentCount++;
			prevParent = parentNode;
            x_min = SCREEN_WIDTH -1;
            x_max = 0;
            y_min = SCREEN_HEIGHT -1;
            y_max = 0;

            for(int k = 0; k < numBins_col; k++)
            {
             histCounter_h[k] = 0;
             histCounter_s[k] = 0;
             histCounter_v[k] = 0;
            }

			for(int k = 0; k < numBins_tex; k++)
			{
			 histCounter_r[k] = 0;
			 histCounter_g[k] = 0;
			 histCounter_b[k] = 0;
			}
	}

    //for a specific number of segments, extract info: color histogram, Texture Gradients, bounding box
    //save those information and update when the end of segment is reached
	if(segmentCount < minimumSize)
	{
       //1. Extract Bounding Box Data ie (x, y)[min, max]
        //1.1 Initial extremities of x and y has been saved
        //1.2 compare x and y values    
        x_pos = nodeID % SCREEN_WIDTH;
        y_pos = nodeID / SCREEN_WIDTH;

        if(x_pos < x_min) x_min = x_pos;
        if(x_pos > x_max) x_max = x_pos;
        if(y_pos < y_min) y_min = y_pos;
        if(y_pos > y_max) y_max = y_pos;

        //2. Extract Color Histogram Data from input hsv image for 16 bins
        col_h = (int)angleImage[4*nodeID + 0];
        col_s = (int)angleImage[4*nodeID + 1];
        col_v = (int)angleImage[4*nodeID + 2];

        histCounter_h[col_h/16]++;
        histCounter_s[col_s/16]++;
        histCounter_v[col_v/16]++;

		//3. Extract Texture Histogram Data of Oriented Gradients from input textures for 8 bins
        col_r = (int)angleImage[4*nodeID + 0];
        col_g = (int)angleImage[4*nodeID + 1];
        col_b = (int)angleImage[4*nodeID + 2];

        histCounter_r[col_r/32] += magnitudeImage[4*nodeID + 0];
        histCounter_g[col_g/32] += magnitudeImage[4*nodeID + 1];
        histCounter_b[col_b/32] += magnitudeImage[4*nodeID + 2];
        

        //assigning colors
   		segmentSpace[4*nodeID + 0] = segmentSpace[4*parentNode + 0];
		segmentSpace[4*nodeID + 1] = segmentSpace[4*parentNode + 1];
		segmentSpace[4*nodeID + 2] = segmentSpace[4*parentNode + 2];
		segmentSpace[4*nodeID + 3] = (255 - segmentCount);

	}
	
	else
	{
		segmentSpace[4*nodeID + 0] = 32;
		segmentSpace[4*nodeID + 1] = 32;
		segmentSpace[4*nodeID + 2] = 32;
		segmentSpace[4*nodeID + 3] = 255;

	}
 }


int numNeighbours = 0;
float simCol;
similarRegions.clear();
similarRegions.reserve((minimumSize * (minimumSize - 1))/2);

	for(int i = 0; i < minimumSize; i++)
	{
		for(int j = i+1; j < minimumSize; j++)
		{
			if(isIntersecting(i, j))
			{
			 simCol = calculateSimilarity(i, j);
			// std::cout<<std::endl<<" regionsF "<<i<<" regionS "<<j<<" simCol "<<simCol;
			 similarRegions.push_back(similarity(i, j, simCol));
			 numNeighbours++;
			}
		}

	}

/*
for(int i = 0; i < similarRegions.size(); i++)
{
std::cout<<" f "<<similarRegions[i].firstRegionID<<" s "<<similarRegions[i].secondRegionID<<" v "<<similarRegions[i].similarityMeasure;
std::cout<<"      ";
if(i % 5 == 0) std::cout<<std::endl;
}
*/

int k = 0;
int fRID, sRID, fPID, sPID;
while(k < minimumSize)
{

 std::vector<similarity>::iterator it = std::max_element(similarRegions.begin(), similarRegions.end());
 fRID = it->firstRegionID;
 sRID = it->secondRegionID;
 fPID = trackParentOf(fRID);
 sPID = trackParentOf(sRID);

 //std::cout<<" itr f "<<fRID<<" fPID "<<fPID<<" s "<<sRID<<" sPID "<<sPID<<" v "<<it->similarityMeasure<<std::endl;
 int mergedRegionID = mergeRegions(fPID, sPID, k);
 int savedRegionID = (mergedRegionID == fPID) ? sPID : fPID;
 //std::cout<<" merged "<<mergedRegionID<<" saved "<<savedRegionID<<std::endl;
 std::cout<<" fID "<<fPID<<" sortID "<<k<<std::endl;

 //update similarity values for mergedRegion pairs
 for(int itr = 0; itr < similarRegions.size(); itr++)
 {
	fRID = similarRegions[itr].firstRegionID;
	sRID = similarRegions[itr].secondRegionID;
	fPID = trackParentOf(fRID);
	sPID = trackParentOf(sRID);

	if(fPID == mergedRegionID)
	{
		if(sPID == savedRegionID || sPID == mergedRegionID)
		 similarRegions[itr].similarityMeasure = 0.0;
		else
		 similarRegions[itr].similarityMeasure = calculateSimilarity(fPID, sPID);
	}
	else if(sPID == mergedRegionID)
	{
		if(fPID == savedRegionID || fPID == mergedRegionID)
		 similarRegions[itr].similarityMeasure = 0.0;
		else
		 similarRegions[itr].similarityMeasure = calculateSimilarity(fPID, sPID);
	}
 }
k++;

}

std::cout<<std::endl;
for(int i = 0; i < minimumSize; i++)
	Regions[i]->Display();

}



int segment::mergeRegions(int firstRegionID, int secondRegionID, int sortCount)
{
//firstRegionID and secondRegionID are regionIDs of topmost parent
int swapID;

 if(Regions[firstRegionID]->getPixelSize() < Regions[secondRegionID]->getPixelSize())
 {
	swapID = firstRegionID;
	firstRegionID = secondRegionID;
	secondRegionID = swapID;
 }

//firstRegion has larger number of pixels
//save first region and update second region as the merged one
 *(Regions[firstRegionID]->getSortedID()) = sortCount;
 *(Regions[firstRegionID]->getParentID()) = secondRegionID;


 rect firstBox = Regions[firstRegionID]->returnBoundingBox();
 rect secondBox = Regions[secondRegionID]->returnBoundingBox();
 Regions[secondRegionID]->getBoundingBox()->x_min = std::min(firstBox.x_min, secondBox.x_min);
 Regions[secondRegionID]->getBoundingBox()->y_min = std::min(firstBox.y_min, secondBox.y_min);
 Regions[secondRegionID]->getBoundingBox()->x_max = std::max(firstBox.x_max, secondBox.x_max);
 Regions[secondRegionID]->getBoundingBox()->y_max = std::max(firstBox.y_max, secondBox.y_max);
	
//Regions[secondRegionID]->getBoundingBox()->Display();

 mergeHistogram(Regions[firstRegionID]->getColorHistogram(), Regions[secondRegionID]->getColorHistogram(), 16, 
				Regions[firstRegionID]->getPixelSize(), Regions[secondRegionID]->getPixelSize());
 mergeHistogram(Regions[firstRegionID]->getTextureHistogram(), Regions[secondRegionID]->getTextureHistogram(), 8, 
				Regions[firstRegionID]->getPixelSize(), Regions[secondRegionID]->getPixelSize());

 *(Regions[secondRegionID]->pixSize()) += (Regions[firstRegionID]->getPixelSize());	
 return secondRegionID;
}

float segment::mergeHistogram(float** firstHist, float** secondHist, int numBins, int firstSize, int secondSize)
{
 for(int i = 0; i < 3; i++)
 {
	for(int j = 0; j < numBins; j++)
		secondHist[i][j] = (float)(firstHist[i][j] * firstSize + secondHist[i][j]* secondSize)/(firstSize + secondSize);
 }
}


int segment::trackParentOf(int regionID) 
{
  int varNode = regionID;

//repeat until parent of node is itself, traverse up the tree
  while (varNode != *(Regions[varNode]->getParentID()))
    varNode = *(Regions[varNode]->getParentID());
//after the parent node is found for inputNode, assign values => varNode
  *(Regions[varNode]->getParentID()) = varNode;
  return varNode;
}

float segment::calcSimOfColor(int firstRegionID, int secondRegionID)
{
 float** firstColHist = Regions[firstRegionID]->getColorHistogram();
 float** secondColHist = Regions[secondRegionID]->getColorHistogram();
 float sum = 0.0;

	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 16; j++)
			sum += std::min(firstColHist[i][j], secondColHist[i][j]);
	}
 return sum/3.0;
}

float segment::calcSimOfTexture(int firstRegionID, int secondRegionID)
{
 float** firstTexHist = Regions[firstRegionID]->getTextureHistogram();
 float** secondTexHist = Regions[secondRegionID]->getTextureHistogram();
 float sum = 0.0;

	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 8; j++)
			sum += std::min(firstTexHist[i][j], secondTexHist[i][j]);
	}
 return sum/3.0;
}

float segment::calcSimOfSize(int firstRegionID, int secondRegionID)
{
 int firstpXSize = Regions[firstRegionID]->getPixelSize();
 int secondpXSize = Regions[secondRegionID]->getPixelSize();
 return (1.0 - (float)(firstpXSize + secondpXSize)/(SCREEN_WIDTH * SCREEN_HEIGHT) );
}


float segment::calcSimOfRect(int firstRegionID, int secondRegionID)
{
 rect firstBox = Regions[firstRegionID]->returnBoundingBox();
 rect secondBox = Regions[secondRegionID]->returnBoundingBox();

 int minX = std::min(firstBox.x_min, secondBox.x_min);
 int minY = std::min(firstBox.y_min, secondBox.y_min);
 int maxX = std::max(firstBox.x_max, secondBox.x_max);
 int maxY = std::max(firstBox.y_max, secondBox.y_max);

 float unionArea = (float)(maxX - minX)*(maxY - minY);
 unionArea -= (Regions[firstRegionID]->getPixelSize() + Regions[secondRegionID]->getPixelSize());
 unionArea /= (SCREEN_WIDTH * SCREEN_HEIGHT);
 return (1 - unionArea);
}

float segment::calculateSimilarity(int firstRegionID, int secondRegionID)
{
float sim = calcSimOfColor(firstRegionID, secondRegionID) + 
			calcSimOfTexture(firstRegionID, secondRegionID) + 
			calcSimOfRect(firstRegionID, secondRegionID)+
			calcSimOfSize(firstRegionID, secondRegionID);

 return sim;
}


bool segment::isIntersecting(int firstRegionID, int secondRegionID)
{
	bool intersectX = false;
	rect firstBox = Regions[firstRegionID]->returnBoundingBox();
	rect secondBox = Regions[secondRegionID]->returnBoundingBox();
	if(secondBox.x_min >= firstBox.x_min && secondBox.x_min < firstBox.x_max) intersectX = true;
	else if(firstBox.x_min >= secondBox.x_min && firstBox.x_min < secondBox.x_max) intersectX = true;
	else return false;

	bool intersectY = false;
	if(secondBox.y_min >= firstBox.y_min && secondBox.y_min < firstBox.y_max) intersectY = true;
	else if(firstBox.y_min >= secondBox.y_min && firstBox.y_min < secondBox.y_max) intersectY = true;
	else return false;

	return true;

}




//draw bounding Box around regions
void segment::drawBoundingBox(int start, int stop, unsigned char* img)
{
	for(int i = 0; i < SCREEN_WIDTH*SCREEN_HEIGHT*4; i++)
	segmentSpace[i] = img[i];

	//x_min, y_min, x_max, y_max
	int vertexID, x_min, x_max, y_min, y_max, ID;
	unsigned char col_r, col_g, col_b;
	for(int i = start; i < stop; i++)
	{
		for(int m = 0; m < segment::minimumSize; m++)
		 {
			if(*(Regions[m]->getSortedID()) == i)
			ID = m;	

		 }

		rect bBox = Regions[ID]->returnBoundingBox();
		x_min = bBox.x_min;
		x_max = bBox.x_max;
		y_min = bBox.y_min;
		y_max = bBox.y_max;
		//std::cout<<std::endl<<" ID "<<i;<<"bBox"<<x_min<<" "<<x_max<<"  "<<y_min<<"  "<<y_max;

		col_r = rand()%255;
		col_g = rand()%255;
		col_b = rand()%255;

		for(int j = x_min; j < x_max; j++)
		{
		 vertexID = y_min * SCREEN_WIDTH + j;
		 segmentSpace[4*vertexID + 0] = col_r;
		 segmentSpace[4*vertexID + 1] = col_g;
		 segmentSpace[4*vertexID + 2] = col_b;
		 segmentSpace[4*vertexID + 3] = 255;

		 vertexID = y_max * SCREEN_WIDTH + j;
		 segmentSpace[4*vertexID + 0] = col_r;
		 segmentSpace[4*vertexID + 1] = col_g;
		 segmentSpace[4*vertexID + 2] = col_b;
		 segmentSpace[4*vertexID + 3] = 255;
		}

		for(int j = y_min; j < y_max; j++)	
		{
		 vertexID = j * SCREEN_WIDTH + x_min;			//(x_min, y_min) to (x_min, y_max) left boundary
		 segmentSpace[4*vertexID + 0] = col_r;
		 segmentSpace[4*vertexID + 1] = col_g;
		 segmentSpace[4*vertexID + 2] = col_b;
		 segmentSpace[4*vertexID + 3] = 255;

		 vertexID = j * SCREEN_WIDTH + x_max;			//(x_max, y_min) to (x_max, y_max) right boundary
		 segmentSpace[4*vertexID + 0] = col_r;
		 segmentSpace[4*vertexID + 1] = col_g;
		 segmentSpace[4*vertexID + 2] = col_b;
		 segmentSpace[4*vertexID + 3] = 255;
		}


	}
}


segment::~segment()
{
 delete threshold;
 delete uni;
 delete pEdge;
 delete[] Regions;
 delete[] bBox;
 delete segmentSpace;

 for(int i = 0; i < minimumSize; i++)
	{
	 delete[] segmentColorHistogram[i];
	 delete[] segmentTextureHistogram[i];
	}
}




