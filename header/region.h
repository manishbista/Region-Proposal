
#ifndef REGION_H
#define REGION_H

struct rect
{
 int x_min, y_min;
 int x_max, y_max;
 rect() : x_min(0), y_min(0), x_max(0), y_max(0) {}
 rect(int a, int b, int c, int d) : x_min(a), y_min(b), x_max(c), y_max(d) {}
 rect operator=(rect& a) { return rect(a.x_min, a.y_min, a.x_max, a.y_max); }
 void Display();
};

class region
{
 private:
 int ID;					//segment identifier
 int parentID;					//ID of its parent
 int sortedID;				//ID of region after it has been sorted
 int pixelSize;				//size of the segment in terms of pixels
 rect* boundingBox;		//info for bounding Box
 float** colorHist;			//3 histogram with 16 bins in each of them
 float** textureHist;		//3 histogram with 8 bins in each

 public:
	

	region();
	region(int ID, int parentID, int sortedID, int size, rect* bBox, float** cHist, float** tHist);
	region operator=(region& a);
	void Display();
	rect returnBoundingBox() { return *(this->boundingBox); }
	float** getColorHistogram() { return (this->colorHist); }
	float** getTextureHistogram() { return (this->textureHist); }
	int getPixelSize() { return pixelSize; }
	rect* getBoundingBox() { return (this->boundingBox);}
	int* getParentID() { return &(this->parentID); }
	int* getSortedID() { return &(this->sortedID); }
	int* pixSize() {return &(this->pixelSize);}
};

#endif
