
#include "../header/region.h"
#include <iostream>

region::region()
{
	region::ID = 0;
	region::parentID = 0;
	region::sortedID = 0;
	region::pixelSize = 0;
	region::boundingBox = 0;
	region::colorHist = 0;
	region::textureHist = 0;
}

void rect::Display()
{
	std::cout<<std::endl<<" bBox "<<this->x_min<<" "<<this->y_min<<"  "<<this->x_max<<"  "<<this->y_max;

}


region::region(int ID, int parentID,int sortedID, int size, rect* bBox, float** colHist, float** texHist)
{
	region::ID = ID;
	region::parentID = parentID;
	region::sortedID = sortedID;
	region::pixelSize = size;
	region::boundingBox = bBox;
	region::colorHist = colHist;
	region::textureHist = texHist;

}

region region::operator=(region& a)
{
	return region(a.ID, a.parentID, a.sortedID, a.pixelSize, a.boundingBox, a.colorHist, a.textureHist);

}

void region::Display()
{
 std::cout<<std::endl<<" ID "<<this->ID<<" pID "<<this->parentID<<" sID "<<this->sortedID<<" size "<<this->pixelSize;
/* this->boundingBox->Display();
 std::cout<<" Histogram Info Color"<<std::endl;
 for(int i = 0; i < 3; i++)
 {
	for(int j = 0; j < 16; j++)
	std::cout<<" Ch "<<i<<" Bin "<<j<<" val "<<colorHist[i][j]<<std::endl;
 }

 std::cout<<" Histogram Info Texture"<<std::endl;
 for(int i = 0; i < 3; i++)
 {
	for(int j = 0; j < 8; j++)
	std::cout<<" Ch "<<i<<" Bin "<<j<<" val "<<textureHist[i][j]<<std::endl;
 }
*/
}





