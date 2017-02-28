#include <GL/glew.h>
#include <SDL2/SDL_opengl.h>
#include "../header/meshLoader.h"
#include "../header/shader.h"
#include "../header/matrices.h"
#include "../header/segment.h"
#include <iostream>
#include <cmath>
#include <time.h>

using namespace std;

float SCREEN_WIDTH = 640.0;
float SCREEN_HEIGHT = 480.0;
SDL_Window* gWindow;
SDL_Surface* gScreenSurface;
matrices pipeline;
meshLoader* scene;
mesh* quad;

shader* convShades;
shader* poolingShades;
shader* grayShades;
shader* displayShades;
shader* differenceShades;
shader* hsvShades;
shader* orientShades;
segment* segmentationGraph;

unsigned int colorFBO, colorImage;
unsigned int grayFBO, grayImage;
unsigned int convFBO, convImage;
unsigned int poolingFBO, poolingImage;
unsigned int differenceIntensityFBO, differenceSquareImage, differenceDiamondImage;
unsigned int orientFBO, orientMagnitudeImage, orientAngleImage;
unsigned int hsvFBO, hsvImage;
int fbStatus;

glm::mat3 kernelMatrix;
unsigned char magnitudeBuffer[640*480*4], angleBuffer[640*480*4];
unsigned char colorBuffer[640*480*4];			//buffer for original Image		
float pixelBuffer[2][640*480*4];		//pixelBuffer for differenceImage
int segmentThresholdConst = 40;			//lesser value gives oversegmentation
int segmentMinSize = 25;			//number of segments preferred

unsigned int createTexture(int,int);
unsigned int createTexture(int, int, unsigned char*);

unsigned int segmentTexture;
unsigned char* segments;

void init()
{
	pipeline.matrixMode(PROJECTION_MATRIX);
	pipeline.loadIdentity();
	pipeline.ortho(-1.0, 1.0, -1.0, 1.0, 1, 100);
	convShades = new shader("../v_shader/convolutionShader.vs","../f_shader/convolutionShader.frag");
	poolingShades = new shader("../v_shader/poolingShader.vs","../f_shader/poolingShader.frag");
	grayShades = new shader("../v_shader/grayShader.vs","../f_shader/grayShader.frag");
	displayShades = new shader("../v_shader/displayShader.vs","../f_shader/displayShader.frag");
	differenceShades = new shader("../v_shader/differenceShader.vs","../f_shader/differenceShader.frag");
	hsvShades = new shader("../v_shader/hsvShader.vs","../f_shader/hsvShader.frag");
	orientShades = new shader("../v_shader/orientShader.vs","../f_shader/orientShader.frag");

	scene = new meshLoader();

	segmentationGraph = new segment((int)SCREEN_WIDTH, (int)SCREEN_HEIGHT, segmentThresholdConst, segmentMinSize);

	colorImage = createTexture(SCREEN_WIDTH, SCREEN_HEIGHT);
	grayImage = createTexture(SCREEN_WIDTH, SCREEN_HEIGHT);
	convImage = createTexture(SCREEN_WIDTH, SCREEN_HEIGHT);
	poolingImage = createTexture(SCREEN_WIDTH/2.0, SCREEN_HEIGHT/2.0);
	differenceSquareImage = createTexture(SCREEN_WIDTH, SCREEN_HEIGHT);
	differenceDiamondImage = createTexture(SCREEN_WIDTH, SCREEN_HEIGHT);
	hsvImage = createTexture(SCREEN_WIDTH, SCREEN_HEIGHT);
	orientMagnitudeImage = createTexture(SCREEN_WIDTH, SCREEN_HEIGHT);
	orientAngleImage = createTexture(SCREEN_WIDTH, SCREEN_HEIGHT);

	glGenFramebuffers(1, &colorFBO);
	glEnable(GL_TEXTURE_2D);
	glBindFramebuffer(GL_FRAMEBUFFER,colorFBO);
	glFramebufferTexture2D(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT0,GL_TEXTURE_2D,colorImage,0);

	fbStatus = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if(fbStatus != GL_FRAMEBUFFER_COMPLETE)
		std::cout << "Framebuffer is not OK, status=" << fbStatus << std::endl;
	glBindFramebuffer(GL_FRAMEBUFFER,0);


	glGenFramebuffers(1, &grayFBO);
	glEnable(GL_TEXTURE_2D);
	glBindFramebuffer(GL_FRAMEBUFFER,grayFBO);
	glFramebufferTexture2D(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT0,GL_TEXTURE_2D,grayImage,0);

	fbStatus = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if(fbStatus != GL_FRAMEBUFFER_COMPLETE)
		std::cout << "Framebuffer is not OK, status=" << fbStatus << std::endl;
	glBindFramebuffer(GL_FRAMEBUFFER,0);


	glGenFramebuffers(1, &convFBO);
	glEnable(GL_TEXTURE_2D);
	glBindFramebuffer(GL_FRAMEBUFFER,convFBO);
	glFramebufferTexture2D(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT0,GL_TEXTURE_2D,convImage,0);

	fbStatus = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if(fbStatus != GL_FRAMEBUFFER_COMPLETE)
		std::cout << "Framebuffer is not OK, status=" << fbStatus << std::endl;
	glBindFramebuffer(GL_FRAMEBUFFER,0);

	glGenFramebuffers(1, &poolingFBO);
	glEnable(GL_TEXTURE_2D);
	glBindFramebuffer(GL_FRAMEBUFFER,poolingFBO);
	glFramebufferTexture2D(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT0,GL_TEXTURE_2D,poolingImage,0);

	fbStatus = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if(fbStatus != GL_FRAMEBUFFER_COMPLETE)
		std::cout << "Framebuffer is not OK, status=" << fbStatus << std::endl;
	glBindFramebuffer(GL_FRAMEBUFFER,0);


	glGenFramebuffers(1, &differenceIntensityFBO);
	glEnable(GL_TEXTURE_2D);
	glBindFramebuffer(GL_FRAMEBUFFER, differenceIntensityFBO);
	glFramebufferTexture2D(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT0,GL_TEXTURE_2D,differenceSquareImage,0);
	glFramebufferTexture2D(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT1,GL_TEXTURE_2D,differenceDiamondImage,0);
	GLenum bufs[2] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};
	glDrawBuffers(2, bufs);

	fbStatus = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if(fbStatus != GL_FRAMEBUFFER_COMPLETE)
		std::cout << "Framebuffer is not OK, status=" << fbStatus << std::endl;
	glBindFramebuffer(GL_FRAMEBUFFER,0);

	glGenFramebuffers(1, &orientFBO);
	glEnable(GL_TEXTURE_2D);
	glBindFramebuffer(GL_FRAMEBUFFER, orientFBO);
	glFramebufferTexture2D(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT0,GL_TEXTURE_2D,orientMagnitudeImage,0);
	glFramebufferTexture2D(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT1,GL_TEXTURE_2D,orientAngleImage,0);
	GLenum buf[2] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};
	glDrawBuffers(2, buf);

	fbStatus = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if(fbStatus != GL_FRAMEBUFFER_COMPLETE)
		std::cout << "Framebuffer is not OK, status=" << fbStatus << std::endl;
	glBindFramebuffer(GL_FRAMEBUFFER,0);


	glGenFramebuffers(1, &hsvFBO);
	glEnable(GL_TEXTURE_2D);
	glBindFramebuffer(GL_FRAMEBUFFER, hsvFBO);
	glFramebufferTexture2D(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT0,GL_TEXTURE_2D, hsvImage,0);

	fbStatus = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if(fbStatus != GL_FRAMEBUFFER_COMPLETE)
		std::cout << "Framebuffer is not OK, status=" << fbStatus << std::endl;
	glBindFramebuffer(GL_FRAMEBUFFER,0);

	{
		std::vector<unsigned int> indices;
		std::vector<vertexData> vertices;
		vertexData tmp;
		//1.
		tmp.position.change(-1.0,1.0,-1.0);
		tmp.U=0;
		tmp.V=1;
		vertices.push_back(tmp);
		//2.
		tmp.position.change(-1.0,-1.0,-1.0);
		tmp.U=0;
		tmp.V=0;
		vertices.push_back(tmp);
		//3.
		tmp.position.change(1.0,-1.0,-1.0);
		tmp.U=1;
		tmp.V=0;
		vertices.push_back(tmp);
		//4.
		tmp.position.change(1.0,1.0,-1.0);
		tmp.U=1;
		tmp.V=1;
		vertices.push_back(tmp);
		
		indices.push_back(0);
		indices.push_back(1);
		indices.push_back(2);		
		
		indices.push_back(0);
		indices.push_back(2);
		indices.push_back(3);
		quad=new mesh(&vertices,&indices);
	}

//gaussian blur
	kernelMatrix[0][0]=1.0/16.0;	kernelMatrix[0][1]=1.0/8.0;	kernelMatrix[0][2]=1.0/16.0;
	kernelMatrix[1][0]=1.0/8.0;	kernelMatrix[1][1]=1.0/4.0;	kernelMatrix[1][2]=1.0/8.0;
	kernelMatrix[2][0]=1.0/16.0;	kernelMatrix[2][1]=1.0/8.0;	kernelMatrix[2][2]=1.0/16.0;
/*
//sharpen
	kernelMatrix[0][0]=0.0;	kernelMatrix[0][1]=-1.0;	kernelMatrix[0][2]=0.0;
	kernelMatrix[1][0]=-1.0;	kernelMatrix[1][1]=5.00;	kernelMatrix[1][2]=-1.0;
	kernelMatrix[2][0]=0.0;	kernelMatrix[2][1]=-1.0;	kernelMatrix[2][2]=0.0;

//edge detection
	kernelMatrix[0][0]=-1.0;	kernelMatrix[0][1]=-1.0;	kernelMatrix[0][2]=-1.0;
	kernelMatrix[1][0]=-1.0;	kernelMatrix[1][1]=8.00;	kernelMatrix[1][2]=-1.0;
	kernelMatrix[2][0]=-1.0;	kernelMatrix[2][1]=-1.0;	kernelMatrix[2][2]=-1.0;

//identity

	kernelMatrix[0][0]=0.0;	kernelMatrix[0][1]=0.0;	kernelMatrix[0][2]=0.0;
	kernelMatrix[1][0]=0.0;	kernelMatrix[1][1]=1.0;	kernelMatrix[1][2]=0.0;
	kernelMatrix[2][0]=0.0;	kernelMatrix[2][1]=0.0;	kernelMatrix[2][2]=0.0;
*/

	//draw Original Scene
	glClearColor(0.25, 0.25, 0.25, 1);
	pipeline.matrixMode(MODEL_MATRIX);
	glBindFramebuffer(GL_FRAMEBUFFER,colorFBO);
		displayShades->useShader();
		glClear(GL_COLOR_BUFFER_BIT);
		pipeline.updateMatrices(displayShades->getProgramId());
		scene->draw(displayShades->getProgramId());
		glBindTexture(GL_TEXTURE_2D, 0);
		displayShades->delShader();
	glBindFramebuffer(GL_FRAMEBUFFER,0);

	//take input image and convert it to grayscale image of required size

	glBindFramebuffer(GL_FRAMEBUFFER,grayFBO);
		grayShades->useShader();
		glClear(GL_COLOR_BUFFER_BIT);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, colorImage);
		glUniform1i(glGetUniformLocation(grayShades->getProgramId(),"texture0"), 0);
		pipeline.updateMatrices(grayShades->getProgramId());
		quad->draw(grayShades->getProgramId());
		glBindTexture(GL_TEXTURE_2D, 0);
		grayShades->delShader();
	glBindFramebuffer(GL_FRAMEBUFFER,0);


	glBindFramebuffer(GL_FRAMEBUFFER,hsvFBO);
		hsvShades->useShader();
		glClear(GL_COLOR_BUFFER_BIT);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, colorImage);
		glUniform1i(glGetUniformLocation(hsvShades->getProgramId(),"rgbImage"), 0);
		pipeline.updateMatrices(hsvShades->getProgramId());
		quad->draw(hsvShades->getProgramId());
		glBindTexture(GL_TEXTURE_2D, 0);		
		hsvShades->delShader();
	glBindFramebuffer(GL_FRAMEBUFFER,0);


	//take grayscale image and convolve
	glBindFramebuffer(GL_FRAMEBUFFER,convFBO);
		convShades->useShader();
		glClear(GL_COLOR_BUFFER_BIT);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, grayImage);

		glUniform1i(glGetUniformLocation(convShades->getProgramId(),"grayInputImage"),0);
		glUniform1f(glGetUniformLocation(convShades->getProgramId(),"SCREEN_WIDTH"), SCREEN_WIDTH);
		glUniform1f(glGetUniformLocation(convShades->getProgramId(),"SCREEN_HEIGHT"), SCREEN_HEIGHT);


		glUniformMatrix3fv(glGetUniformLocation(convShades->getProgramId(),"kernelMatrix"),1,GL_FALSE,&kernelMatrix[0][0]);

		pipeline.updateMatrices(convShades->getProgramId());
		quad->draw(convShades->getProgramId());
		glBindTexture(GL_TEXTURE_2D, 0);		
		convShades->delShader();
	glBindFramebuffer(GL_FRAMEBUFFER,0);

	//take color image and extract Oriented Gradients
	glBindFramebuffer(GL_FRAMEBUFFER,orientFBO);
		orientShades->useShader();
		glClear(GL_COLOR_BUFFER_BIT);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, colorImage);

		glUniform1i(glGetUniformLocation(orientShades->getProgramId(),"colorImage"),0);
		glUniform1f(glGetUniformLocation(orientShades->getProgramId(),"SCREEN_WIDTH"), SCREEN_WIDTH);
		glUniform1f(glGetUniformLocation(orientShades->getProgramId(),"SCREEN_HEIGHT"), SCREEN_HEIGHT);
		glBindFragDataLocation(orientShades->getProgramId(), 0, "outMagnitude");
		glBindFragDataLocation(orientShades->getProgramId(), 1, "outAngle");

		pipeline.updateMatrices(orientShades->getProgramId());
		quad->draw(orientShades->getProgramId());
		glBindTexture(GL_TEXTURE_2D, 0);		
		orientShades->delShader();
	glBindFramebuffer(GL_FRAMEBUFFER,0);


	//take smoothed image and encode intensity difference in two output textures
	glBindFramebuffer(GL_FRAMEBUFFER,differenceIntensityFBO);
		differenceShades->useShader();
		glClear(GL_COLOR_BUFFER_BIT);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, convImage);

		glUniform1i(glGetUniformLocation(differenceShades->getProgramId(),"smoothImage"),0);
		glUniform1f(glGetUniformLocation(differenceShades->getProgramId(),"SCREEN_WIDTH"), SCREEN_WIDTH);
		glUniform1f(glGetUniformLocation(differenceShades->getProgramId(),"SCREEN_HEIGHT"), SCREEN_HEIGHT);

		glBindFragDataLocation(differenceShades->getProgramId(), 0, "outSquare");
		glBindFragDataLocation(differenceShades->getProgramId(), 1, "outDiamond");

		pipeline.updateMatrices(differenceShades->getProgramId());
		quad->draw(differenceShades->getProgramId());
		glBindTexture(GL_TEXTURE_2D, 0);		
		differenceShades->delShader();
	glBindFramebuffer(GL_FRAMEBUFFER,0);



	//read from frame buffer
	glBindFramebuffer(GL_FRAMEBUFFER, hsvFBO);
	glReadBuffer(GL_COLOR_ATTACHMENT0);
	glReadPixels(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT, GL_RGBA, GL_UNSIGNED_BYTE, &colorBuffer[0]);	
	glReadBuffer(GL_BACK);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	glBindFramebuffer(GL_FRAMEBUFFER, orientFBO);
	glReadBuffer(GL_COLOR_ATTACHMENT0);
	glReadPixels(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT, GL_RGBA, GL_UNSIGNED_BYTE, &magnitudeBuffer[0]);
	glReadBuffer(GL_COLOR_ATTACHMENT1);
	glReadPixels(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT, GL_RGBA, GL_UNSIGNED_BYTE, &angleBuffer[0]);	
	glReadBuffer(GL_BACK);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	glBindFramebuffer(GL_FRAMEBUFFER, differenceIntensityFBO);
	for(int m = 0; m < 2; m++)
	 {
		glReadBuffer(GL_COLOR_ATTACHMENT0 + m);
		glReadPixels(0, 0, SCREEN_WIDTH, SCREEN_HEIGHT, GL_RGBA, GL_FLOAT, &pixelBuffer[m][0]);	
	 }
	glReadBuffer(GL_BACK);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

//segmentation of image
clock_t timer;
timer = clock();

	segmentationGraph->fillEdges(&(pixelBuffer[0][0]), &(pixelBuffer[1][0]));
	segmentationGraph->segmentGraph();
	segmentationGraph->fillRegions(colorBuffer, magnitudeBuffer, angleBuffer);

timer = clock() - timer;
std::cout<<std::endl<<" "<<segmentMinSize<<" segments found in "<<((float)timer)/CLOCKS_PER_SEC<<" seconds "<<std::endl;
	 segments = segmentationGraph->getSegmentedSpace();
 
	 segmentTexture = createTexture(SCREEN_WIDTH, SCREEN_HEIGHT, segments);

}


int bbStart = 0;
int bbStop = 0;

int pStart = 0, pStop = 0;

void display()
{
	if((pStart != bbStart || pStop != bbStop) && (bbStart < segmentMinSize && bbStop < segmentMinSize))
	{
	 std::cout<<std::endl<<" bbStart "<<bbStart<<" bbStop "<<bbStop<<std::endl;
 	 segmentationGraph->drawBoundingBox(bbStart, bbStop, colorBuffer);
	 segments = segmentationGraph->getSegmentedSpace();
 
	 segmentTexture = createTexture(SCREEN_WIDTH, SCREEN_HEIGHT, segments);
	 pStart = bbStart;
	 pStop = bbStop;
	}


/*
	//use convolved and RELued image, to do max-pooling
	glBindFramebuffer(GL_FRAMEBUFFER,poolingFBO);
		poolingShades->useShader();
		glClear(GL_COLOR_BUFFER_BIT);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, convImage);

		glUniform1f(glGetUniformLocation(poolingShades->getProgramId(),"SCREEN_WIDTH"), SCREEN_WIDTH);
		glUniform1f(glGetUniformLocation(poolingShades->getProgramId(),"SCREEN_HEIGHT"), SCREEN_HEIGHT);

		glUniform1i(glGetUniformLocation(poolingShades->getProgramId(),"textureUnit"),0);
		pipeline.updateMatrices(poolingShades->getProgramId());
		quad->draw(poolingShades->getProgramId());

		glBindTexture(GL_TEXTURE_2D, 0);
		poolingShades->delShader();
	glBindFramebuffer(GL_FRAMEBUFFER,0);
*/

	//display ouput, less changes on produced results
	glClearColor(1, 0, 0, 1);
	displayShades->useShader();
	glClear(GL_COLOR_BUFFER_BIT);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, segmentTexture);
		glUniform1i(glGetUniformLocation(displayShades->getProgramId(),"texture0"),0);
	pipeline.updateMatrices(displayShades->getProgramId());
	quad->draw(displayShades->getProgramId());
	displayShades->delShader();

}


int main()
{

	SDL_Init(SDL_INIT_EVERYTHING);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
	SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);

	gWindow = SDL_CreateWindow("SDL_COLLIDE", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH, SCREEN_HEIGHT, SDL_WINDOW_SHOWN | SDL_WINDOW_OPENGL);
	SDL_GLContext gContext = SDL_GL_CreateContext(gWindow);
	glewExperimental = GL_TRUE;
	glewInit();
	SDL_GL_SetSwapInterval( 1 );
	gScreenSurface = SDL_GetWindowSurface( gWindow );

	bool running=true;
	SDL_Event event;	
	init();

	while(running)
	{
		while(SDL_PollEvent(&event))
		{
			switch(event.type)
			{
				case SDL_QUIT:
				running = false;
				break;
	
				case SDL_KEYDOWN:
				switch(event.key.keysym.sym)
					{
						case SDLK_ESCAPE:
							running=false;
							break;
						case SDLK_RIGHT:	
							bbStart++;
							break;
						case SDLK_UP:
							bbStop++;
							break;
						
					}
	
			}
		}

		display();
		SDL_GL_SwapWindow(gWindow);

	}

	
	delete hsvShades;
	delete grayShades;
	delete convShades;
	delete poolingShades;
	delete displayShades;
	delete differenceShades;
	delete orientShades;
	delete segmentationGraph;
	delete scene;
	delete quad;

	SDL_FreeSurface(gScreenSurface);
	SDL_GL_DeleteContext(gContext);
	SDL_DestroyWindow(gWindow);
	SDL_Quit();
	return 0;
}

unsigned int createTexture(int w,int h)
{
	unsigned int textureId;
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1,&textureId);
	glBindTexture(GL_TEXTURE_2D,textureId);
	glTexImage2D(GL_TEXTURE_2D,0, GL_RGBA8, w, h, 0, GL_RGBA, GL_FLOAT, NULL);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_BORDER);
	
	int i;
	i=glGetError();
	if(i!=0)
		std::cout << "Error happened while loading the texture: " << gluErrorString(i) << std::endl;
	glBindTexture(GL_TEXTURE_2D,0);
	return textureId;
}

unsigned int createTexture(int w, int h, unsigned char *pixels)
{
	unsigned int textureId;
	glEnable(GL_TEXTURE_2D);
	glGenTextures(1,&textureId);
	glBindTexture(GL_TEXTURE_2D,textureId);
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA8,w,h,0,GL_RGBA,GL_UNSIGNED_BYTE,pixels);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	
	//int i;
	//i=glGetError();
	//if(i!=0)
	//	std::cout << "Error happened while loading the texture: " << gluErrorString(i) << std::endl;
	glBindTexture(GL_TEXTURE_2D,0);
	return textureId;

}
