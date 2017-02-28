#version 130
uniform sampler2D smoothImage;
varying vec2 outUV;
uniform float SCREEN_WIDTH;
uniform float SCREEN_HEIGHT;

out vec4 outSquare;
out vec4 outDiamond;

float texelWidthOffset;
float texelHeightOffset;

vec3 center, positionA, positionB, positionC, positionD;
vec4 intensityDiff;

void main()
{
 texelWidthOffset = 1.0/SCREEN_WIDTH;
 texelHeightOffset = 1.0/SCREEN_HEIGHT;

 center = texture2D(smoothImage, outUV).rgb;
 positionA = texture2D(smoothImage, vec2(outUV.x - texelWidthOffset, outUV.y + texelHeightOffset)).rgb;
 positionB = texture2D(smoothImage, vec2(outUV.x + texelWidthOffset, outUV.y + texelHeightOffset)).rgb;
 positionC = texture2D(smoothImage, vec2(outUV.x + texelWidthOffset, outUV.y - texelHeightOffset)).rgb;
 positionD = texture2D(smoothImage, vec2(outUV.x - texelWidthOffset, outUV.y - texelHeightOffset)).rgb;
 intensityDiff.x = distance(center, positionA);
 intensityDiff.y = distance(center, positionB);
 intensityDiff.z = distance(center, positionC);
 intensityDiff.w = distance(center, positionD);
 outSquare = intensityDiff;

 positionA = texture2D(smoothImage, vec2(outUV.x, outUV.y + texelHeightOffset)).rgb;
 positionB = texture2D(smoothImage, vec2(outUV.x + texelWidthOffset, outUV.y)).rgb;
 positionC = texture2D(smoothImage, vec2(outUV.x, outUV.y - texelHeightOffset)).rgb;
 positionD = texture2D(smoothImage, vec2(outUV.x - texelWidthOffset, outUV.y)).rgb;
 intensityDiff.x = distance(center, positionA);
 intensityDiff.y = distance(center, positionB);
 intensityDiff.z = distance(center, positionC);
 intensityDiff.w = distance(center, positionD);
 outDiamond = intensityDiff;
}


