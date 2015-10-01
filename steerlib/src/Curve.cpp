//
// Copyright (c) 2015 Mahyar Khayatkhoei
// Copyright (c) 2009-2014 Shawn Singh, Glen Berseth, Mubbasir Kapadia, Petros Faloutsos, Glenn Reinman
// See license.txt for complete license.
//

#include <algorithm>
#include <vector>
#include <util/Geometry.h>
#include <util/Curve.h>
#include <util/Color.h>
#include <util/DrawLib.h>
#include "Globals.h"

using namespace Util;

static bool bRet = false;
bool compareTime(Util::CurvePoint a, Util::CurvePoint b)
{
	return (a.time < b.time);
}

Curve::Curve(const CurvePoint& startPoint, int curveType) : type(curveType)
{
	controlPoints.push_back(startPoint);
}

Curve::Curve(const std::vector<CurvePoint>& inputPoints, int curveType) : type(curveType)
{
	controlPoints = inputPoints;
	sortControlPoints();
}

// Add one control point to the vector controlPoints
void Curve::addControlPoint(const CurvePoint& inputPoint)
{
	controlPoints.push_back(inputPoint);
	sortControlPoints();
}

// Add a vector of control points to the vector controlPoints
void Curve::addControlPoints(const std::vector<CurvePoint>& inputPoints)
{
	for (int i = 0; i < inputPoints.size(); i++)
		controlPoints.push_back(inputPoints[i]);
	sortControlPoints();
}

// Draw the curve shape on screen, usign window as step size (bigger window: less accurate shape)
void Curve::drawCurve(Color curveColor, float curveThickness, int window)
{
#ifdef ENABLE_GUI
	Util::Point outPoint, previousPoint;
	bRet = true;
	int index = 0;
	if (type == catmullCurve)
	{
		for (float i = 0.0; i < 1.0f; i += 0.001)
		{
			outPoint = useCatmullCurve(index, i);
			DrawLib::drawLine(previousPoint, outPoint, Util::Color(1.0, 0.0, 0.0), curveThickness);
			previousPoint = outPoint;
		}
	}
	if (type == hermiteCurve)
	{
		for (index = 0; index < controlPoints.size() - 1; ++index)
		{
			for (float i = 0.0; i < 1.0f; i += 0.001)
			{
				outPoint = useHermiteCurve(index, i);
				DrawLib::drawLine(previousPoint, outPoint, Util::Color(1.0, 0.0, 0.0), curveThickness);
				previousPoint = outPoint;
			}
		}
	}
	return;
#endif
}

// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints()
{

	std::sort(controlPoints.begin()+1, controlPoints.end(), &compareTime);

	return;
}

// Calculate the position on curve corresponding to the given time, outputPoint is the resulting position
bool Curve::calculatePoint(Point& outputPoint, float time)
{
	// Robustness: make sure there is at least two control point: start and end points
	if (!checkRobust())
		return false;
	bRet = false;
	// Define temporary parameters for calculation
	unsigned int nextPoint;
	float normalTime, intervalTime;
	// Find the current interval in time, supposing that controlPoints is sorted (sorting is done whenever control points are added)
	if (!findTimeInterval(nextPoint, time))
		return false;

	// Calculate position at t = time on curve
	if (type == hermiteCurve)
	{
		outputPoint = useHermiteCurve(nextPoint, time);
	}
	else if (type == catmullCurve)
	{
		outputPoint = useCatmullCurve(nextPoint, time);
	}

	// Return
	return true;

}

// Check Roboustness
bool Curve::checkRobust()
{

	if (controlPoints.size() < 2)
		return false;

	return true;
}

// Find the current time interval (i.e. index of the next control point to follow according to current time)
bool Curve::findTimeInterval(unsigned int& nextPoint, float time)
{
	//std::cout << "5" << std::endl;
	static int index = 0;
	if(index < controlPoints.size()-1)
	{
		if (time > controlPoints[index+1].time && type == hermiteCurve)
			++index;
		nextPoint = index;
		return true;
	}
	return false;
}

// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	float normalTime, intervalTime;
	float nTime;
	if (nextPoint >= controlPoints.size() - 1)
		return newPosition;
	if (!bRet)
	{
		nTime = (time-controlPoints[nextPoint].time) / (controlPoints[nextPoint+1].time - controlPoints[nextPoint].time);
	}
	else
	{
		nTime = time;
	}
	
	float timeSquare = nTime*nTime;
	float timeCube = timeSquare*nTime;
	
	newPosition.x =(((2*timeCube) -  (3*timeSquare) + 1)*controlPoints[nextPoint].position.x) +
		(((3*timeSquare) -(2*timeCube))*controlPoints[nextPoint +1].position.x) +
		((timeCube - (2*timeSquare) + nTime)*controlPoints[nextPoint].tangent.x) + 
		((timeCube - timeSquare)*controlPoints[nextPoint +1].tangent.x);
	newPosition.y = (((2 * timeCube) - (3 * timeSquare) + 1)*controlPoints[nextPoint].position.y) +
		(((3 * timeSquare) - (2 * timeCube))*controlPoints[nextPoint + 1].position.y) +
		((timeCube - (2 * timeSquare) + nTime)*controlPoints[nextPoint].tangent.y) + 
		((timeCube - timeSquare)*controlPoints[nextPoint + 1].tangent.y);
	newPosition.z = (((2 * timeCube) - (3 * timeSquare) + 1)*controlPoints[nextPoint].position.z) +
		(((3 * timeSquare) - (2 * timeCube))*controlPoints[nextPoint + 1].position.z) +
		((timeCube - (2 * timeSquare) + nTime)*controlPoints[nextPoint].tangent.z) + 
		((timeCube - timeSquare)*controlPoints[nextPoint + 1].tangent.z);
	
	return newPosition;
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;

	float nTime;
	if (nextPoint > 0)
		return newPosition;
	if (!bRet)
	{
		nTime = (time - controlPoints[nextPoint].time) / (controlPoints[nextPoint + controlPoints.size()-1].time - controlPoints[nextPoint].time);
	}
	else
	{
		nTime = time;
	}
	float timeSquare = nTime*nTime;
	float timeCube = timeSquare*nTime;
	
	/*newPosition.x = .5 * ((2 * controlPoints[nextPoint + 1].position.x) +
		(nTime * (-1 * controlPoints[nextPoint].position.x + controlPoints[nextPoint + 2].position.x)) +
		(timeSquare * (2 * controlPoints[nextPoint].position.x - 5 * controlPoints[nextPoint + 1].position.x + 4 * controlPoints[nextPoint + 2].position.x - controlPoints[nextPoint + 3].position.x)) +
		(timeCube * (-1 * controlPoints[nextPoint].position.x + 3 * controlPoints[nextPoint + 1].position.x - 3 * controlPoints[nextPoint + 2].position.x + controlPoints[nextPoint + 3].position.x)));

	newPosition.y = .5 * ((2 * controlPoints[nextPoint + 1].position.y) +
		(nTime * (-1 * controlPoints[nextPoint].position.y + controlPoints[nextPoint + 2].position.y)) +
		(timeSquare * (2 * controlPoints[nextPoint].position.y - 5 * controlPoints[nextPoint + 1].position.y + 4 * controlPoints[nextPoint + 2].position.y - controlPoints[nextPoint + 3].position.y)) +
		(timeCube * (-1 * controlPoints[nextPoint].position.y + 3 * controlPoints[nextPoint + 1].position.y - 3 * controlPoints[nextPoint + 2].position.y + controlPoints[nextPoint + 3].position.y)));

	newPosition.z = .5 * ((2 * controlPoints[nextPoint + 1].position.z) +
		(nTime * (-1 * controlPoints[nextPoint].position.z + controlPoints[nextPoint + 2].position.z)) +
		(timeSquare * (2 * controlPoints[nextPoint].position.z - 5 * controlPoints[nextPoint + 1].position.z + 4 * controlPoints[nextPoint + 2].position.z - controlPoints[nextPoint + 3].position.z)) +
		(timeCube * (-1 * controlPoints[nextPoint].position.z + 3 * controlPoints[nextPoint + 1].position.z - 3 * controlPoints[nextPoint + 2].position.z + controlPoints[nextPoint + 3].position.z)));
*/
	newPosition.x = 0.5*((2 * controlPoints[nextPoint + 1].position.x) + (-1 * controlPoints[nextPoint].position.x + controlPoints[nextPoint + 2].position.x)*nTime + (2 * controlPoints[nextPoint].position.x - 5 * controlPoints[nextPoint + 1].position.x + 4 * controlPoints[nextPoint + 2].position.x - controlPoints[nextPoint + 3].position.x)*timeSquare + (controlPoints[nextPoint + 3].position.x - 3 * controlPoints[nextPoint + 2].position.x + 3 * controlPoints[nextPoint + 1].position.x - controlPoints[nextPoint].position.x)*timeCube);
	newPosition.y = 0.5*((2 * controlPoints[nextPoint + 1].position.y) + (-1 * controlPoints[nextPoint].position.y + controlPoints[nextPoint + 2].position.y)*nTime + (2 * controlPoints[nextPoint].position.y - 5 * controlPoints[nextPoint + 1].position.y + 4 * controlPoints[nextPoint + 2].position.y - controlPoints[nextPoint + 3].position.y)*timeSquare + (controlPoints[nextPoint + 3].position.y - 3 * controlPoints[nextPoint + 2].position.y + 3 * controlPoints[nextPoint + 1].position.y - controlPoints[nextPoint].position.y)*timeCube);
	newPosition.z = 0.5*((2 * controlPoints[nextPoint + 1].position.z) + (-1 * controlPoints[nextPoint].position.z + controlPoints[nextPoint + 2].position.z)*nTime + (2 * controlPoints[nextPoint].position.z - 5 * controlPoints[nextPoint + 1].position.z + 4 * controlPoints[nextPoint + 2].position.z - controlPoints[nextPoint + 3].position.z)*timeSquare + (controlPoints[nextPoint + 3].position.z - 3 * controlPoints[nextPoint + 2].position.z + 3 * controlPoints[nextPoint + 1].position.z - controlPoints[nextPoint].position.z)*timeCube);
	
	return newPosition;

}
