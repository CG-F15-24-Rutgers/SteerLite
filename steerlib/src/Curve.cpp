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
#include <unordered_set>

using namespace Util;

static bool bRet = false;

bool IsEqualTime(Util::CurvePoint a, Util::CurvePoint b)
{
	return (a.time == b.time);
}
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
	Util::Point outPoint, previousPoint(controlPoints[0].position);
	bRet = true;
	int index = 0;
	for (index = 0; index < controlPoints.size() - 1; ++index)
	{
		for (float i = 0.0; i < 1.0f; i += 0.001)
		{
			if (type == Util::hermiteCurve)
				outPoint = useHermiteCurve(index, i);
			else
				outPoint = useCatmullCurve(index, i);
			DrawLib::drawLine(previousPoint, outPoint, curveColor, curveThickness);
			previousPoint = outPoint;
		}
	}

	return;
#endif
}

// Sort controlPoints vector in ascending order: min-first
void Curve::sortControlPoints()
{
	/*std::unordered_set<CurvePoint> s;

	for (int i = 0; i < controlPoints.size(); ++i)
	s.insert(controlPoints[i]);
	controlPoints.assign(s.begin(), s.end());
	std::sort(controlPoints.begin() + 1, controlPoints.end(), &compareTime);
	*/
	std::sort(controlPoints.begin() + 1, controlPoints.end(), &compareTime);
	controlPoints.erase(std::unique(controlPoints.begin(), controlPoints.end(), &IsEqualTime), controlPoints.end());
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
	int index =  0;
	for (int i = 0; i < controlPoints.size() - 1; )
	{
		if (time > controlPoints[i + 1].time) {
			++i;
			++index;
			continue;
		}
		nextPoint = index;
		return true;
	}
	return false;
}

// Implement Hermite curve
Point Curve::useHermiteCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	///float normalTime, intervalTime;
	float nTime;
	if (nextPoint >= controlPoints.size() - 1)
		return newPosition;
	float timeDuration = controlPoints[nextPoint + 1].time - controlPoints[nextPoint].time;
	if (!bRet)
	{
		nTime = (time - controlPoints[nextPoint].time) / timeDuration;
	}
	else
	{
		nTime = time;
	}

	float timeSquare = nTime*nTime;
	float timeCube = timeSquare*nTime;

	newPosition.x = (((2 * timeCube) - (3 * timeSquare) + 1)*controlPoints[nextPoint].position.x) +
		(((3 * timeSquare) - (2 * timeCube))*controlPoints[nextPoint + 1].position.x) +
		((timeCube - (2 * timeSquare) + nTime)*controlPoints[nextPoint].tangent.x)*timeDuration +
		((timeCube - timeSquare)*controlPoints[nextPoint + 1].tangent.x)*timeDuration;
	newPosition.y = (((2 * timeCube) - (3 * timeSquare) + 1)*controlPoints[nextPoint].position.y) +
		(((3 * timeSquare) - (2 * timeCube))*controlPoints[nextPoint + 1].position.y) +
		((timeCube - (2 * timeSquare) + nTime)*controlPoints[nextPoint].tangent.y)*timeDuration +
		((timeCube - timeSquare)*controlPoints[nextPoint + 1].tangent.y)*timeDuration;
	newPosition.z = (((2 * timeCube) - (3 * timeSquare) + 1)*controlPoints[nextPoint].position.z) +
		(((3 * timeSquare) - (2 * timeCube))*controlPoints[nextPoint + 1].position.z) +
		((timeCube - (2 * timeSquare) + nTime)*controlPoints[nextPoint].tangent.z)*timeDuration +
		((timeCube - timeSquare)*controlPoints[nextPoint + 1].tangent.z)*timeDuration;
	//std::cout << newPosition.x << " " << newPosition.y << " " << newPosition.z << std::endl;
	return newPosition;
}

// Implement Catmull-Rom curve
Point Curve::useCatmullCurve(const unsigned int nextPoint, const float time)
{
	Point newPosition;
	//float normalTime, intervalTime;
	float nTime;
	if (nextPoint >= controlPoints.size() - 1)
		return newPosition;
	//std::cout << nextPoint << std::endl;
	float timeDuration = controlPoints[nextPoint + 1].time - controlPoints[nextPoint].time;
	if (!bRet)
	{
		nTime = (time - controlPoints[nextPoint].time) / timeDuration;
	}
	else
	{
		nTime = time;
	}
	static float s0x, s0y, s0z, s1x, s1y, s1z, s2x, s2y, s2z, s3x, s3y, s3z;
	float timeSquare = nTime*nTime;
	float timeCube = timeSquare*nTime;

	float t0, t1, t2, t3;
	Point *prevTan, *nextTan;

	if (nextPoint < controlPoints.size() - 3)
	{
		t0 = controlPoints[nextPoint].time;
		t1 = controlPoints[nextPoint + 1].time;
		t2 = controlPoints[nextPoint + 2].time;
		t3 = controlPoints[nextPoint + 3].time;
	}

	if (controlPoints.size() - nextPoint > 3)
	{
		s0x = ((t2 - t0) / (t2 - t1))* ((controlPoints[nextPoint + 1].position.x - controlPoints[nextPoint].position.x) / (t1 - t0)) - ((controlPoints[nextPoint + 2].position.x - controlPoints[nextPoint].position.x) / (t2 - t0))* ((t2 - t0) / (t2 - t1));
		s0y = ((t2 - t0) / (t2 - t1))* ((controlPoints[nextPoint + 1].position.y - controlPoints[nextPoint].position.y) / (t1 - t0)) - ((controlPoints[nextPoint + 2].position.y - controlPoints[nextPoint].position.y) / (t2 - t0))* ((t2 - t0) / (t2 - t1));
		s0z = ((t2 - t0) / (t2 - t1))* ((controlPoints[nextPoint + 1].position.z - controlPoints[nextPoint].position.z) / (t1 - t0)) - ((controlPoints[nextPoint + 2].position.z - controlPoints[nextPoint].position.z) / (t2 - t0))* ((t2 - t0) / (t2 - t1));

		s1x = ((t1 - t0) / (t2 - t0))* ((controlPoints[nextPoint + 2].position.x - controlPoints[nextPoint + 1].position.x) / (t2 - t1)) + ((controlPoints[nextPoint + 1].position.x - controlPoints[nextPoint].position.x) / (t1 - t0))* ((t2 - t1) / (t2 - t0));
		s1y = ((t1 - t0) / (t2 - t0))* ((controlPoints[nextPoint + 2].position.y - controlPoints[nextPoint + 1].position.y) / (t2 - t1)) + ((controlPoints[nextPoint + 1].position.y - controlPoints[nextPoint].position.y) / (t1 - t0))* ((t2 - t1) / (t2 - t0));
		s1z = ((t1 - t0) / (t2 - t0))* ((controlPoints[nextPoint + 2].position.z - controlPoints[nextPoint + 1].position.z) / (t2 - t1)) + ((controlPoints[nextPoint + 1].position.z - controlPoints[nextPoint].position.z) / (t1 - t0))* ((t2 - t1) / (t2 - t0));

		s2x = ((t2 - t1) / (t3 - t1))* ((controlPoints[nextPoint + 3].position.x - controlPoints[nextPoint + 2].position.x) / (t3 - t2)) + ((controlPoints[nextPoint + 2].position.x - controlPoints[nextPoint + 1].position.x) / (t2 - t1))* ((t3 - t2) / (t3 - t1));
		s2y = ((t2 - t1) / (t3 - t1))* ((controlPoints[nextPoint + 3].position.y - controlPoints[nextPoint + 2].position.y) / (t3 - t2)) + ((controlPoints[nextPoint + 2].position.y - controlPoints[nextPoint + 1].position.y) / (t2 - t1))* ((t3 - t2) / (t3 - t1));
		s2z = ((t2 - t1) / (t3 - t1))* ((controlPoints[nextPoint + 3].position.z - controlPoints[nextPoint + 2].position.z) / (t3 - t2)) + ((controlPoints[nextPoint + 2].position.z - controlPoints[nextPoint + 1].position.z) / (t2 - t1))* ((t3 - t2) / (t3 - t1));

		s3x = ((t3 - t1) / (t3 - t2))* ((controlPoints[nextPoint + 2].position.x - controlPoints[nextPoint + 1].position.x) / (t2 - t1)) - ((controlPoints[nextPoint + 3].position.x - controlPoints[nextPoint + 1].position.x) / (t3 - t1))* ((t2 - t1) / (t3 - t2));
		s3y = ((t3 - t1) / (t3 - t2))* ((controlPoints[nextPoint + 2].position.y - controlPoints[nextPoint + 1].position.y) / (t2 - t1)) - ((controlPoints[nextPoint + 3].position.y - controlPoints[nextPoint + 1].position.y) / (t3 - t1))* ((t2 - t1) / (t3 - t2));
		s3z = ((t3 - t1) / (t3 - t2))* ((controlPoints[nextPoint + 2].position.z - controlPoints[nextPoint + 1].position.z) / (t2 - t1)) - ((controlPoints[nextPoint + 3].position.z - controlPoints[nextPoint + 1].position.z) / (t3 - t1))* ((t2 - t1) / (t3 - t2));
	}
	if (nextPoint % 3 == 0)
	{
		prevTan = new Point(s0x, s0y, s0z);
		nextTan = new Point(s1x, s1y, s1z);
	}
	else if (nextPoint % 2 == 1)
	{
		prevTan = new Point(s1x, s1y, s1z);
		nextTan = new Point(s2x, s2y, s2z);
	}
	else
	{
		prevTan = new Point(s2x, s2y, s2z);
		nextTan = new Point(s3x, s3y, s3z);
	}

	newPosition.x = (((2 * timeCube) - (3 * timeSquare) + 1)*controlPoints[nextPoint].position.x) +
		(((3 * timeSquare) - (2 * timeCube))*controlPoints[nextPoint + 1].position.x) +
		((timeCube - (2 * timeSquare) + nTime)*(prevTan->x))*timeDuration +
		((timeCube - timeSquare)*nextTan->x)*timeDuration;
	newPosition.y = (((2 * timeCube) - (3 * timeSquare) + 1)*controlPoints[nextPoint].position.y) +
		(((3 * timeSquare) - (2 * timeCube))*controlPoints[nextPoint + 1].position.y) +
		((timeCube - (2 * timeSquare) + nTime)*prevTan->y)*timeDuration +
		((timeCube - timeSquare)*nextTan->y)*timeDuration;
	newPosition.z = (((2 * timeCube) - (3 * timeSquare) + 1)*controlPoints[nextPoint].position.z) +
		(((3 * timeSquare) - (2 * timeCube))*controlPoints[nextPoint + 1].position.z) +
		((timeCube - (2 * timeSquare) + nTime)*prevTan->z)*timeDuration +
		((timeCube - timeSquare)*nextTan->z)*timeDuration;

	delete prevTan;
	delete nextTan;
	return newPosition;

}
