// FloodFill-Algorithmus zum Finden der Hotspots
// Autor: Michael Bothmann
// Datum: Juli 2013

#include "FloodFill.h"
#include <utility>
/*
#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
*/

CFloodFill::CFloodFill(double* pImage , int xDimension, int yDimension)
	:m_pImage(pImage)
	,m_xDimension(xDimension)
	,m_yDimension(yDimension)
	,m_stack()
{
}

CFloodFill::~CFloodFill()
{
}

void CFloodFill::Fill4( int xStart, int yStart, int iOldColor, int iNewColor )
{
	if (m_pImage == 0)
	{
		return;
	}
	m_stack.push_back(std::make_pair(xStart, yStart));
	
	
	while (m_stack.empty() == false)
	{
		std::pair<int, int> point = m_stack.back();
		m_stack.pop_back();

		int x = point.first;
		int y = point.second;
		
		if (GetPixel(x, y) == iOldColor) 
		{
			MarkPixel(x, y, iNewColor);

			m_stack.push_back(std::make_pair(x, y + 1));
			m_stack.push_back(std::make_pair(x, y - 1));
			m_stack.push_back(std::make_pair(x - 1, y));
			m_stack.push_back(std::make_pair(x + 1, y));
		}
	}
	return;
}

double CFloodFill::GetPixel( int x, int y )
{
	if (x < 0 || x >= m_xDimension)
	{
		return -1;
	}
	if (y < 0 || y >= m_yDimension)
	{
		return -1;
	}

	return m_pImage[ (y * m_xDimension) + x];
}

bool CFloodFill::MarkPixel( int x, int y, int newColor )
{
	if (x < 0 || x >= m_xDimension)
	{
		return false;
	}
	if (y < 0 || y >= m_yDimension)
	{
		return false;
	}

	m_pImage[ (y * m_xDimension) + x] = newColor;
	return true;
}

bool CFloodFill::FloodCluster(const int iOldColor)
{
	int iColor = 1;

	for (int y = 0; y < m_yDimension; y++)
	{	
		for (int x = 0; x < m_xDimension; x++)
		{
			if (GetPixel(x,y) == iOldColor)
			{
				iColor++;
				Fill4(x, y, iOldColor, iColor);
			}
		}
	}
	return true;
}
