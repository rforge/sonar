// FloodFill-Algorithmus zum Finden der Hotspots
// Autor: Michael Bothmann
// Datum: Juli 2013

#ifndef __FLOODFILL_H_
#define __FLOODFILL_H_
#include <vector>
#include <R.h>
#include <Rmath.h>
//#include <R_ext/Lapack.h>
//#include <R_ext/BLAS.h>

class CFloodFill
{
public:
	CFloodFill(double* pImage, int xDimension, int yDimension);
	virtual ~CFloodFill();

	bool	FloodCluster(const int iOldColor);

private:
	void	Fill4(int x, int y, int iOldColor, int iNewColor);
	double	GetPixel(int x, int y);
	bool	MarkPixel(int x, int y, int newColor);

	double*	m_pImage;
	int		m_xDimension;
	int		m_yDimension;

	std::vector<std::pair<int, int> >		m_stack;

	//lint -e1704, diese Elemente sind privat, da sie nicht benoetigt werden
	CFloodFill(const CFloodFill& rSrc);
	CFloodFill &operator=(const CFloodFill& rSrc);
	//lint +e1704
};
#endif //__FLOODFILL_H_