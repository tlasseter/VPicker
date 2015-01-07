
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System


package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

/** A general scratch class for information about a point relative to a section.
 *  Used in construction of surfacePoints and surface edges (intersection of
 *  surface and sections).
 */
public class StsSectionPoint implements Cloneable
{
    /** @param i row location if on grid                        */
	public int i;
    /** @param j col location if on grid                        */
	public int j;
    /** @param gridRowF float values of grid row                */
	public float gridRowF;
    /** @param gridColF float values of grid col                      */
	public float gridColF;
    /** @param sectionRowF float value of section row           */
	public float sectionRowF;
    /** @param sectionColF float value of section col           */
	public float sectionColF = -1.0f;
    /** @param rowOrCol alignment with sectionGrid ROW, COL, or NONE    */
	public int rowOrCol;
    /** @param point location of this point                     */
	public StsPoint point;
    /** @param nearestPoint location of nearest section point   */
	public StsPoint nearestPoint;
    /** @param section nearest section                          */
	public StsSection section;
    /** @param side side of section we're on                    */
	public int side;
	/** @param sideSet indicates side has been initialized      */
	public boolean sideSet = false;
	/** @param offSection indicates off ends of the section */
	public boolean offSection = false;
    /** @param next queue link: used by faultSurface        */
	public StsSectionPoint next = null;
    /** @param directionToSection direction on row/col towards section */
	public int directionToSection;

	/** @param gridType  If gridType is X, X can vary and point and nearestPoint
	 *  have the same Y and Z; if Y, Y can vary and they have the same X and Z;
	 *  if XY, X and Y can vary and they have the same Z.
	 */
    public int gridType = StsPoint.XY;

    /** @param distance Distance from point to nearestPoint. */
	public float distance;

	/** @param distanceType coordinates used to measure distance between point and
	 *  nearestPoint.  If X, distance is in X-direction; if Y, distance is in Y-direction;
	 *  if XY, distance is measured in X and Y directions.  Later could add XYZ, XZ, etc.
	 */
    public int distanceType = StsPoint.XY;

    /** @param connect direction this point is connected: MINUS or PLUS */
	public int connect = 0;
    static public final int MINUS = StsParameters.MINUS;
    static public final int PLUS = StsParameters.PLUS;

	/** @param noGapPoints number of points in this row/col which are skipped in the
	 *  connect direction before good points are found.  This is the number of points
	 *  which need to be subsequently interpolated.
	 */
    public int noGapPoints;

	/** @param slope slope along the row or column away from the fault used to compute
	 *  z values for the gap points; since gap points are crossed by row and column,
	 *  the weights will be used to compute a weighted value of z
	 */
    public float slope;

	/** angleWeight is the cosine of the angle from the horizontal surface normal
	 *  to the grid line if on a grid line or the line from point to nearestPoint
	 *  if not on a grid line.
	 */
    public float angleWeight = 0.0f;

    /** @param sigma standard-deviation of Z for this point */
	public float sigma = StsParameters.largeFloat;
    /** @param weight used in smoothing edge */
	public float weight = 0.0f;
    /** @param maxGapMul How far we should search as mult times gap size */
	static final float maxGapMult = 4.0f;

    /** Statics used to move additional info in/out of section calculations.
     *  We are assuming that once retrieved, the user doesn't need them associated
     *  with sectionPoint anymore.
     */
    /** @param tangent horizontal tangent vector on section     */
	static public StsPoint tangent;
//    static public boolean wantTangent = false;  @param wantTangent compute tangent vector and side


    // State flags for rowOrCol
    static final int ROW = StsParameters.ROW;
    static final int COL = StsParameters.COL;
    static final int ROWCOL = StsParameters.ROWCOL;
    static final int NONE= StsParameters.NONE;

    static final float nullValue = StsParameters.nullValue;
    static final float roundOff = StsParameters.roundOff;

    public StsSectionPoint()
    {
        point = new StsPoint(3);
    }

    public StsSectionPoint(StsPoint point)
    {
        if (point.getLength() == 3)
            this.point = point;
        else
            this.point = new StsPoint(point.getXYZorT());

    }

    public StsSectionPoint(StsPoint point, float colF)
    {
        this.point = point;
//        this.point = new StsPoint(point);
        sectionColF = colF;
    }
/*
    public StsSectionPoint(StsPoint point, float sectionF, StsSection section)
    {
        this.point = new StsPoint(point);
        this.sectionF = sectionF;
        checkSetRowCol(section);
    }

    public StsSectionPoint(int i, int j, StsSurface surface)
    {
        this.i = i;
        this.j = j;
        this.gridRowF = (float)i;
        this.gridColF = (float)j;
        this.point = new StsPoint(surface.getGrid().getPointXYZ(i, j));

    }
*/
    public StsSectionPoint(int i, int j, StsSection section, StsModelSurface surface)
    {
        this.i = i;
        this.j = j;
        this.gridRowF = (float)i;
        this.gridColF = (float)j;
        this.section = section;
        this.point = new StsPoint(surface.getXYZorT(i, j));
    }

    public StsSectionPoint(int i, int j, StsSection section, StsModelSurface surface, int rowOrCol)
    {
        this.i = i;
        this.j = j;
        this.gridRowF = (float)i;
        this.gridColF = (float)j;
        this.section = section;
        this.point = new StsPoint(surface.getXYZorT(i, j));
        float z = point.getZorT();
        point.setZorT(StsMath.checkRoundOffInteger(z));
        this.rowOrCol = rowOrCol;
    }
/*
    public StsSectionPoint(float gridRowF, float gridColF, float z, float sectionF, int gridType,
                           int connect, StsSurface surface, StsSection section)
    {
        this.gridRowF = gridRowF;
        this.gridColF = gridColF;
        this.sectionF = sectionF;
        this.gridType = gridType;
        this.connect = connect;
        this.section = section;
        point = surface.getGrid().getXYPoint(gridRowF, gridColF);
    }

    public StsSectionPoint(StsPoint point, float sectionF, float gridRowF, float gridColF,
                           int connect, int gridType, StsSection section)
    {
        this.point = point;
        this.sectionF = sectionF;
        this.gridRowF = gridRowF;
        this.gridColF = gridColF;
        this.i = (int)gridRowF;
        this.j = (int)gridColF;
        this.connect = connect;
        this.gridType = gridType;
        this.section = section;
    }

    public StsSectionPoint(StsPoint point, float sectionF, StsSection section, StsGridable grid)
    {
        this.point = point;
        this.sectionF = sectionF;
        this.section = section;
        gridRowF = grid.getRowCoor(point.getY());
        gridColF = grid.getColCoor(point.getX());
        i = (int)gridRowF;
        j = (int)gridColF;
    }

    public StsSectionPoint(StsSectionPoint sectionPoint0, StsSectionPoint sectionPoint1, float f)
    {
        point = new StsPoint();
        point.interpolatePoints(sectionPoint0.point, sectionPoint1.point, f);
        section = sectionPoint0.section;
        gridRowF = sectionPoint0.gridRowF + f*(sectionPoint1.gridRowF - sectionPoint0.gridRowF);
        gridColF = sectionPoint0.gridColF + f*(sectionPoint1.gridColF - sectionPoint0.gridColF);
        i = (int)gridRowF;
        j = (int)gridColF;
    }

    public StsSectionPoint(StsGridSectionPoint gridPoint)
    {
        point = gridPoint.getPoint();
        section = gridPoint.getSection();
        gridRowF = gridPoint.getGridRowF();
        gridColF = gridPoint.getGridColF();
        i = (int)gridRowF;
        j = (int)gridColF;
        sectionF = gridPoint.getSectionF();
    }
*/
    public float getGridRowF() { return gridRowF; }
    public float getGridColF() { return gridColF; }
    public float getSectionRowF() { return sectionRowF; }
    public float getSectionColF() { return sectionColF; }

    public int getRowOrCol() { return rowOrCol; }
    public int getConnect() { return connect; }

    public void checkSetRowCol(StsSection section)
    {
        sectionRowF = section.getRowF(point.getZ());

        if(nearestPoint == null)
			sectionColF = section.getColF(point);

        boolean isRow = isIntegralSectionRowF();
        boolean isCol = isIntegralSectionColF();

        if(isRow)
        {
            if(isCol) rowOrCol = ROWCOL;
            else
                rowOrCol = ROW;
        }
        else if(isCol)
            rowOrCol = COL;
        else
            rowOrCol = NONE;
    }

    public boolean isIntegralSectionRowF()
    {
        float value = sectionRowF;
        int ivalue;
        float remainder;

        if(value == 0.0f) return true;

        ivalue = (int)value;
        remainder = Math.abs(value - (float)ivalue);
        if(remainder <= roundOff )
        {
            sectionRowF = (float)ivalue;
            return true;
        }
        else if(1.0f - remainder <= roundOff)
        {
            sectionRowF = (float)(ivalue+1);
            return true;
        }
        else
            return false;
    }

    public boolean isIntegralSectionColF()
    {
        float value = sectionColF;
        int ivalue;
        float remainder;

        if(value == 0.0f) return true;

        ivalue = (int)value;
        remainder = Math.abs(value - (float)ivalue);
        if(remainder <= roundOff )
        {
            sectionColF = (float)ivalue;
            return true;
        }
        else if(1.0f - remainder <= roundOff)
        {
            sectionColF = (float)(ivalue+1);
            return true;
        }
        else
            return false;
    }

    public void setSide(int side)
    {
        this.side = side;
        sideSet = true;
    }

    public boolean insideGap(float gridInc)
    {
        if(section == null) return false;
        return distance < section.getFaultGap(side)*gridInc;
    }

    /** Given two sectionPoints in a row or column interpolate a point
     *  in between which is estimated to be on the fault section using
     *  the distances of the two points to the fault section.  This is
     *  only a starting point and will be later adjusted for fault throw.
     */
/*
    public void interpolateSectionPointsToSection(StsSectionPoint point0,
                                         StsSectionPoint point1, int gridType)
    {
        float d0, d1, f, colf0, colf1;

        if(point0.section != point1.section)
        {
            StsException.outputException(new StsException("StsSectionPoint.interpolateSectionPoints(...)",
                                         "sections don't match for two section points."));
            return;
        }
        else
            this.section = point0.section;

        d0 = point0.distance;
        d1 = point1.distance;

        f = d0/(d0 + d1);

        point.interpolatePoints(point0.point, point1.point, f);

        colf0 = point0.sectionColF;
        colf1 = point1.sectionColF;
        sectionColF = colf0 + f*(colf1 - colf0);

        this.gridType = gridType;

        if(gridType == StsPoint.X)      // points have the same Y (row) values: we can adjust X
        {
            gridColF = point0.gridColF + f*(point1.gridColF - point0.gridColF);
            this.j = (int)gridColF;
            this.i = point0.i;
            gridRowF = (float)i;
        }
        else if(gridType == StsPoint.Y) // points have the same X (col) values: we can adjust Y
        {
            gridRowF = point0.gridRowF + f*(point1.gridRowF - point0.gridRowF);
            this.i = (int)gridRowF;
            this.j = point0.j;
            gridColF = (float)j;
        }
        else                          // points not constant in X or Y: we can adjust X & Y
            gridType = StsPoint.XY;
    }


    public void interpolatePoints(StsPoint point0, StsPoint point1, float f)
    {
        point.interpolatePoints(point0, point1, f);
    }

    public void interpolateSectionPoints(StsSectionPoint point0, StsSectionPoint point1, StsSection section, float f)
    {
        point.interpolatePoints(point0.point, point1.point, f);
        sectionColF = point0.sectionColF + f*(point1.sectionColF - point0.sectionColF);
        checkSetRowCol(section);
    }
*/
    static public StsSectionPoint getSectionPointOnSection(int row, int col, float sectionColF, StsSection section, StsModelSurface surface)
    {
        StsSectionPoint sectionPoint = new StsSectionPoint(row, col, section, surface);
        sectionPoint.sectionColF = sectionColF;
        section.computeNearestPoint(sectionPoint);

        return sectionPoint;
    }

    /** Called from StsSection when a point has been adjusted on the section
     *  in order to save all the useful things we need to know.
     */
    public void computeSectionValues(StsSection section, int index, float f, StsPoint[] linePoints)
    {
        if(linePoints == null) return;

        try
        {
            int nLinePoints = linePoints.length;

            /**Limit the nearest point to the extent of the line, however */

            if(nearestPoint == null) nearestPoint = new StsPoint(3);

            /** Save sectionF as the actual position relative to the line;
             *  which may extend off the line (<0.0 or > 1.0).
             *  Set the flag offSection = true in this case.
             */
            sectionColF = index + f;

            if(sectionColF < 0.0f)
			{
                nearestPoint.setValues(linePoints[0]);
				offSection = true;
				sectionColF = 0.0f;
			}
            else if(sectionColF > nLinePoints-1)
			{
                nearestPoint.setValues(linePoints[nLinePoints-1]);
				offSection = true;
				sectionColF = (float)(nLinePoints-1);
			}
            else
                nearestPoint.interpolatePoints(linePoints[index], linePoints[index+1], f);

            /** sideSet is initially false.  After computing, it is set to true
             *  so it won't be changed in subsequent operations.
             *  If you want the tangent vector computed, set wantVectors = true.
             */

            tangent = StsPoint.subPointsStatic(linePoints[index+1], linePoints[index]);
            tangent.normalize();

            if(!sideSet)
            {
                StsPoint normal = StsPoint.subPointsStatic(point, nearestPoint);
                distance = normal.normalizeReturnLength();

                /** AngleWeight is probably not used by points which are not on
                 *  a gridLine, but the calculation is included here for completeness
                 */
                float sideValue = StsPoint.crossProduct2D(normal, tangent);
                if(sideValue > 0.0f)
                {
                    side = StsSection.RIGHT;
                    angleWeight = sideValue;
                }
                else
                {
                    side = StsSection.LEFT;
                    angleWeight = -sideValue;
                }

                sideSet = true;
            }

            /** angleWeight is the cosine of the angle from the surface normal
             *  to the grid line.
             */
            if(gridType == StsPoint.X)
                angleWeight = Math.abs(tangent.v[1]);
            else if(gridType == StsPoint.Y)
                angleWeight = Math.abs(tangent.v[0]);

            /** set sectionRowF, sectionColF, and rowOrCol according to
             *  nearestPoint.z
             */
             checkSetRowCol(section);
        }
        catch(Exception e)
        {
            StsException.outputException("StsSectionPoint.computeSectionValues() failed.",
                e, StsException.WARNING);
        }
    }

    /** Construct a scratch sectionPoint and get horizontal tangent
     *  on section at this point and return it
     */
    public static StsPoint getTangentAtPoint(StsPoint point, StsSection section)
    {
        StsSectionPoint sectionPoint = new StsSectionPoint(point);
        section.computeNearestPoint(sectionPoint);
        return sectionPoint.tangent;
    }
/*
    public float computeNearestVerticalGridPoint(StsGridable grid)
    {
        nearestPoint = grid.getNearestXYZ(point.v, false);
        return nearestPoint.getZ();
    }
*/
    public float getDZ()
    {
        return point.getZ() - nearestPoint.getZ();
    }
/*
    public static StsList getGridCrossings(StsGridable grid, StsSection section,
                          StsSectionPoint sectionPoint0, StsSectionPoint sectionPoint1)
    {
        return getGridCrossings(grid, section, sectionPoint0, sectionPoint1, StsParameters.NONE);
    }
*/
    /** Creates an StsList of StsGridSectionPoints on gridCrossings and then
     *  converts them to StsSectionPoints for convenience.
     */
/*
    public static StsList getGridCrossings(StsGridable grid, StsSection section,
                          StsSectionPoint sectionPoint0, StsSectionPoint sectionPoint1,
                          int extrapolate)
    {
        StsList gridPoints = StsGridSectionPoint.getGridCrossings(grid, section,
                                sectionPoint0, sectionPoint1, extrapolate, false);


        int nPoints = gridPoints.getSize();

        StsList sectionPoints = new StsList(nPoints);
        for(int n = 0; n < nPoints; n++)
        {
            StsGridSectionPoint gridPoint = (StsGridSectionPoint)gridPoints.getElement(n);
            StsSectionPoint sectionPoint = new StsSectionPoint(gridPoint);
            sectionPoint.computeNearestVerticalGridPoint(grid);
            sectionPoints.add(sectionPoint);
        }
        return sectionPoints;
    }

    public void addToPointLinkedGrid(StsPointLinkedGrid pointLinkedGrid)
    {
        try
        {
            if(rowOrCol == NONE) return;
            StsLink link = new StsLink(this, rowOrCol, sectionRowF, sectionColF, 0, 0);
            pointLinkedGrid.insertLink((StsLink)link);
        }
        catch(Exception e)
        {
            StsException.outputException("StsSectionPoint.addToPointLinkedGrid() failed",
                e, StsException.WARNING);
        }
    }

    public float compare(StsGridPointLinkable object)
    {
        float compareValue, otherCompareValue;

        if(!(object instanceof StsSectionPoint))
        {
            StsException.outputException(new StsException(StsException.WARNING,
                "StsSectionPoint.compare() failed.",
                "Object: " + object.getClass().toString() + " not an StsSectionPoint."));
            return StsParameters.largeFloat;
        }
        return getCompareValue() - object.getCompareValue();
    }

    public float getCompareValue()
    {
        if(rowOrCol == ROW)
            return sectionRowF;
        else if(rowOrCol == COL)
            return sectionColF;
        else
            return nullValue;
    }
*/
    public void setDirectionToSection(int insideDirection, int correctSide)
    {
        if(side == correctSide) directionToSection = -insideDirection;
        else directionToSection = insideDirection;
    }

    public boolean isDirection(int rowOrCol, int directionToSection)
    {
        return rowOrCol == this.rowOrCol && this.directionToSection == directionToSection;
    }

    public String toString()
    {
        return new String(" StsSectionPoint - gridRowF: " + gridRowF + " gridColF: " + gridColF +
            " sectionColF: " + sectionColF);
    }
}
