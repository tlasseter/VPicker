
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

public class StsLineZone implements StsEdgeLinkable
{
    // instance fields
    protected StsLine line;
    protected StsSurfaceVertex topVertex = null;
    protected StsSurfaceVertex botVertex = null;

    protected StsList edgePoints;

    static final int NONE = StsParameters.NONE;

    /** constructors */

	public StsLineZone()
	{
	}

    public StsLineZone(StsModel model, StsSurfaceVertex topVertex, StsSurfaceVertex botVertex, StsLine line)
    {
        this.topVertex = topVertex;
        this.botVertex = botVertex;
        this.line = line;
        constructEdgePoints(model);
    }

    public StsLineZone(StsLine line, int type, String name,
                       StsSurfaceVertex topVertex, StsSurfaceVertex botVertex)
    {
        if (name==null)
        {
            StsException.systemError("StsLineZone.StsLineZone:" +
                                     " Cannot create a line zone with a null name.");
        }
        if (line==null)
        {
            StsException.systemError("StsLineZone.StsLineZone:" +
                                     " Cannot create a zone for a null line.");
        }

	    this.line = line;
	    setTop(topVertex);
	    setBot(botVertex);
    }

    // included for interface compatibility
    public int getRowOrCol() { return -1; }
    public int getRowCol() { return -1; }
    public StsEdgeLinkable getNextEdge() { return null; }
    public StsEdgeLinkable getPrevEdge() { return null; }
    public StsSurfaceVertex getNextVertex() { return botVertex; }
    public StsSurfaceVertex getPrevVertex() { return topVertex; }

    // Accessors
    public StsLine getLine(){ return line; }
    public void setTop(StsSurfaceVertex topVertex)
    {
        if (topVertex == null)
        {
            StsException.systemError("StsLineZone.setTop:"
                      + " Cannot create a zone with a null vertex");
        }
        this.topVertex = topVertex;
        if (botVertex==null) return;

        float topZ = getTopZ();
        float baseZ = getBotZ();
    }
    public StsSurfaceVertex getTop(){ return topVertex; }
    public float getTopZ()
    {
        if (topVertex==null) StsException.systemError("StsLineZone.getTopZ:"
                          + " Null top vertex.");
        return topVertex.getPoint().getZorT();
    }
    public void setBot(StsSurfaceVertex botVertex)
    {
        if (botVertex == null)
        {
            StsException.systemError("StsLineZone.setBot:"
                      + " Cannot create a zone with a null vertex");
        }
        this.botVertex = botVertex;
        if (topVertex==null) return;

        float topZ = getTopZ();
        float baseZ = getBotZ();
    }
    public StsSurfaceVertex getBot(){ return botVertex; }
    public float getBotZ()
    {
        if (botVertex==null) StsException.systemError("StsLineZone.getBotZ:"
                          + " Null top vertex.");
        return botVertex.getPoint().getZorT();
    }

    private void constructEdgePoints(StsModel model)
    {
        StsGridSectionPoint topSectionPoint, botSectionPoint;
        StsObjectList sectionPoints;
        StsGridSectionPoint sectionPoint;

        try
        {
            boolean persistent = false;

            topSectionPoint = topVertex.getSurfacePoint();
            botSectionPoint = botVertex.getSurfacePoint();

            float zMin = topSectionPoint.getXYZorT()[2];
            float zMax = botSectionPoint.getXYZorT()[2];

            StsProject project = model.getProject();

            int rowMin = project.getIndexBelow(zMin);
            int rowMax = project.getIndexAbove(zMax);

            int nRows = 0;
            if(rowMax >= rowMin) nRows = rowMax - rowMin + 1;

            sectionPoints = new StsObjectList(nRows + 2, 10);

            // add top and bot points and any section row crossings
            sectionPoints.add(topSectionPoint);

			if(!line.areAllSectionsVertical())
			{
				for(int row = rowMin; row <= rowMax; row++)
				{
					float z = project.getZAtIndex(row);
					StsPoint point = line.getXYZPointAtZorT(z, true);
					sectionPoint = new StsGridSectionPoint(point, line, false);
					sectionPoint.setSectionRow();
					sectionPoints.add(sectionPoint);
				}
			}
            sectionPoints.add(botSectionPoint);

            // add section col crossings for section this line is on and sectionRowCols for connected sections
            sectionPoints = addSectionColCrossings(sectionPoints,  persistent);

			// Add points on grid crossings
			sectionPoints = addGridCrossings(sectionPoints, persistent);

            edgePoints = sectionPoints;
        }
        catch(Exception e)
        {
            StsException.outputException("StsLineZone.constructEdgePoints() failed.",
                e, StsException.WARNING);
        }
    }

    private StsObjectList addSectionColCrossings(StsObjectList sectionPoints, boolean persistent)
    {
        StsGridSectionPoint prevSectionPoint, sectionPoint;

        try
        {
            StsSection onSection = line.getOnSection();
            if(onSection != null)
			{
				int nSectionPoints = sectionPoints.getSize();
				StsObjectList newSectionPoints = new StsObjectList(nSectionPoints, 10);

				sectionPoint = (StsGridSectionPoint)sectionPoints.getElement(0);
				newSectionPoints.add(sectionPoint);

				for(int n = 1; n < nSectionPoints; n++)
				{
					prevSectionPoint = sectionPoint;
					sectionPoint = (StsGridSectionPoint)sectionPoints.getElement(n);
					StsList gridPoints = getColCrossings(prevSectionPoint, sectionPoint, onSection, persistent);
					newSectionPoints.addList(gridPoints);
					newSectionPoints.add(sectionPoint);
				}
				sectionPoints = newSectionPoints;
			}

			StsSection[] connectedSections = line.getConnectedSections();
			for(int n = 0; n < connectedSections.length; n++)
			{
				StsSection section = connectedSections[n];
				int col = section.getLineIndex(line);
				if(col < 0) continue;
				int nSectionPoints = sectionPoints.getSize();
				for(int i = 0; i < nSectionPoints; i++)
				{
					sectionPoint = (StsGridSectionPoint)sectionPoints.getElement(i);
					sectionPoint.setRowOrColIndex(section, StsParameters.COL, col, false);
				}
			}
            return sectionPoints;
        }
        catch(Exception e)
        {
            StsException.outputException("StsLineZone.addSectionRowColCrossings() failed.",
                e, StsException.WARNING);
            return sectionPoints;
        }
    }

    private StsList getColCrossings(StsGridSectionPoint point0, StsGridSectionPoint point1,
                                    StsSection section, boolean persistent)
    {
        float colF0, colF1;
        int col0, col1, inc, nCols, col, n;
        StsGridSectionPoint point;

        colF0 = point0.getSectionRowCol(section).getColF();
        colF1 = point1.getSectionRowCol(section).getColF();

        if(colF1 > colF0)
        {
            col0 = StsMath.above(colF0);
            col1 = StsMath.below(colF1);
            inc = 1;
            nCols = col1 - col0 + 1;
        }
        else
        {
            col0 = StsMath.below(colF0);
            col1 = StsMath.above(colF1);
            inc = -1;
            nCols = col0 - col1 + 1;
        }

        if(nCols <= 0) return null;

        StsList colPoints = new StsList(nCols);
        for(n = 0, col = col0; n < nCols; n++, col += inc)
        {
            double f = (double)((col - colF0)/(colF1 - colF0));
            point = StsGridSectionPoint.sectionInterpolate(point0, point1, f, this, section, persistent);
			point.setRowOrColIndex(section, StsParameters.COL, col);
//            point.getSectionRowCol(section).setColF((float)col);
            colPoints.add(point);
        }
        return colPoints;
    }


    private StsObjectList addGridCrossings(StsObjectList sectionPoints, boolean persistent)
    {
        StsGridSectionPoint sectionPoint, prevSectionPoint;
        StsList gridPoints;

        try
        {
            if(sectionPoints == null) return null;
			if(!line.hasConnectedSection()) return sectionPoints;

            int nSectionPoints = sectionPoints.getSize();
            StsObjectList newSectionPoints = new StsObjectList(nSectionPoints, 10);

            sectionPoint = (StsGridSectionPoint)sectionPoints.getElement(0);
            newSectionPoints.add(sectionPoint);

			StsSection onSection = line.getOnSection();
			StsSection[] connectedSections = line.getAllSections();
			int nConnectedSections = connectedSections.length;

            for(int n = 1; n < nSectionPoints; n++)
            {
                prevSectionPoint = sectionPoint;
                sectionPoint = (StsGridSectionPoint)sectionPoints.getElement(n);

                gridPoints = StsGridSectionPoint.getGridCrossings(prevSectionPoint, sectionPoint, onSection, persistent, NONE);

				if(gridPoints.getSize() > 0)
				{
                    /*
                    for(int s = 0; s < nConnectedSections; s++)
					{
						StsSection connectedSection = connectedSections[s];
						addConnectedSectionRowCols(gridPoints, connectedSection);
					}
					*/
					newSectionPoints.addList(gridPoints);
				}
                newSectionPoints.add(sectionPoint);
            }
            return newSectionPoints;
        }
        catch(Exception e)
        {
            StsException.outputException("StsLineZone.addGridCrossings() failed.",
                e, StsException.WARNING);
            return sectionPoints;
        }
    }

    private void addConnectedSectionRowCols(StsList gridPoints, StsSection section)
    {
        StsGridSectionPoint gridPoint;

        if(section == null || gridPoints == null) return;

        try
        {
            int nGridPoints = gridPoints.getSize();
            for(int n = 0; n < nGridPoints; n++)
            {
                gridPoint = (StsGridSectionPoint)gridPoints.getElement(n);
                gridPoint.addSectionRowCol(section, line);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsLineZone.addSectionRowCols() failed.",
                e, StsException.WARNING);
        }
    }

    public StsList getGridEdgePointsList()
    {
        return getEdgePointsList();
    }

    public StsList getEdgePointsList()
    {
        return edgePoints;
    }

    public StsGridSectionPoint getFirstGridPoint()
    {
        if(edgePoints == null) return null;
        else return (StsGridSectionPoint)edgePoints.getFirst();
    }

    public StsGridSectionPoint getLastGridPoint()
    {
        if(edgePoints == null) return null;
        else return (StsGridSectionPoint)edgePoints.getLast();
    }

    public float[][] getXYZPoints()
    {
        StsList points = getEdgePointsList();
        if(points == null) return null;
        int nPoints = points.getSize();
        float[][] xyzPoints = new float[nPoints][];
        for(int n = 0; n < nPoints; n++)
        {
            StsGridSectionPoint gridPoint = (StsGridSectionPoint)points.getElement(n);
            xyzPoints[n] = gridPoint.getPoint().getPointValues();
        }
        return xyzPoints;
    }

    public boolean delete()
    {
        // topVertex.deleteEdgeFromVertex(this);
       //  botVertex.deleteEdgeFromVertex(this);
        if(edgePoints != null) edgePoints.deleteAll();
        return true;
    }

    public String getLabel()
    {
        return new String("LineZone-" + line.getName() + "-" + topVertex.getName());
    }
}














