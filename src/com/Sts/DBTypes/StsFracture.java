package com.Sts.DBTypes;

import com.Sts.Interfaces.StsFractureDisplayable;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Nov 9, 2007
 * Time: 8:17:37 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsFracture extends StsSection implements StsFractureDisplayable
{
    StsFractureSet fractureSet;
    protected boolean drawStimulated = false;  
    protected int areaScale = 1;
    protected float currentZ = StsParameters.nullValue;

    public StsFracture()
    {
    }

    public StsFracture(boolean persistent)
    {
        super(persistent);
    }

    public StsFracture(StsFractureSet fractureSet)
	{
		super(StsSection.FRACTURE, false);
        this.fractureSet = fractureSet;
        setStsColor(fractureSet.getStsColor());
        drawSurface = true;
    }
	
    public void setAreaScale(int scale) 
    { 
    	areaScale = scale; 
    	dbFieldChanged("areaScale", areaScale);
    }
    public int getAreaScale() { return areaScale; }	
    
	public void addToSet()
	{
		fractureSet.addFracture(this);
	}
	
	public boolean deleteSectionEdge(StsSectionEdge edge)
	{
		for(int i=0; i<sectionEdges.getSize(); i++)
		{
			if(edge.equals(sectionEdges.getElement(i)))
			{
				edge.delete();
				sectionEdges.delete(i);
				return true;
			}
		}
		return false;
	}
	
    public void setDrawStimulated(boolean draw)
    {
    	drawStimulated = draw;
    }
    public boolean getDrawStimulated() { return drawStimulated; }
    
	public void display(StsGLPanel3d glPanel3d, boolean displaySectionEdges)
	{
		super.display(glPanel3d, displaySectionEdges);
		if(drawStimulated) drawStimulatedArea();
	}
	
	private void drawStimulatedArea()
	{
		int scale = 1;
	    float zCoor = currentModel.win3d.getCursor3d().getCurrentDirCoordinate(StsCursor3d.ZDIR);

	    StsPoint[] points = getSectionIntersectionAtZ(zCoor);
	    if(points == null)	return;
	    if(points.length <= 0) return;
	    for(int i=0; i<points.length; i++)
	    {
	    	if(points[i] == null) continue;
	    	scale = getAreaScale();  
    		StsGLDraw.drawEllipse(points[i].getXYZ(), StsColor.BLUE, currentModel.getGlPanel3d(), scale, scale, 0.0f, 2.0f);
    		StsGLDraw.drawFilledEllipse(points[i].getXYZ(), StsColor.BLUE, currentModel.win3d.getGlPanel3d(), scale, scale, 0.0f, 2.0f);
	    }
	}

	public StsBoundingBox getBoundingBox()
	{
		if(boundingBox != null) return boundingBox;
		boundingBox = new StsBoundingBox();
		int nSectionEdges = sectionEdges.getSize();
		for(int n = 0; n < nSectionEdges; n++)
		{
			StsSectionEdge sectionEdge = (StsSectionEdge)sectionEdges.getElement(n);
			StsPoint[] points = sectionEdge.getPoints();
			boundingBox.addPoints(points);
		}
		float stimulatedRadius = fractureSet.getStimulatedRadius();
		boundingBox.addBorderXY(1.1f*stimulatedRadius);
		return boundingBox;
	}

	public boolean intersects(StsGolderFracture fracture)
	{
		StsPoint[][] fractureLines = fracture.getTopAndBottomLines();
		StsPoint[] topLine = fractureLines[0];
		StsPoint[] botLine = fractureLines[1];
		if(intersectsLine(topLine)) return true;
		return intersectsLine(botLine);
	}

	public boolean intersectsLine(StsPoint[] line)
	{
		StsPoint[] sectionEdgePoints = getSectionIntersectionAtZ(line[0].getZ());
		if(sectionEdgePoints == null) return false;
		StsPoint point1 = sectionEdgePoints[0];
		for(int n = 1; n < sectionEdgePoints.length; n++)
		{
			StsPoint point0 = point1;
			point1 = sectionEdgePoints[n];
			if(StsMath.lineIntersectXY(line[0], line[1], point0, point1)) return true;
		}
		return false;
	}
}
