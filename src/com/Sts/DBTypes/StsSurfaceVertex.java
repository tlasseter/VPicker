
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;

public class StsSurfaceVertex extends StsObject implements Cloneable, Comparable<StsSurfaceVertex>
{
    StsObjectList edges;             /** @param edges edges connected to this vertex */
    StsXYSurfaceGridable surface;                 /** @param surface surface associated (optional) */
    protected StsGridSectionPoint surfacePoint;   /** point and normal at this vertex */
    protected StsObject association;    /** @param association: well or edge         */
	protected StsBlock block;           /** block this vertex belongs to */
//    protected boolean isFault = false; /** @param isFault indicates in faultVertices list for well */

    // Convenience copies of various flags
    static final int PLUS = StsParameters.PLUS;
    static final int MINUS = StsParameters.MINUS;
    static final float nullValue = StsParameters.nullValue;
	static final int dotWidth = StsGraphicParameters.vertexDotWidth;

// Constructors

	public StsSurfaceVertex()
	{
	}

    public StsSurfaceVertex(boolean persistent)
	{
        super(persistent);
    }
/*
	public StsSurfaceVertex(StsPoint point, StsWell well, StsSurface surface)
	{
        this(point, well, surface, false);
    }
*/
	// construct for a sectionEdge vertex; vertex is stored in wellLine vertices
	public StsSurfaceVertex(StsPoint point, StsLine line, StsXYSurfaceGridable surface, StsBlock block, boolean persistent)
	{
        super(persistent);
        this.association = line;
        this.surface = surface;
		this.block = block;
        surfacePoint = new StsGridSectionPoint(point, line, this, persistent);
	}

	// constructor for surfaceEdge vertex; vertex is stored in one of well vertices list depending on isFault
/*
	public StsSurfaceVertex(StsPoint point, StsWell well, StsSurface surface, boolean isFault)
	{
		this(point, well, surface);
		this.isFault = isFault;
	}
*/
    public StsSurfaceVertex(StsGridSectionPoint surfacePoint)
    {
        this.surfacePoint = surfacePoint;
        surfacePoint.setVertex(this);
    }

    public StsSurfaceVertex(StsGridSectionPoint surfacePoint, StsObject association)
    {
        this(surfacePoint, association, true);
    }

    public StsSurfaceVertex(StsGridSectionPoint surfacePoint, StsObject association, boolean persistent)
    {
        super(persistent);
        this.surfacePoint = surfacePoint;
        surfacePoint.setVertex(this);
        this.association = association;
    }

	public StsSurfaceVertex(StsPoint point, StsLine line)
	{
		this.association = line;
        surfacePoint = new StsGridSectionPoint(point, line, this);
    }

    public boolean initialize(StsModel model) { return true; }

    public StsSurfaceVertex copy()
    {
        StsGridSectionPoint newSurfacePoint = surfacePoint.copy();
        return new StsSurfaceVertex(newSurfacePoint);
    }

    public boolean delete()
    {
        edges = null;
        if(surfacePoint != null) surfacePoint.delete();
        return super.delete();
    }

// Accessors

//    public StsSurfaceEdge getPrevEdge() { return prevEdge; }
//    public StsSurfaceEdge getNextEdge() { return nextEdge; }
    public StsGridSectionPoint getSurfacePoint() { return surfacePoint; }
    public void setSurfacePoint(StsGridSectionPoint surfacePoint) { this.surfacePoint = surfacePoint; }
//    public void setPrevEdge(StsSurfaceEdge edge) { prevEdge = edge; }
//    public void setNextEdge(StsSurfaceEdge edge) { nextEdge = edge; }
    public float getSectionColF(StsSection section)
    {
        if(surfacePoint == null) return nullValue;
        else return surfacePoint.getColF(section);
    }

    public float getSectionRowF(StsSection section)
    {
        if(surfacePoint == null) return nullValue;
        else return surfacePoint.getRowF(section);
    }

    public float getGridColF()
    {
        if(surfacePoint == null) return nullValue;
        else return surfacePoint.getColF(null);
    }

    public float getGridRowF()
    {
        if(surfacePoint == null) return nullValue;
        else return surfacePoint.getRowF(null);
    }

	public float getColF(StsSurfaceGridable grid) { return surfacePoint.getColF(grid); }
	public float getRowF(StsSurfaceGridable grid) { return surfacePoint.getRowF(grid); }

	//    public void setColF(StsSection section, float colF) { getSurfacePoint().setColF(section, sectionF); }
//    public void setSectionF(StsSection section) { surfacePoint.setSectionF(section); }
    public void setSurface(StsXYSurfaceGridable surface) { this.surface = surface; }
    public StsXYSurfaceGridable getSurface() { return surface; }
    // public void setBlock(StsBlock block) { this.block = block; }
    public StsBlock getBlock() { return block; }

    public void addEdge(StsEdgeLinkable edge)
    {
        if(edges == null) edges = new StsObjectList(2, 2);
        if(edges.contains((StsObject)edge))
        {
            StsException.systemError("StsSurfaceVertex.addEdge(): " + getLabel() + " already has this edge: " + edge.getLabel());
            return;
        }
        edges.add((StsObject)edge);
    }

    public StsEdgeLinkable getNextEdge()
    {
        int nEdges;

        if(edges == null || (nEdges = edges.getSize()) == 0) return null;
        for(int n = 0; n < nEdges; n++)
        {
            StsEdgeLinkable edge = (StsEdgeLinkable)edges.getElement(n);
            if(edge.getPrevVertex() == this) return edge;
        }
        return null;
    }

    public StsEdgeLinkable getPrevEdge()
    {
        int nEdges;

        if(edges == null || (nEdges = edges.getSize()) == 0) return null;
        for(int n = 0; n < nEdges; n++)
        {
            StsEdgeLinkable edge = (StsEdgeLinkable)edges.getElement(n);
            if(edge.getNextVertex() == this) return edge;
        }
        return null;
    }

	public void setAssociation(StsObject object) { this.association = object; }


    public StsLine getSectionLine()
    {
        if(association instanceof StsLine)
            return (StsLine)association;
        else
            return null;
    }
/*
    public StsWell getWell()
    {
        if(association instanceof StsWell)
            return (StsWell)association;
        else
            return null;
    }
*/
	public void setPoint(StsPoint point)
	{
		if(surfacePoint == null) return;
		surfacePoint.setPoint(point);
	}

    public StsPoint getPoint()
    {
        if(surfacePoint == null) return null;
        else
            return surfacePoint.getPoint();
    }

    /** For this vertex, check if it is on a well(fault) which is either not
     *  on a section or is on a non-fault section.
     */
    public boolean onDyingFault()
    {
        StsLine line = getSectionLine();
        if(line == null) return false;
        return line.isDyingFault();
    }

    public boolean onFault()
    {
        StsLine line = getSectionLine();
        if(line == null) return false;
        return line.isFault();
    }

    public StsEdgeLinkable getPrevConnectedSurfaceEdge(StsEdgeLinkable edge)
    {
        if(edges == null) return null;
//        if(onDyingFault()) return null;
        int nEdges = edges.getSize();
        for(int n = 0; n < nEdges; n++)
        {
            StsEdgeLinkable prevEdge = (StsEdgeLinkable)edges.getElement(n);
            if(edge == prevEdge) continue;
            if(prevEdge.getNextVertex() == this) return prevEdge;
        }
        return null;
    }

    public StsEdgeLinkable getNextConnectedSurfaceEdge(StsEdgeLinkable edge)
    {
        if(edges == null) return null;
//        if(onDyingFault()) return null;
        int nEdges = edges.getSize();
        for(int n = 0; n < nEdges; n++)
        {
            StsEdgeLinkable nextEdge = (StsEdgeLinkable)edges.getElement(n);
            if(nextEdge == edge) continue;
            if(nextEdge instanceof StsSurfaceEdge && nextEdge.getPrevVertex() == this) return nextEdge;
        }
        return null;
    }

    public StsSurfaceEdge getOtherSurfaceEdge(StsEdgeLinkable edge)
    {
        if(edges == null) return null;
        int nEdges = edges.getSize();
        for(int n = 0; n < nEdges; n++)
        {
            StsEdgeLinkable otherEdge = (StsEdgeLinkable)edges.getElement(n);
            if(otherEdge == edge) continue;
            if(otherEdge instanceof StsSurfaceEdge)
            {
                if(otherEdge.getPrevVertex() == this) return (StsSurfaceEdge)otherEdge;
                if(otherEdge.getNextVertex() == this) return (StsSurfaceEdge)otherEdge;
            }
        }
        return null;
    }

	public StsSectionEdge getSectionEdge()
	{
		if(association instanceof StsSectionEdge) return (StsSectionEdge)association;
		else return null;
	}

	public void setSectionEdge(StsEdge sectionEdge)
	{
		if(association != null)
		{
			StsException.systemError("StsSurfaceVertex.setSectionEdge() failed." +
				" surfaceVertex: " + getLabel() + " already has association: " + association.getLabel());
			return;
		}
		association = sectionEdge;
	}

	public void setSectionLine(StsLine sectionLine)
	{
    /*
		if(association != null)
		{
			StsException.systemError("StsSurfaceVertex.setSectionEdge() failed." +
				" surfaceVertex: " + getLabel() + " already has association: " + association.getLabel());
			return;
		}
    */
		association = sectionLine;
	}

	public StsObject getAssociatedSectionOrLine()
	{
		if(association == null)
			return null;
		else if(association instanceof StsLine)
			return association;
		else if(association instanceof StsSectionEdge)
			return ((StsSectionEdge)association).getSection();
		else
			return null;
	}

	public StsSection getAssociatedSection()
	{
		if(association instanceof StsSectionEdge)
			return ((StsSectionEdge)association).getSection();
		else
			return null;
	}

    public void deleteEdgeFromVertex(StsEdgeLinkable edge)
    {
        if(edges == null) return;
        int nEdges = edges.getSize();
        for(int n = 0; n < nEdges; n++)
        {
            if(edge == (StsEdgeLinkable)edges.getElement(n))
                edges.delete((StsObject)edge);
        }

        if(edges.getSize() == 0)
        {
            if(association != null && association instanceof StsLine)
            {
                StsLine line = (StsLine) association;
                line.deleteSectionVertex(this);
            }
            surfacePoint.delete();
            delete();
        }
    }


    public boolean isConnectedToEdge(StsEdge edge)
    {
        return hasEdge(edge) || isOnEdge(edge);
    }

   public boolean hasEdge(StsEdge edge)
    {
        if(edges == null) return false;
        int nEdges = edges.getSize();
        for(int n = 0; n < nEdges; n++)
        {
            StsEdge connectedEdge = (StsEdge)edges.getElement(n);
            if(connectedEdge == edge) return true;
        }
        return false;
    }

    public StsEdge getSurfaceEdge(StsSurface surface)
    {
        if(edges == null) return null;
        int nEdges = edges.getSize();
        for(int n = 0; n < nEdges; n++)
        {
            StsEdge edge = (StsEdge)edges.getElement(n);
            if(edge.surface == surface) return edge;
        }
        return null;
    }

    public boolean isOnEdge()
    {
		return association instanceof StsSectionEdge;
    }

    public boolean isOnLine()
    {
		return association instanceof StsLine;
    }

    public boolean isOnEdge(StsEdge edge)
    {
		return edge == association;
    }

    public boolean isConnectedToSection(StsSection section)
    {
		if(section == null) return false;
		return onSection() == section;
    }

    public StsSection onSection()
    {
		if(association instanceof StsSectionEdge)
			return ((StsSectionEdge)association).getSection();
		else if(association instanceof StsLine)
			return ((StsLine)association).getOnSection();
		else
			return null;
    }

    public boolean isOnSection()
    {
        return onSection() != null;
    }

    public boolean isFullyConnected()
    {
        if(edges == null) return false;
        int nEdges = edges.getSize();
        if(nEdges == 0) return false;
        if(nEdges >= 2) return true;
		// nEdges == 1
		return onSection() != null;
    }

    public StsEdge getOnlyConnectedEdge()
    {
        if(edges == null) return null;
        if(edges.getSize() != 1) return null;
        return (StsEdge)edges.getElement(0);
    }

    public boolean adjustToGrid()
    {
        return adjustToGrid(surface);
    }
    public boolean adjustToGrid(StsXYSurfaceGridable grid)
    {
        StsPoint point;

        StsLine line = getSectionLine();
        if(line == null) return false;

        if(line.isVertical)
        {
            point = surfacePoint.getPoint();
            grid.interpolateBilinearZ(point, true, true);
        }
        else
        {
            StsGridPoint gridPoint;
            float startZ = surfacePoint.getPoint().getZorT();

            if(onDyingFault() || grid == null)
                gridPoint = line.computeGridIntersect(surface, startZ);
            else
                gridPoint = line.computeGridIntersect(grid, startZ);

            if(gridPoint == null) return false;
            point = gridPoint.getPoint();
            if(point == null) return false;
            // if surface extends below surface below, truncate against surface below, setting z to 0.5 units above surface below
            StsSurfaceVertex vertexBelow = line.getVertexBelow(this);
            if(vertexBelow != null)
            {
                float z = point.getZ();
                float zBelow = vertexBelow.getPoint().getZorT();
                if(z >= zBelow)
                {
                    float adjustedZ = zBelow - 0.5f;
                    point = line.getXYZPointAtZorT(adjustedZ, true);
                }
            }
        }
        float z = point.getZorT();
        point.setZorT(StsMath.checkRoundOffInteger(z));
        surfacePoint.setPoint(point);
        return true;
    }

    public void checkAddTimeOrDepth()
    {
        try
        {
            StsPoint point = getPoint();
            float x = point.getX();
            float y = point.getY();
            float zt = point.getZorT();
            point.checkExtendVector(5);
            StsSeismicVelocityModel velocityModel = currentModel.getProject().getSeismicVelocityModel();

            if (isDepth)
            {
                if(velocityModel != null)
                {
                    float t = (float) velocityModel.getT(x, y, zt, zt);
                    point.setT(t);
                }
            }
            else
            {
                if(velocityModel != null)
                {
                    float z = (float) velocityModel.getZ(x, y, zt);
                    point.setZ(z);
                }
                point.setT(zt);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsSurfaceVertex.checkAddTimeOrDepth() failed.", e, StsException.WARNING);
        }
    }
/*
    public boolean connectedOK()
    {
        if(association instanceof StsEdge) return true;
		StsWell well = getWell();
        if(well != null && well.getSection() != null) return true;
        if(edges != null && edges.getSize() == 2) return true;
        return false;
    }
*/
    public void replaceConnectedEdge(StsSectionEdge oldEdge, StsSectionEdge newEdge)
    {
        if(edges == null)
		    edges = new StsObjectList(2, 2);
        else
            edges.delete(oldEdge);
        edges.add(newEdge);
    }

	// if vertex is already associated with a well we can't add association
	// but an associated edge is ok if well is not already attached to a section
    public boolean addEdgeAssociation(StsEdge pickedEdge)
    {
		if(association == null)
		{
            association = pickedEdge;
            return true;
		}
		return onSection() == null;
    }

	public void deleteEdge(StsEdge pickedEdge)
	{
		if(association == pickedEdge)
		{
			pickedEdge = null;
			return;
		}
		if(edges == null) return;
		int nEdges = edges.getSize();
		for(int n = 0; n < nEdges; n++)
		{
			if(edges.getElement(n) == pickedEdge)
			{
				edges.delete(n);
				return;
			}
		}
	}

    public StsEdge getOtherEdge(StsEdge edge)
    {
        StsEdge connectedEdge;

        if(edges == null) return null;

        int nEdges = edges.getSize();

        if(nEdges == 0)
        {
            errorEdgeNotConnected(edge);
            return null;
        }
        if(nEdges == 1)
        {
            if(edge != (StsEdge)edges.getFirst()) errorEdgeNotConnected(edge);
            return null;
        }
        else
        {
            if(edge == (StsEdge)edges.getFirst())
                return (StsEdge)edges.getLast();
            else if(edge == (StsEdge)edges.getLast())
                return (StsEdge)edges.getFirst();
            else
            {
                errorEdgeNotConnected(edge);
                return null;
            }
        }
    }

	public int getSectionPosition(StsSection section)
	{
		return surfacePoint.getSectionPosition(section);
	}

    private void errorEdgeNotConnected(StsEdge edge)
    {
        StsException.systemError("StsSurfaceVertex.getOtherEdge() failed. " +
                " Edge: " + edge.getLabel() + " not connected to vertex: " + getLabel());
    }

    public StsPoint getRotatedPoint(StsProject project, float dXOrigin, float dYOrigin)
    {
        StsPoint unrotatedPoint = getPoint();
        float[] xy = project.getRotatedRelativeXYFromUnrotatedRelativeXY(dXOrigin + unrotatedPoint.v[0], dYOrigin + unrotatedPoint.v[1]);
        StsPoint rotatedPoint = new StsPoint(5, unrotatedPoint);
        rotatedPoint.v[0] = xy[0];
        rotatedPoint.v[1] = xy[1];
        return rotatedPoint;
    }

    public void pick(GL gl, StsGLPanel glPanel)
    {
        StsGLDraw.drawPoint(getPoint().getXYZorT(), gl, dotWidth);
    }

    public String getLabel()
    {
		String associationLabel = "", blockLabel = "";

		if(association != null) associationLabel = association.getLabel();
		if(block != null) blockLabel = block.getLabel();
		return new String("vertex-" + getIndex() + " on: " + associationLabel + blockLabel);
    }

    public void display(StsGLPanel3d glPanel3d)
    {
		display(glPanel3d, StsColor.WHITE, dotWidth, StsGraphicParameters.vertexShift);
    }

    public void display(StsGLPanel3d glPanel3d, StsColor color)
    {
        display(glPanel3d, color, dotWidth, StsGraphicParameters.vertexShift);
    }

    public void display(StsGLPanel3d glPanel3d, double viewShift)
    {
        display(glPanel3d, StsColor.WHITE, dotWidth, viewShift);
    }

    public void display(StsGLPanel3d glPanel3d, StsColor color, int size)
    {
		display(glPanel3d, color, size, StsGraphicParameters.vertexShift);
    }

    public void display(StsGLPanel3d glPanel3d, double viewShift, int zIndex)
    {
        display(glPanel3d, StsColor.WHITE, dotWidth, viewShift, zIndex);
    }

    public void display(StsGLPanel3d glPanel3d, StsColor color, int size, double viewShift)
    {
        if (glPanel3d == null || size==0) return;
        float[] xyz = surfacePoint.getXYZorT(isDepth);
        StsGLDraw.drawPoint(xyz, (color==null) ? StsColor.WHITE : color, glPanel3d, size, size, viewShift);
    }

    public void display(StsGLPanel3d glPanel3d, StsColor color, int size, double viewShift, int zIndex)
    {
        if (glPanel3d == null || size==0) return;
        float[] xyz = surfacePoint.getXYZorT(isDepth);
        StsGLDraw.drawPoint(xyz, zIndex, (color==null) ? StsColor.WHITE : color, glPanel3d, size, viewShift);
    }

    public int compareTo(StsSurfaceVertex otherVertex)
    {
        float z = getPoint().getZ();
        float otherZ = otherVertex.getPoint().getZ();
        return Float.compare(z, otherZ);
    }

    public String toString() { return  "association: " + toString(association) + " point: " + getPoint().toString(); }
}










