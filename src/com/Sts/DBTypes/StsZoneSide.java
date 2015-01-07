

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
import java.util.*;

public class StsZoneSide
{
    transient protected StsZone zone;
    transient protected StsZoneBlock zoneBlock;
    transient protected StsSurfaceEdge topEdge, botEdge;
    transient protected StsLineZone prevLineZone, nextLineZone;
    transient protected StsEdgeLoop edgeLoop;
    transient protected StsSection section;
    transient protected int side;

    transient StsEdgeLoopRadialLinkGrid linkedGrid;
    transient boolean gridConstructed = false;


    // transient float[][][] zoneGrid;
    // transient float[][][][] subZoneGrids;
    // private transient int[][] rowCols;
    // private transient float[][] normals;

    transient int nRows, nCols, nSubZones;

//	transient private int surfDisplayListNum = 0; // display list number (> 0)
//	transient private int gridDisplayListNum = 0; // display list number (> 0)

    public static final int RIGHT = StsParameters.RIGHT;
    public static final int LEFT = StsParameters.LEFT;

    public static final int PLUS = StsParameters.PLUS;
    public static final int MINUS = StsParameters.MINUS;

    public StsZoneSide()
    {
    }

    public StsZoneSide(StsZone zone, StsZoneBlock zoneBlock, StsSurfaceEdge topEdge, StsSurfaceEdge botEdge)
    {
        this.zone = zone;
        this.zoneBlock = zoneBlock;
        this.topEdge = topEdge;
        this.botEdge = botEdge;
        this.section = topEdge.getSection();
        this.side = topEdge.getSide();
        section.addZoneSide(this);
        topEdge.setBotZoneSide(this);
        botEdge.setTopZoneSide(this);
        constructEdgeLoopLinks();
        constructLinkGrid();
        //computeGridBoundingBox();
        edgeLoop.constructPolygonsAndGridLines(this);
    }

    public int getSide() { return side; }
    public int getNRows() { return nRows; }
    public int getNCols() { return nCols; }
    public StsSection getSection() { return section; }
    public StsEdgeLoop getEdgeLoop() { return edgeLoop; }
    public StsEdgeLinkable getTopEdge() { return topEdge; }
    public StsZoneBlock getZoneBlock() { return zoneBlock; }
    public StsBlock getBlock() { return zoneBlock.getBlock(); }

    public StsGridBoundingBox computeGridBoundingBox()
    {
        return linkedGrid.computeGridBoundingBox(edgeLoop);
    }
    /*
    private float[] getValues(int dim)
    {
        float[] values = new float[nCols*nRows];
        int n=0;
        for(int sz = 0; sz < nSubZones; sz++)
        {
            float[][][] subZoneGrid = subZoneGrids[sz];
            int nLayers = subZoneGrid.length;
            if(sz == nSubZones-1) nLayers++;

            for(int i=0; i <= nLayers; i++)
                for(int j=0; j<nCols; j++, n++)
                    values[n] = subZoneGrid[i][j][dim];
        }
        return values;
    }
    public float[] getXValues() { return getValues(0); }
    public float[] getYValues() { return getValues(1); }
    public float[] getZValues() { return getValues(2); }
    */
    // initialization is done by StsZoneBlock.classInitialize
    public boolean initialize(StsModel model)
    {
        if(edgeLoop == null) return false;
//        if(edgeLoop == null) return true;
        if(!edgeLoop.constructEdgeLinks()) return false;
        constructLinkGrid();
        // computeGridBoundingBox();
        return edgeLoop.constructPolygonsAndGridLines(this);
    }

    /** For all points on and inside the edgeLoop for this zoneSide, compute a loopBoundingBox
     *  which bounds the XYZ and row-column range of this zoneSide.
     *  We first put all the edgeLoop points in the boundingBox.  If the section this zoneSide
     *  is on is not vertical, it may "belly" out and extend past the edgeLoop in X and/or Y,
     *  so we need to add all the interior points as well.
     */
    private StsGridBoundingBox computeGridBoundingBox(StsEdgeLoop edgeLoop)
    {
        StsGridBoundingBox gridBoundingBox = computeGridBoundingBox(edgeLoop);
        StsBlockGrid xyGrid = zoneBlock.getTopGrid();
        if(!section.isVertical())
            edgeLoop.addInsidePointsToGridBoundingBox(gridBoundingBox, section, xyGrid);       
        gridBoundingBox.computeRowColRanges(xyGrid);
        return gridBoundingBox;
    }
    /*
    public void reconstructSubZoneGrid(StsSubZone subZone)
    {
        if(subZoneGrids == null) return;
        int subZoneNumber = subZone.getSubZoneNumber();
        subZoneGrids[subZoneNumber] = null;
    }
    */
/*
    public void constructGrid()
    {
        StsObjectList botEdgePoints, topEdgePoints;
        StsGridSectionPoint botEdgePoint, topEdgePoint;
        float[] botXYZ, topXYZ;
        float[][] rowColF;
        int i, j;

        try
        {
            StsXYSurfaceGridable surface = botEdge.getSurface();
            if(!(surface instanceof StsModelSurface)) return;
            StsModelSurface modelSurface = (StsModelSurface)surface;
            StsZone zone = modelSurface.getZoneAbove();

            botEdgePoints = botEdge.getEdgePoints();
            if(botEdgePoints == null) return;
            int nBotEdgePoints = botEdgePoints.getSize();

            topEdgePoints = topEdge.getEdgePoints();
            if(topEdgePoints == null) return;
            int nTopEdgePoints = topEdgePoints.getSize();

            if(nBotEdgePoints != nTopEdgePoints)
            {
                StsException.systemError("Number of bottom and top edges points are different." +
                    " for: " + getLabel());
                return;
            }

            nCols = nTopEdgePoints;
            nSubZones = zone.getNSubZones();
            int nRows = zone.getNLayers()+1; // not used here, but set for other routines

            zoneGrid = new float[nSubZones+1][nCols][];
            rowCols = new int[nCols][];

            for(j = 0; j < nCols; j++)
            {
                topEdgePoint = (StsGridSectionPoint)topEdgePoints.getElement(j);
                zoneGrid[0][j] = topEdgePoint.getCoordinates();

                botEdgePoint = (StsGridSectionPoint)botEdgePoints.getElement(j);
                zoneGrid[nSubZones][j] = botEdgePoint.getCoordinates();

                rowCols[j] = botEdgePoint.getLowerLeftRowCol(section);

                float rowF = topEdgePoint.getRowF(null);
                float colF = topEdgePoint.getColF(null);
                topXYZ = zoneGrid[0][j];
                botXYZ = zoneGrid[nSubZones][j];
                for(i = 1; i < nSubZones; i++)
                {
                    float z = zoneBlock.getZoneGridZ(i, rowF, colF);
                    zoneGrid[i][j] = new float[] { botXYZ[0], botXYZ[1], z };
                }

                constructNormals(zoneGrid[0]);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneSide.constructGrid() failed.",
                e, StsException.WARNING);
        }
    }
*/
    /*
    public void constructSubZoneGrids()
    {
        subZoneGrids = new float[nSubZones][][][];
        StsObjectRefList subZones = zone.getSubZones();
        for(int sz = 0; sz < nSubZones; sz++)
        {
            StsSubZone subZone = (StsSubZone)subZones.getElement(sz);
            constructSubZoneGrid(subZone);
        }
    }

    private float[][][] getSubZoneGrid(int subZoneNumber)
    {
        if(subZoneNumber >= nSubZones) return null;
        if(subZoneGrids == null || subZoneGrids[subZoneNumber] == null)
            return constructSubZoneGrid(subZoneNumber);
        else
            return subZoneGrids[subZoneNumber];
    }

    private float[][][] constructSubZoneGrid(int subZoneNumber)
    {
        StsObjectRefList subZones = zone.getSubZones();
        StsSubZone subZone = (StsSubZone)subZones.getElement(subZoneNumber);
        return constructSubZoneGrid(subZone);
    }

    private float[][][] constructSubZoneGrid(StsSubZone subZone)
    {
        int subZoneNumber = subZone.getSubZoneNumber();
        int nLayers = subZone.getNLayers();
        float[][][] subZoneGrid = new float[nLayers+1][nCols][];

        if(subZoneGrids == null) subZoneGrids = new float[nSubZones][][][];
        subZoneGrids[subZoneNumber] = subZoneGrid;

        subZoneGrid[0] = zoneGrid[subZoneNumber];
        subZoneGrid[nLayers] = zoneGrid[subZoneNumber+1];
        if(nLayers <= 1) return subZoneGrid;
        int subZoneType = subZone.getSubZoneType();
        float thickness = subZone.getLayerThickness();
        constructSubZoneGrid(subZoneGrid, nLayers, nCols, subZoneType, thickness);
        return subZoneGrid;
    }

    private void constructSubZoneGrid(float[][][] subZoneGrid, int nLayers, int nCols,
                                      int subZoneType, float thickness)
    {
        float dZ, z;

        for(int j = 0; j < nCols; j++)
        {
            float x = subZoneGrid[0][j][0];
            float y = subZoneGrid[0][j][1];
            float topZ = subZoneGrid[0][j][2];
            float botZ = subZoneGrid[nLayers][j][2];

            switch(subZoneType)
            {
                case StsZone.SUBZONE_UNIFORM:
                    dZ = (botZ - topZ)/nLayers;
                    z = topZ + dZ;
                    for(int k = 1; k < nLayers; k++)
                    {
                        subZoneGrid[k][j] = new float[] { x, y, z };
                        z += dZ;
                    }
                    break;
                case StsZone.SUBZONE_OFFLAP:
                    dZ = thickness;
                    z = botZ - dZ;
                    for(int k = nLayers-1; k >= 1 ; k--)
                    {
                        z = Math.max(z, topZ);
                        subZoneGrid[k][j] = new float[] { x, y, z };
                        z -= dZ;
                    }
                    break;
                case StsZone.SUBZONE_ONLAP:
                    dZ = zone.getSubZoneThickness();
                    z = topZ + dZ;
                    for(int k = 1; k < nLayers; k++)
                    {
                        z = Math.min(z, botZ);
                        subZoneGrid[k][j] = new float[] { x, y, z };
                        z += dZ;
                    }
                    break;
            }
        }
    }

    private void deleteSubZoneGrids()
    {
        if(subZoneGrids == null) return;
        for(int sz = 0; sz < nSubZones; sz++)
            subZoneGrids[sz] = null;
        subZoneGrids = null;
    }
    */
    private void constructEdgeLoopLinks()
    {
        StsSurfaceVertex topPrevVertex, botPrevVertex, topNextVertex, botNextVertex;
        StsLine prevLine, nextLine;

        topPrevVertex = topEdge.getPrevVertex();
        botPrevVertex = botEdge.getPrevVertex();
        prevLine = topPrevVertex.getSectionLine();
        prevLineZone = prevLine.getLineZone(topPrevVertex, botPrevVertex, section, side, PLUS);

        topNextVertex = topEdge.getNextVertex();
        botNextVertex = botEdge.getNextVertex();
        nextLine = topNextVertex.getSectionLine();
        nextLineZone = nextLine.getLineZone(topNextVertex, botNextVertex, section, side, PLUS);

        edgeLoop = new StsEdgeLoop(section);
        if(section.isPlanar())
        {
            float[] normal = section.getPlanarNormal();
            edgeLoop.setNormal(normal, true);
        }

        if(side == RIGHT)
        {
            edgeLoop.addEdge(botEdge, PLUS);
            edgeLoop.addEdge(nextLineZone, MINUS);
            edgeLoop.addEdge(topEdge, MINUS);
            edgeLoop.addEdge(prevLineZone, PLUS);
        }
        else
        {
            edgeLoop.addEdge(botEdge, MINUS);
            edgeLoop.addEdge(prevLineZone, MINUS);
            edgeLoop.addEdge(topEdge, PLUS);
            edgeLoop.addEdge(nextLineZone, PLUS);
        }
        edgeLoop.constructEdgeLinks();
    }

    private void constructLinkGrid()
    {
        linkedGrid = new StsEdgeLoopRadialLinkGrid(section, edgeLoop);
        //linkedGrid.addEdgeLoop(edgeLoop);
        // edgeLoop.addToLinkedGrid(linkedGrid);
        //linkedGrid.resetMinMaxLimits();
        //linkedGrid.orderLinks();
    }

    public void deleteDisplayLists(GL gl)
    {
        edgeLoop.deleteDisplayLists(gl);
//        gridConstructed = false;
    }
/*
    public void deleteGridDisplayList(GL gl)
    {
        if(gridDisplayListNum <= 0) return;

        gl.glDeleteLists(gridDisplayListNum, 1 );
        gridDisplayListNum = 0;
    }

    public void constructNonVerticalGrid()
    {
        deleteDisplayLists(currentModel.win3d.glPanel3d.getGL());
        edgeLoop.constructGrid();
        gridConstructed = true;
    }
*/
/*
    public void constructVerticalGrid()
    {
        StsObjectList botEdgePoints, topEdgePoints;
        StsGridSectionPoint botEdgePoint, topEdgePoint;
        float[] botXYZ, topXYZ;
        float[][] rowColF;
        int i, j;

        try
        {
            StsZone zone = botEdge.getSurface().getZoneAbove();

            botEdgePoints = botEdge.getEdgePoints();
            if(botEdgePoints == null) return;
            int nBotEdgePoints = botEdgePoints.getSize();

            topEdgePoints = topEdge.getEdgePoints();
            if(topEdgePoints == null) return;
            int nTopEdgePoints = topEdgePoints.getSize();

            if(nBotEdgePoints != nTopEdgePoints)
            {
                StsException.systemError("Number of bottom and top edges points are different." +
                    " for: " + getLabel());
                return;
            }

            nCols = nTopEdgePoints;
            nSubZones = zone.getNSubZones();
            int nRows = zone.getNLayers()+1; // not used here, but set for other routines

            zoneGrid = new float[nSubZones+1][nCols][];
            rowCols = new int[nCols][];

            for(j = 0; j < nCols; j++)
            {
                topEdgePoint = (StsGridSectionPoint)topEdgePoints.getElement(j);
                zoneGrid[0][j] = topEdgePoint.getCoordinates();

                botEdgePoint = (StsGridSectionPoint)botEdgePoints.getElement(j);
                zoneGrid[nSubZones][j] = botEdgePoint.getCoordinates();

                rowCols[j] = botEdgePoint.getLowerLeftRowCol(section);

                float rowF = topEdgePoint.getRowF(null);
                float colF = topEdgePoint.getColF(null);
                topXYZ = zoneGrid[0][j];
                botXYZ = zoneGrid[nSubZones][j];
                for(i = 1; i < nSubZones; i++)
                {
                    float z = zoneBlock.getZoneGridZ(i, rowF, colF);
                    zoneGrid[i][j] = new float[] { botXYZ[0], botXYZ[1], z };
                }

                constructNormals(zoneGrid[0]);
                gridConstructed = true;
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneSide.constructGrid() failed.",
                e, StsException.WARNING);
        }
    }
*/
    public void debugDisplay(StsWin3d win3d)
    {
/*
        StsEdgeLoopRadialLinkGrid gridLinkedGrid = edgeLoop.getGridLinkedGrid();
        if(gridLinkedGrid == null) return;
        gridLinkedGrid.debugDisplay(currentModel, win3d);
*/
    }

    private void debugDisplayEdgePoints(StsObjectList edgePoints, StsGLPanel3d glPanel3d)
    {
        GL gl = glPanel3d.getGL();
        if(gl == null) return;

        gl.glDisable(GL.GL_LIGHTING);
        glPanel3d.setViewShift(gl, StsGraphicParameters.vertexShift);

        int nPnts = edgePoints.getSize();
        for(int n = 0; n < nPnts; n++)
        {
            StsGridSectionPoint edgePoint = (StsGridSectionPoint)edgePoints.getElement(n);
            edgePoint.debugDisplayGridPoint(glPanel3d);
        }

        glPanel3d.resetViewShift(gl);
        gl.glEnable(GL.GL_LIGHTING);
    }



    public void display(StsModel model, StsGLPanel3d glPanel3d, String displayMode, boolean displayLayers, boolean displayProperties, StsColor stsColor, StsPropertyVolume propertyVolume)
    {
        if(section.getType() == StsParameters.BOUNDARY && side == LEFT) return;
        if(!section.getIsVisible()) return;

        if(!zone.getIsVisible() ) return;
        int sectionSide = section.getDisplaySide();
        if(side != sectionSide) return;

        if(!edgeLoop.subPolygonsConstructed)
        {
            edgeLoop.deleteDisplayLists(glPanel3d.getGL());
            constructLayerPolygons();
        }

        boolean displayFill = zone.getDisplayFill();
        boolean displayGrid = zone.getDisplayGrid();
        edgeLoop.checkSetPropertyChanged(displayMode);
 //       boolean useAdjustedPoints = currentModel.project.getUseAdjustedPoints();
        edgeLoop.display(model, glPanel3d, displayFill, displayGrid, displayLayers, displayProperties, stsColor, propertyVolume, true);
    }

    private StsColor getBlockColor()
    {
        return zoneBlock.getStsColor();
    }
/*
    public void display(StsGLPanel3d glPanel3d, StsColorScheme colorScheme)
    {
        if(section.getType() == StsParameters.BOUNDARY && side == LEFT) return;

        if(!zone.getIsVisible() ) return;
        int sectionSide = section.getDisplaySide();
        if(side != sectionSide) return;

        if(currentModel.getDisplaySubZoneColors() && !subPolygonsConstructed) constructSubPolygons();

        boolean displayFill = zone.getDisplayFill();
        boolean displayGrid = zone.getDisplayGrid();
        boolean displayGaps = true;
        edgeLoop.display(glPanel3d, displayFill, displayGrid, displayGaps, colorScheme);
    }
*/
    private void constructLayerPolygons()
    {
        edgeLoop.constructZoneSideLayerPolygons(this);
    }

    public boolean checkConstructLayerPolygons()
    {
        if(edgeLoop == null) return true;
        return edgeLoop.constructZoneSideLayerPolygons(this);
    }
    public StsGridBoundingBox getSectionBoundingBox()
    {
        return new StsGridBoundingBox(0, section.nRows-1, 0, section.nCols-1, false);
    }
/*
    public void addBlockGridPolygons(byte geometryType)
    {
        Iterator gridPolygonIterator = getGridPolygonIterator();
        StsBlock block = zoneBlock.getBlock();
        while(gridPolygonIterator.hasNext())
        {
            StsPolygon polygon = (StsPolygon)gridPolygonIterator.next();
            block.addCellColumnZoneSidePolygon(polygon, geometryType);
        }
    }
*/
    private Iterator getGridPolygonIterator()
    {
        return edgeLoop.getZoneSideGridPolygonIterator();
    }
/*
    private void drawNonVerticalZoneSurface(GL gl)
    {
        if(!gridConstructed) constructNonVerticalGrid();

        StsColor color = topEdge.getSurface().getStsColor();
        StsColor.setGLColor(gl, color);
        edgeLoop.displaySurface(model, gl, currentModel);
//        if(side == LEFT) gl.glFrontFace(GL.GL_CCW);
    }

    private void displayZoneSurface(GL gl, boolean displaySubZoneColors, boolean useDisplayLists)
    {
		//timer.start();

        try
        {
            if (surfDisplayListNum == 0 && useDisplayLists)  // build display list
            {
                surfDisplayListNum = gl.glGenLists(1);
                if(surfDisplayListNum == 0)
                {
                    currentModel.model.statusBarMessage("System Error in StsZoneSide.displaySurface: " +
                                           "Failed to allocate a display list");
                    return;
                }

                gl.glNewList(surfDisplayListNum, GL.GL_COMPILE_AND_EXECUTE);
                drawZoneSurface(gl, displaySubZoneColors);
                gl.glEndList();
                //timer.stop("display list surface setup: ");

            }
            else if(useDisplayLists) // use existing display list
            {
                gl.glCallList( surfDisplayListNum );
                //timer.stop("display list surface draw: ");
            }
            else // immediate mode draw
            {
                if(surfDisplayListNum > 0)
                {
                    gl.glDeleteLists(surfDisplayListNum, 1 );
                    surfDisplayListNum = 0;
                }
                drawZoneSurface(gl, displaySubZoneColors);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneSide.displayZoneSurface() failed.",
                                         e, StsException.WARNING);
        }
    }

    private void drawZoneSurface(GL gl, boolean displaySubZoneColors)
    {
        StsPropertyVolumeOld pv = zone.getPropertyVolume();

        if(!gridConstructed)
        {
            constructNonVerticalGrid();
            if(!gridConstructed) return;
        }
        drawNonVerticalZoneSurface(gl);
    }
*/
/*
    private void drawVerticalZoneSurface1(GL gl, boolean displaySubZoneColors)
    {
        float[][] topXYZs, botXYZs;

        try
        {
            if(zoneGrid == null) constructVerticalGrid();

            if(!displaySubZoneColors)
            {
                topXYZs = zoneGrid[0];
                botXYZs = zoneGrid[nSubZones];

                StsColor color = zone.getStsColor();
		        StsColor.setGLColor(gl, color);
                gl.glBegin(GL.GL_TRIANGLE_STRIP);
                for(int n = 0; n < nCols; n++)
                {
                    gl.glNormal3fv(normals[n], 0);
                    gl.glVertex3fv(topXYZs[n]);
                    gl.glVertex3fv(botXYZs[n]);
                }
                gl.glEnd();
            }
            else
            {
                StsObjectRefList subZones = zone.getSubZones();
                for(int sz = 0; sz < nSubZones; sz++)
                {
                    StsSubZone subZone = (StsSubZone)subZones.getElement(sz);
                    StsColor color = subZone.getStsColor();
		            StsColor.setGLColor(gl, color);

                    float[][][] subZoneGrid = getSubZoneGrid(sz);
                    if(subZoneGrid == null) continue;

                    int nLayers = subZoneGrid.length;
                    botXYZs = subZoneGrid[0];
                    for(int k = 0; k < nLayers-1; k++)
                    {
                        topXYZs = botXYZs;
                        botXYZs = subZoneGrid[k+1];
                        gl.glBegin(GL.GL_TRIANGLE_STRIP);
                        for(int n = 0; n < nCols; n++)
                        {
                            gl.glNormal3fv(normals[n], 0);
                            gl.glVertex3fv(topXYZs[n]);
                            gl.glVertex3fv(botXYZs[n]);
                        }
                        gl.glEnd();
                    }
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneSide.displayVerticalZoneSurface() failed.",
                                         e, StsException.WARNING);
            gl.glEnd();
        }
    }
*/
    /*
    private void constructNormals(float[][] xYZs)
    {
        int nPnts = xYZs.length;
        normals = new float[nPnts][];
        for(int n = 0; n < nPnts-1; n++)
            normals[n] = StsMath.horizontalNormal(xYZs[n], xYZs[n+1], 1);
        normals[nPnts-1] = normals[nPnts-2];
    }
    */
   /*
    private void drawVerticalZoneSurface1(GL gl, StsPropertyVolumeOld pv)
    {
        float[][] botXYZs, topXYZs;
        try
        {
            if(subZoneGrids == null) subZoneGrids = new float[nSubZones][][][];

            int nZoneLayer = 0;
            for(int sz = 0; sz < nSubZones; sz++)
            {
                float[][][] subZoneGrid = getSubZoneGrid(sz);
                if(subZoneGrid == null) continue;
                int nLayers = subZoneGrid.length-1;
                float[][] values = pv.getValuesArray(rowCols, nZoneLayer, nLayers);
                botXYZs = subZoneGrid[0];
                for(int k = 0; k < nLayers; k++, nZoneLayer++)
                {
                    topXYZs = botXYZs;
                    botXYZs = subZoneGrid[k+1];
                    gl.glBegin(GL.GL_QUAD_STRIP);
                    for(int n = 0; n < nCols-1; n++)
                    {
                        pv.getScaledColor(values[n][k]).setGLColor(gl);
//                        StsColor.setGLColor(gl, pv.getScaledColor(rowCols[n][0], rowCols[n][1], k));

                        gl.glNormal3fv(normals[n], 0);
                        gl.glVertex3fv(topXYZs[n], 0);
                        gl.glVertex3fv(botXYZs[n], 0);

                        gl.glNormal3fv(normals[n+1], 0);
                        gl.glVertex3fv(topXYZs[n+1], 0);
                        gl.glVertex3fv(botXYZs[n+1], 0);
                    }
                    gl.glEnd();
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsZoneSide.displayVerticalZoneSurface() failed.",
                                         e, StsException.WARNING);
            gl.glEnd();
        }
    }
    */
/*
    private void displayVerticalZoneGrid(StsWin3d win3d, GL gl, boolean useDisplayLists)
    {
		//timer.start();

        if (gridDisplayListNum == 0 && useDisplayLists)  // build display list
		{
            gridDisplayListNum = gl.glGenLists(1);
            if(gridDisplayListNum == 0)
            {
                currentModel.logMessage("System Error in StsGrid.displayGrid: " +
                                       "Failed to allocate a display list");
                return;
            }

            win3d.glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);
            gl.glNewList(gridDisplayListNum, GL.GL_COMPILE_AND_EXECUTE);
            drawVerticalZoneGrid(win3d, gl);
			gl.glEndList();
            win3d.glPanel3d.resetViewShift(gl);

			//timer.stop("display list surface setup: ");

		}
		else if(useDisplayLists) // use existing display list
		{
            win3d.glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);
			gl.glCallList( gridDisplayListNum );
            win3d.glPanel3d.resetViewShift(gl);
 			//timer.stop("display list surface draw: ");
		}
        else
        {
            if(gridDisplayListNum > 0)
            {
                gl.glDeleteLists(gridDisplayListNum, 1 );
                gridDisplayListNum = 0;
            }
            win3d.glPanel3d.setViewShift(gl, StsGraphicParameters.gridShift);
            drawVerticalZoneGrid(win3d, gl);
            win3d.glPanel3d.resetViewShift(gl);
        }
    }
*/
/*
    private void drawVerticalZoneGrid(StsWin3d win3d, GL gl)
    {
        float[] botXYZ = null;
        float[] topXYZ = null;

        try
        {
            gl.glDisable(GL.GL_LIGHTING);
            gl.glLineWidth(StsGraphicParameters.gridLineWidth);

            StsObjectList topEdgePoints = topEdge.getEdgePoints();

            gl.glBegin(GL.GL_LINES);
            for(int j = 0; j < nCols; j++)
            {
                topXYZ = zoneGrid[0][j];
                botXYZ = zoneGrid[nSubZones][j];
                StsGridSectionPoint topEdgePoint = (StsGridSectionPoint)topEdgePoints.getElement(j);
                if(topEdgePoint.isRowOrCol(null)) // don't draw section col lines
                {
                    gl.glVertex3fv(topXYZ);
                    gl.glVertex3fv(botXYZ);
                }
            }
            gl.glEnd();

            if(subZoneGrids == null) subZoneGrids = new float[nSubZones][][][];

            float[][][] subZoneGrid;
            for(int sz = 0; sz < nSubZones; sz++)
            {
                subZoneGrid = getSubZoneGrid(sz);
                if(subZoneGrid == null) continue;
                int nLayers = subZoneGrid.length;
                if(sz < nSubZones-1) nLayers--;

                for(int k = 0; k < nLayers; k++)
                {
                    gl.glBegin(GL.GL_LINE_STRIP);
                    for(int j = 0; j < nCols; j++)
                        gl.glVertex3fv(subZoneGrid[k][j]);
                    gl.glEnd();
                }
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsBlockSide.drawVerticalZoneGrid() failed.",
                                         e, StsException.WARNING);
            gl.glEnd();
        }
        finally
        {
            gl.glEnable(GL.GL_LIGHTING);
        }
    }
*/
    public void drawGridLines(int rowOrCol, int rowCol, GL gl)
    {
        if(edgeLoop == null) return;
        edgeLoop.drawGridLines(rowOrCol, rowCol, gl);
    }

    public byte getGeometryType()
    {
        return section.getGeometryType();
    }

    public String getLabel()
    {
        return new String(zone.getLabel() + section.getLabel());
    }

    public String toString()
    {
        return new String(StsObject.toString(zone) + " " + StsObject.toString(section) + StsParameters.sideLabel(side));
    }

    public void propertyTypeChanged()
    {
        edgeLoop.propertyTypeChanged();
    }
}
