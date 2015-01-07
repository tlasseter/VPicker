package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.StsClassFractureDisplayable;
import com.Sts.Interfaces.StsFractureDisplayable;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.util.*;

public class StsFractureCursorSection extends StsCursor3dVolumeTexture
{
	/** fracture instances displayable on cursor section */
	transient ArrayList<StsFractureDisplayable> displayableFractures = new ArrayList<StsFractureDisplayable>();
	/** bounding box around all fractures used for texture definition */
	transient StsRotatedGridBoundingBox fracturesBoundingBox;
    /** convenience copy of class manager */
    transient StsFractureSetClass fractureSetClass;
    /** index in volume of this cursorSection */
    transient public int nPlane = -1;
    /** texture data */
    transient byte[] data;
    /** subVolumed data (same size as data with some points removed */
    transient byte[] subVolumeData;
    transient boolean isDisplayingVolume = false;
	transient boolean isDisplayingRadius = true;
    transient boolean isDisplayingOnOtherCursorSections = false;
    transient boolean wasDragging = false;
    transient boolean isDragging = false;
    transient boolean displayVolume = false;
    transient boolean displayRadius = false;
    transient boolean displayOnOtherCursorSections = false;

    static final float largeFloat = StsParameters.largeFloat;
    /** colors used in drawing background and contour lines; fracture set colors follow these in the colorList. */
    static public final int nBasicColors = 2;
    static public final byte clearColorByte = 0;
    static public final byte lineColorByte = 1;

    static final boolean debug = false;

    public StsFractureCursorSection()
    {
    }

    public StsFractureCursorSection(StsModel model, StsCursor3d cursor3d, int dir)
    {
        initialize(model, cursor3d, dir);
    }

    public boolean initialize(StsModel model, StsCursor3d cursor3d, int dir)
    {
        super.initialize(model, cursor3d, dir);
		initializeDisplayableFractures();
        nullByte = 0;
        // isDisplayingOnOtherCursorSections = getDisplayedOnOtherCursorSections(glPanel3d);
//        initializeTextureTiles();
        return true;
    }

	private void initializeDisplayableFractures()
	{
		zDomain = StsProject.TD_DEPTH;
		StsArrayList fractureDisplayableClasses = model.getFractureDisplayableClasses();
		int nClasses = fractureDisplayableClasses.size();
		for(int n = 0; n < nClasses; n++)
		{
			StsClassFractureDisplayable stsClass = (StsClassFractureDisplayable)fractureDisplayableClasses.get(n);
			ArrayList<StsFractureDisplayable> objects = stsClass.getDisplayableFractures();
			displayableFractures.addAll(objects);
		}
		fracturesBoundingBox = new StsRotatedGridBoundingBox();
		for(StsFractureDisplayable displayableFracture : displayableFractures)
			fracturesBoundingBox.addUnrotatedBoundingBox(displayableFracture.getBoundingBox());
		model.getProject().adjustRotatedBoundingBoxGrid(fracturesBoundingBox, zDomain);
        fractureSetClass = (StsFractureSetClass)model.getStsClass(StsFractureSet.class);
	}

    public boolean isVisible()
    {
        if(!isVisible) return false;
        return fractureSetClass.displayStimulatedVolume();
    }

    protected void checkTextureAndGeometryChanges(StsGLPanel3d glPanel3d, boolean is3d)
    {
        StsCursor3d cursor3d = glPanel3d.getCursor3d();
        isDragging = cursor3d.getIsDragging(dirNo);
        displayVolume = (dirNo == ZDIR && fractureSetClass.isDisplayingVolume() && !isDragging);
        displayRadius = fractureSetClass.isDisplayingRadius() && !isDragging;
        displayOnOtherCursorSections = getDisplayedOnOtherCursorSections();
        if(debug) System.out.println("dir: " + dirNo + " isDragging: " + isDragging + " displayVolume " + displayVolume +
                " displayRadius " + displayRadius + " displayOnOtherCursorSections " + displayOnOtherCursorSections);

        if(displayVolume)
        {
            StsCursorSection cursorSection = cursor3d.cursorSections[dirNo];
            StsCursor3dTexture[] visibleCursorSections = cursorSection.getVisibleDisplayableSections();
            displayOnOtherCursorSections = visibleCursorSections.length > 1;
        }

        if(displayVolume != isDisplayingVolume)
        {
            textureChanged = true;
            isDisplayingVolume = displayVolume;
        }

        if(wasDragging)
        {
            if(!isDragging)
            {
                textureChanged = true;
                if(debug) System.out.println("dirNo: " + dirNo + " was dragging, now is not: clearing subVolumeData");
                subVolumeChanged();
                wasDragging = false;
            }
        }
        else
        {
            wasDragging = isDragging;
        }

        if(displayVolume != isDisplayingVolume)
        {
            textureChanged = true;
            isDisplayingVolume = displayVolume;
        }

        if(displayOnOtherCursorSections != isDisplayingOnOtherCursorSections)
        {
            textureChanged = true;
            isDisplayingOnOtherCursorSections = displayOnOtherCursorSections;
        }

        byte projectZDomain = StsObject.getCurrentModel().getProject().getZDomain();
        if (projectZDomain != zDomain)
        {
            if(debug) System.out.println("dir: " + dirNo + " domainChanged: changing texture, geometry, subVolume ");
            geometryChanged = true;
            textureChanged = true;
            zDomain = projectZDomain;
            subVolumeChanged();
        }
    }

    protected boolean enableGLState(StsGLPanel3d glPanel3d, GL gl, boolean is3d)
    {
        super.enableGLState(glPanel3d, gl, is3d);
        if(is3d)
        {
            gl.glDepthFunc(GL.GL_LEQUAL);
            glPanel3d.setViewShift(gl, 1.0f);
            if(dirNo == ZDIR)
            {
                gl.glEnable(GL.GL_POLYGON_STIPPLE);
                gl.glPolygonStipple(StsGraphicParameters.halftone, 0);
            }
        }
        return fractureSetClass.setColorList(gl, textureTiles.shader, displayOnOtherCursorSections);
    }

    protected void disableGLState(StsGLPanel3d glPanel3d, GL gl, boolean is3d)
    {
        super.disableGLState(glPanel3d, gl, is3d);
        if(is3d)
        {
            gl.glDepthFunc(GL.GL_LESS);
            glPanel3d.resetViewShift(gl);
            if(dirNo == ZDIR)
                gl.glDisable(GL.GL_POLYGON_STIPPLE);
        }
    }

    protected boolean initializeTextureTiles(StsGLPanel3d glPanel3d)
    {
        if(textureTiles == null)
        {
            if(debug) System.out.println("initialize textureTiles");
            StsCropVolume cropVolume = model.getProject().getCropVolume();
            textureTiles = StsTextureTiles.constructor(model, this, dirNo, fracturesBoundingBox, isPixelMode, cropVolume);
            if(textureTiles == null) return false;
            nTextureRows = textureTiles.nTotalRows;
            nTextureCols = textureTiles.nTotalCols;
            textureChanged = true;
        }
        else if (!textureTiles.isSameSize(fracturesBoundingBox))
        {
            if(debug) System.out.println("textureTiles size changed: rebuilding tiles");
            textureTiles.constructTiles(fracturesBoundingBox);
            nTextureRows = textureTiles.nTotalRows;
            nTextureCols = textureTiles.nTotalCols;
            textureChanged = true;
        }

        if (textureTiles.shaderChanged())
        {
            textureChanged = true;
        }
        textureTiles.setTilesDirCoordinate(dirCoordinate);

        return true;
     }

    protected void computeTextureData()
    {
        if(debug) System.out.println("Computing fractureCursorSection bytePlane for dir: " + dirNo + " dirCoor: " + dirCoordinate);
        data = computeBytePlaneData(displayVolume, displayRadius);

        if(debug) textureDebug(data, nTextureRows, nTextureCols, this, "displayTexture", "compute fracture texture data", nullByte);

        if(data == null)
        {
            textureChanged = false;
            return;
        }
        if(!isDragging)
        {
            subVolumeData = applySubVolume();
            if(subVolumeData == null)
            {
                StsException.systemError(this, "computeTextureData", "Failed to compute subVolumeData");
                subVolumeData = data;
            }
            if(debug) textureDebug(subVolumeData, nTextureRows, nTextureCols, this, "displayTexture", "not dragging: fracture texture subVolume applied", nullByte);
        }
        else
        {
            subVolumeData = data;
            if(debug) textureDebug(subVolumeData, nTextureRows, nTextureCols, this, "displayTexture", "dragging: fracture texture subVolume not applied", nullByte);
        }
    }

    protected void displayTiles3d(StsGLPanel3d glPanel3d, GL gl)
    {
        if(textureChanged)
        {
            if(debug) StsException.systemDebug(this, "displayTexture", "Displaying changed texture 3d in dirNo " + dirNo + ".");
            textureTiles.displayTiles(this, gl, isPixelMode, subVolumeData, nullByte);
        }
        else
        {
            if(debug) StsException.systemDebug(this, "displayTexture", "Displaying current texture 3d in dirNo " + dirNo + ".");
            textureTiles.displayTiles(this, gl, isPixelMode, (byte[])null, nullByte);
        }
    }

    protected void displayTiles2d(StsGLPanel3d glPanel3d, GL gl)
    {
        StsViewCursor viewCursor = (StsViewCursor)glPanel3d.getView();
        if(textureChanged)
            textureTiles.displayTiles2d(this, gl, viewCursor.axesFlipped, isPixelMode, subVolumeData, nullByte);
       	else
            textureTiles.displayTiles2d(this, gl, viewCursor.axesFlipped, isPixelMode, (byte[])null, nullByte);
    }

    private void initializeDataArray()
    {
        data = new byte[nTextureRows * nTextureCols];
    }

    public int getDefaultShader() { return StsJOGLShader.NONE; }
    public boolean getUseShader() { return false; }
/*
    private void initShader()
	{
		if(StsJOGLShader.canUseShader)
			shader = StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS;
		else
			shader = StsJOGLShader.NONE;
	}
*/
    /** This puts texture display on delete list.  Operation is performed at beginning of next draw operation. */
    public boolean textureChanged()
    {
        nPlane = -1;
        textureChanged = true;
        if(debug) StsException.systemDebug(this, "textureChanged");
        return true;
    }

    public boolean dataChanged()
    {
        data = null;
        textureChanged();
        return true;
    }

    private byte[] computeBytePlaneData(boolean displayVolume, boolean displayRadius)
    {
        if(!checkComputeBytePlaneData(displayVolume, displayRadius)) return null;
        return data;
    }

    private boolean checkComputeBytePlaneData(boolean displayVolume, boolean displayRadius)
    {
//        constructTemplate();
        initializeDataArray();
        if(dirNo == ZDIR)
        {
            if(!computeZBytePlaneData(fracturesBoundingBox)) return false;
            if(!displayRadius) return true;
            return computeDistanceTransformZBytePlaneData(fracturesBoundingBox);
        }
        else if(dirNo == XDIR)
            return computeXBytePlaneData(fracturesBoundingBox);
        else
            return computeYBytePlaneData(fracturesBoundingBox);
    }

    public int getNumberOfFilledCells(float z)
    {
    	dirCoordinate = z;
        computeZBytePlaneData(fracturesBoundingBox);
        if(data == null) return 0;
        byte[] subData = applySubVolume();
        if(subData == null) return 0;
    	int count = 0;
    	for(int i=0; i<subData.length; i++)
    		if(subData[i] != 0) count++;
    	return count;
    }


    private boolean computeDistanceTransformZBytePlaneData(StsRotatedGridBoundingBox fractureBoundingBox)
    {
        return computeDistanceTransformZBytePlaneData(dirCoordinate, null);
//        StsProgressBarDialog progressDialog = StsProgressBarDialog.constructor(null, "Distance Transform");
//        return computeDistanceTransformZBytePlaneData(fractureBoundingBox, dirCoordinate, progressDialog);
    }

    public boolean computeDistanceTransformZBytePlaneData(float z)
    {
        return computeDistanceTransformZBytePlaneData(z, null);
    }
    
    private boolean computeDistanceTransformZBytePlaneData(float z, StsProgressBarDialog progressDialog)
    {
        final StsDistanceTransform distanceTransform = new StsDistanceTransform(fracturesBoundingBox);
        boolean hasPoints = false;
		StsProject project = model.getProject();
		float rowGridSize = project.getRowGridSize();
        float colGridSize = project.getColGridSize();

        for(StsFractureDisplayable fracture : displayableFractures)
        {
            if(!fracture.getIsVisible()) continue;
			StsPoint[] zLinePoints = fracture.getSectionIntersectionAtZ(z);
			if(zLinePoints == null) continue;
			ArrayList<StsGridCrossingPoint> gridCrossingPoints = StsGridCrossings.computeGridCrossingPoints(zLinePoints, fracturesBoundingBox);
			int nCrossingPoints = gridCrossingPoints.size();
			if(nCrossingPoints > 0) hasPoints = true;
			for(int p = 0; p < nCrossingPoints; p++)
			{
				StsGridCrossingPoint gridPoint = gridCrossingPoints.get(p);
				int rowOrCol = gridPoint.rowOrCol;
				if(rowOrCol == StsParameters.ROW)
				{
					int row = (int)gridPoint.iF;
					int col = Math.round(gridPoint.jF);
					float distance = Math.abs(col - gridPoint.jF)*colGridSize;
					/*
					if(distance > stimulatedRadius/100.0f)
					{
					   System.out.println("Setting distance to 0.0, outside radius.");
					   distance = largeFloat;
					}
					*/
					distanceTransform.setDistance(row, col, distance);
				}
				else if(rowOrCol == StsParameters.COL)
				{
				   int row = Math.round(gridPoint.iF);
				   int col = (int)gridPoint.jF;
				   float distance = Math.abs(row - gridPoint.iF)*rowGridSize;
					/*
				   if(distance > stimulatedRadius/100.0f)
				   {
					   System.out.println("Setting distance to 0.0, outside radius.");
					   distance = largeFloat;
				   }
				   */
				   distanceTransform.setDistance(row, col, distance);
				}
				else if(rowOrCol == StsParameters.ROWCOL)
				{
				   int row = (int)gridPoint.iF;
				   int col = (int)gridPoint.jF;
				   distanceTransform.setDistance(row, col, 0.0f);
				}
				else
				{
				   int row = StsMath.floor(gridPoint.iF);
				   float dRowF = gridPoint.iF - row;
				   if(dRowF > 0.5f)
				   {
					   row++;
					   dRowF = 1.0f - dRowF;
				   }
				   int col = StsMath.floor(gridPoint.jF);
				   float dColF = gridPoint.jF - col;
				   if(dColF > 0.5f)
				   {
					   col++;
					   dColF = 1.0f - dColF;
				   }

				   float rowDistance = dRowF*rowGridSize;
				   float colDistance = dColF*colGridSize;
				   float distance = (float)Math.sqrt(rowDistance*rowDistance + colDistance*colDistance);
					/*
				   if(distance > stimulatedRadius/100.0f)
				   {
					   System.out.println("Setting distance to 0.0, outside radius.");
					   distance = largeFloat;
				   }
				   */
				   distanceTransform.setDistance(row, col, distance);
				}
			}
        }
        if(!hasPoints) return false;
        runDistanceTransform(distanceTransform, progressDialog);
        return true;
    }

    private void runDistanceTransform(StsDistanceTransform distanceTransform, StsProgressBarDialog progressDialog)
    {
//        StsModel model = StsModel.getCurrentModel();
//        model.disableDisplay();
        if(progressDialog == null)
            distanceTransform.distanceTransform(null);
        else
        {
            distanceTransform.distanceTransform(progressDialog.progressPanel);
            progressDialog.finished();
        }
//        fractureSetClass.initializeDistanceColorscale();
        constructDistanceBytePlaneData(distanceTransform.distances);
//        model.enableDisplay();
    }

    private void constructDistanceBytePlaneData(float[][] distances)
    {
        float maxGridSize = Math.max(fractureSetClass.getRowGridSize(), fractureSetClass.getColGridSize());
        float contourInterval = fractureSetClass.getContourInterval();
        if(contourInterval == 0.0f)
            contourInterval = getMinimumStimulatedRadius()/2;
        int n = 0;
        float contourLineHalfWidth = maxGridSize/2;
        for(int row = 0; row < nTextureRows; row++)
            for(int col = 0; col < nTextureCols; col++, n++)
            {
                float remainder = distances[row][col]%contourInterval;
                if(remainder < contourLineHalfWidth || (contourInterval - remainder) <= contourLineHalfWidth)
                    data[n] = lineColorByte;  
            }
    }

    private float getMinimumStimulatedRadius()
    {
        float radius = largeFloat;
        int nFractureSets = fractureSetClass.getSize();
        for(int s = 0; s < nFractureSets; s++)
        {
            StsFractureSet fractureSet = (StsFractureSet)fractureSetClass.getElement(s);
            if(!fractureSet.getIsVisible()) continue;
            radius = Math.min(radius, fractureSet.getStimulatedRadius());
        }
        return radius;
    }

    private boolean computeZBytePlaneData(StsRotatedGridBoundingBox fractureBoundingBox)
    {
        return computeZBytePlaneData(fractureBoundingBox, dirCoordinate);
    }

    public boolean computeZBytePlaneData(StsRotatedGridBoundingBox fractureBoundingBox, float z)
    {
        if(data == null) return false;
        Arrays.fill(data, (byte)0);
        StsFractureSet fractureSet = (StsFractureSet)fractureSetClass.getElement(0);
		byte fractureByteColorIndex = (byte)(fractureSet.getIndex() + nBasicColors);
        float stimulatedRadius = fractureSet.getStimulatedRadius();
		float rowGridSize = fracturesBoundingBox.yInc;
		float colGridSize = fractureBoundingBox.xInc;
		for(StsFractureDisplayable displayableFracture : displayableFractures)
		{
			StsPoint[] zLinePoints = displayableFracture.getSectionIntersectionAtZ(z);
			if(zLinePoints == null) continue;
			int nPoints = zLinePoints.length;
			float[][] centerPoints = new float[nPoints][2];
			for(int p = 0; p < nPoints; p++)
			{
				StsPoint point = zLinePoints[p];
				centerPoints[p][0] = fractureBoundingBox.getCursorCol(ZDIR, point.getX());
				centerPoints[p][1] = fractureBoundingBox.getCursorRow(ZDIR, point.getY());
			}
			constructFractureVolumePlane(centerPoints, fractureByteColorIndex, stimulatedRadius, rowGridSize, colGridSize);
		}
        return true;
    }

    private boolean computeXBytePlaneData(StsRotatedGridBoundingBox fractureBoundingBox)
    {
        int nFractureSets = fractureSetClass.getSize();
        for(int s = 0; s < nFractureSets; s++)
        {
            StsFractureSet fractureSet = (StsFractureSet)fractureSetClass.getElement(s);
            byte fractureByteIndex = (byte)(fractureSet.getIndex() + nBasicColors);
            float stimulatedRadius = fractureSet.getStimulatedRadius();
            float rowGridSize = fractureSet.getRowGridSize();
            float colGridSize = fractureSet.getColGridSize();
            StsObjectRefList fractureList = fractureSet.getFractureList();
            int nFractures = fractureList.getSize();
            float zMin = fractureBoundingBox.zMin;
            float zInc = fractureBoundingBox.zInc;
            int nSlices = nTextureCols;
            float boxColF = fractureBoundingBox.getColCoor(dirCoordinate);
            for(int f = 0; f < nFractures; f++)
            {
                StsFracture fracture = (StsFracture)fractureList.getElement(f);
                float z = zMin;
                for(int k = 0; k < nSlices; k++, z += zInc)
                {
                    StsPoint[] zLinePoints = fracture.getSectionIntersectionAtZ(z);
                    if(zLinePoints == null) continue;
                    int nPoints = zLinePoints.length;
                    float[][] centerPoints = new float[nPoints][2];
                    for(int p = 0; p < nPoints; p++)
                    {
                        StsPoint point = zLinePoints[p];
                        centerPoints[p][0] = fractureBoundingBox.getCursorCol(ZDIR, point.getX());
                        centerPoints[p][1] = fractureBoundingBox.getCursorRow(ZDIR, point.getY());
                    }
                    constructFractureVolumeLine(centerPoints, fractureByteIndex, k, boxColF, stimulatedRadius, rowGridSize, colGridSize);
                }
            }
        }
        return true;
    }

    private boolean computeYBytePlaneData(StsRotatedGridBoundingBox fractureBoundingBox)
    {
         int nFractureSets = fractureSetClass.getSize();
        for(int s = 0; s < nFractureSets; s++)
        {
            StsFractureSet fractureSet = (StsFractureSet)fractureSetClass.getElement(s);
            byte fractureByteIndex = (byte)(fractureSet.getIndex() + nBasicColors);
            float stimulatedRadius = fractureSet.getStimulatedRadius();
            float rowGridSize = fractureSet.getRowGridSize();
            float colGridSize = fractureSet.getColGridSize();
            StsObjectRefList fractureList = fractureSet.getFractureList();
            int nFractures = fractureList.getSize();
            float zMin = fractureBoundingBox.zMin;
            float zInc = fractureBoundingBox.zInc;
            int nSlices = nTextureCols;
            float boxRowF = fractureBoundingBox.getRowCoor(dirCoordinate);
            for(int f = 0; f < nFractures; f++)
            {
                StsFracture fracture = (StsFracture)fractureList.getElement(f);
                float z = zMin;
                for(int k = 0; k < nSlices; k++, z += zInc)
                {
                    StsPoint[] zLinePoints = fracture.getSectionIntersectionAtZ(z);
                    if(zLinePoints == null) continue;
                    int nPoints = zLinePoints.length;
                    float[][] centerPoints = new float[nPoints][2];
                    for(int p = 0; p < nPoints; p++)
                    {
                        StsPoint point = zLinePoints[p];
                        centerPoints[p][0] = fractureBoundingBox.getCursorCol(ZDIR, point.getX());
                        centerPoints[p][1] = fractureBoundingBox.getCursorRow(ZDIR, point.getY());
                    }
                    constructFractureVolumeLine(centerPoints, fractureByteIndex, k, boxRowF, stimulatedRadius, rowGridSize, colGridSize);
                }
            }
        }
        return true;
    }


    private void constructFractureVolumePlane(float[][] centerPoints, byte fractureByteIndex, float stimulatedRadius, float rowGridSize, float colGridSize)
    {
        int nPoints = centerPoints.length;        
        float[][] corners = new float[4][];
        for(int n = 0; n < nPoints-1; n++)
        {
            float[] normal = StsMath.horizontalNormal2D(centerPoints[n], centerPoints[n+1]);
            radiusScaleNormal(normal, stimulatedRadius, rowGridSize, colGridSize);
            corners[0] = StsMath.subtract(centerPoints[n], normal);
            corners[1] = StsMath.add(centerPoints[n], normal);
            corners[2] = StsMath.add(centerPoints[n+1], normal);
            corners[3] = StsMath.subtract(centerPoints[n+1], normal);
            BoundingBox boundingBox = new BoundingBox();
            for(int i = 0; i < 4; i++)
                boundingBox.addPoint(corners[i]);
            boundingBox.computeRange();
            for(int i = 0; i < 3; i++)
                boundingBox.addSegment(corners[i], corners[i+1]);
            boundingBox.addSegment(corners[3], corners[0]);
            boundingBox.fillDataPlane(fractureByteIndex, nTextureCols);
        }
        for(int n = 0; n < nPoints; n++)
            fillCircle(centerPoints[n], nTextureCols, fractureByteIndex, stimulatedRadius, rowGridSize, colGridSize);
    }


    private void constructFractureVolumeLine(float[][] centerPoints, byte fractureByteIndex, int nSlice, float boxRowColF, float stimulatedRadius, float rowGridSize, float colGridSize)
    {
        int nPoints = centerPoints.length;
        float[][] corners = new float[4][];
        for(int n = 0; n < nPoints-1; n++)
        {
            float[] normal = StsMath.horizontalNormal2D(centerPoints[n], centerPoints[n+1]);
            radiusScaleNormal(normal, stimulatedRadius, rowGridSize, colGridSize);
            corners[0] = StsMath.subtract(centerPoints[n], normal);
            corners[1] = StsMath.add(centerPoints[n], normal);
            corners[2] = StsMath.add(centerPoints[n+1], normal);
            corners[3] = StsMath.subtract(centerPoints[n+1], normal);

            BoundingBox boundingBox = new BoundingBox();
            boundingBox.initializeLine();
            for(int i = 0; i < 3; i++)
                boundingBox.addIntersection(corners[i], corners[i+1], boxRowColF);
            boundingBox.addIntersection(corners[3], corners[0], boxRowColF);
            boundingBox.fillDataPlaneLine(fractureByteIndex, nTextureCols, nSlice);
        }
        for(int n = 0; n < nPoints; n++)
            circleLine(centerPoints[n], nTextureCols, fractureByteIndex, boxRowColF, nSlice, stimulatedRadius, rowGridSize, colGridSize);
    }

    private void radiusScaleNormal(float[] normal, float radius, float rowGridSize, float colGridSize)
    {
        normal[0] *= radius/colGridSize;
        normal[1] *= radius/rowGridSize;
    }

    class BoundingBox
    {
        float rowMinF = largeFloat;
        float rowMaxF = -largeFloat;
        float colMinF = largeFloat;
        float colMaxF = -largeFloat;
        int rowMin, rowMax, colMin, colMax;
        int nRows;
        float[][] rowColPoints = null;
        float[] colRowPoints;

        BoundingBox()
        {
        }

        void addPoint(float[] colRowF)
        {
            rowMinF = Math.min(rowMinF, colRowF[1]);
            rowMaxF = Math.max(rowMaxF, colRowF[1]);
            colMinF = Math.min(colMinF, colRowF[0]);
            colMaxF = Math.max(colMaxF, colRowF[0]);
        }

        void computeRange()
        {
            rowMin = StsMath.floor(rowMinF);
            rowMax = StsMath.ceiling(rowMaxF);
            colMin = StsMath.floor(colMinF);
            colMax = StsMath.ceiling(colMaxF);
            nRows = rowMax - rowMin + 1;
            rowColPoints = new float[nRows][0];
        }

        void initializeLine()
        {
            colRowPoints = null;
        }

        void addSegment(float[] point0, float[] point1)
        {
            float colF0, colF1;

            float rowF0 = point0[1];
            float rowF1 = point1[1];
            if(rowF0 == rowF1)
                return;
            else if(rowF0 < rowF1)
            {
                colF0 = point0[0];
                colF1 = point1[0];
            }
            else
            {
                rowF0 = point1[1];
                rowF1 = point0[1];
                colF0 = point1[0];
                colF1 = point0[0];
            }
            int rowStart = StsMath.ceiling(rowF0);
            int rowEnd = StsMath.floor(rowF1);
            float df = 1/(rowF1 - rowF0);
            float f = (rowStart - rowF0)*df;
            float colF = colF0 + f*(colF1 - colF0);
            float dColF = df*(colF1 - colF0);
            for(int row = rowStart; row <= rowEnd; row++, colF += dColF)
                addRowColPoint(row, colF);
        }

        void addIntersection(float[] point0, float[] point1, float planeRowColF)
        {
            float planeRowF;
            if(dirNo == XDIR)
            {
                planeRowF = StsMath.interpolateValueNoExtrapolation(planeRowColF, point0, point1, 0, 1);
            }
            else if(dirNo == YDIR)
            {
                planeRowF = StsMath.interpolateValueNoExtrapolation(planeRowColF, point0, point1, 1, 0);
            }
            else
                return;
            if(planeRowF == StsParameters.nullValue) return;
            addColRowPoint(planeRowF);
        }

        /** inserts this col number in a sorted list of column numbers for this row */
        void addRowColPoint(int row, float colF)
        {
            rowColPoints[row - rowMin] = StsMath.floatListAddSortedValue(rowColPoints[row - rowMin], colF);
        }

        void addColRowPoint(float rowF)
        {
            colRowPoints = StsMath.floatListAddSortedValue(colRowPoints, rowF);
        }

        void fillDataPlane(byte fractureByteIndex, int nDataCols)
        {
            for(int row = rowMin; row <= rowMax; row++)
            {
                float[] colFs = rowColPoints[row - rowMin];
                int nPoints = colFs.length;
                if(nPoints > 0)
                {
                    if(nPoints%2 != 0)
                    {
                        nPoints--;
                        // StsException.systemError(this, "fillDataPlane", "Not even number of points");
                        // continue;
                    }
                    for(int n = 0; n < nPoints; n += 2)
                        fillDataRow(row, colFs[n], colFs[n+1], nDataCols, fractureByteIndex);
                }
            }
        }


        void fillDataPlaneLine(byte fractureByteIndex, int nDataCols, int col)
        {
            if(colRowPoints == null) return;
            int nPoints = colRowPoints.length;
            if(nPoints > 0)
            {
                if(nPoints%2 != 0)
                {
                    StsException.systemError(this, "fillDataPlane", "Not even number of points");
                    return;
                }
                for(int n = 0; n < nPoints; n += 2)
                    fillDataCol(col, colRowPoints[n], colRowPoints[n+1], nDataCols, fractureByteIndex);
            }
        }
    }

    private void fillDataRow(int row, float colStartF, float colEndF, int nDataCols, byte fillValue)
    {
        int colStart = StsMath.ceiling(colStartF);
        if(colStart >= nTextureCols) return;
        colStart = Math.max(colStart, 0);
        int colEnd = StsMath.floor(colEndF);
        colEnd = Math.min(colEnd, nTextureCols -1);
        int iStart = row*nDataCols + colStart;
        if(iStart < 0) return;
        int iEnd = row*nDataCols + colEnd;
        if(iEnd < iStart) return;
        Arrays.fill(data, iStart, iEnd+1, fillValue);
    }

    private void fillDataCol(int col, float planeStartRowF, float planeRowEndF, int nDataCols, byte fillValue)
    {
        int rowStart = StsMath.ceiling(planeStartRowF);
        if(rowStart >= nTextureRows) return;
        rowStart = Math.max(rowStart, 0);
        int rowEnd = StsMath.floor(planeRowEndF);
        rowEnd = Math.min(rowEnd, nTextureRows -1);
        int i = rowStart*nDataCols + col;
        for(int row = rowStart; row <= rowEnd; row++, i += nDataCols)
            data[i] = fillValue;
    }

    private void fillCircle(float[] centerPoint, int nDataCols, byte fillByteIndex, float radius, float rowGridSize, float colGridSize)
    {
        float colMinF, colMaxF;
        float rowRadius = radius/rowGridSize;
        int rowMin = StsMath.floor(centerPoint[1] - rowRadius);
        rowMin = Math.max(rowMin, 0);
        int rowMax = StsMath.ceiling(centerPoint[1] + rowRadius);
        rowMax = Math.min(rowMax, nTextureRows -1);
        float radiusSq = radius*radius;
        for(int row = rowMin; row <= rowMax; row++)
        {
            float dRow = rowGridSize*(row - centerPoint[1]);
            if(radiusSq < dRow*dRow) continue;
            float dColF = (float)Math.sqrt(radiusSq - dRow*dRow)/colGridSize;
            colMinF = centerPoint[0] - dColF;
            colMaxF = centerPoint[0] + dColF;
            fillDataRow(row, colMinF, colMaxF, nDataCols, fillByteIndex);
        }
    }

    private void circleLine(float[] centerPoint, int nDataCols, byte fillByteIndex, float boxRowColF, int planeCol, float radius, float rowGridSize, float colGridSize)
    {
        if(dirNo == XDIR)
        {
            float radiusSq = radius*radius;
            float dColF = colGridSize*(boxRowColF - centerPoint[0]);
            float dColFSq = dColF*dColF;
            if(radiusSq < dColF*dColF) return;
            float dRowF = (float)Math.sqrt(radiusSq - dColFSq)/rowGridSize;
            float boxRowMinF = centerPoint[1] - dRowF;
            float boxRowMaxF = centerPoint[1] + dRowF;
            fillDataCol(planeCol, boxRowMinF, boxRowMaxF, nDataCols, fillByteIndex);
        }
        else if(dirNo == YDIR)
        {
            float radiusSq = radius*radius;
            float dRowF = rowGridSize*(boxRowColF - centerPoint[1]);
            float dRowFSq = dRowF*dRowF;
            if(radiusSq < dRowF*dRowF) return;
            float dColF = (float)Math.sqrt(radiusSq - dRowFSq)/colGridSize;
            float boxColMinF = centerPoint[0] - dColF;
            float boxColMaxF = centerPoint[0] + dColF;
            fillDataCol(planeCol, boxColMinF, boxColMaxF, nDataCols, fillByteIndex);
        }
    }

    public boolean setDirCoordinate(float dirCoordinate)
    {
        if(this.dirCoordinate == dirCoordinate) return false;
        this.dirCoordinate = dirCoordinate;
        textureChanged();
        if(textureTiles != null) textureTiles.setTilesDirCoordinate(dirCoordinate);
        return true;
    }
    public boolean isDisplayableObject(Object object) { return object instanceof StsFractureSet; }

    public boolean isDisplayingObject(Object object)
    {
        return true;
    }

    public boolean setObject(Object object)
	{
       if(object == null) isVisible = false;
       textureChanged();
       return true;
    }

    public Object getObject() { return null; }

    public Class getDisplayableClass() { return StsFracture.class; }

    public boolean canDisplayClass(Class c) { return c == StsFractureSet.class; }

    private boolean getDisplayedOnOtherCursorSections()
    {
        StsCursorSection cursorSection = cursor3d.cursorSections[dirNo];
        StsCursor3dTexture[] visibleCursorSections = cursorSection.getVisibleDisplayableSections();
        return visibleCursorSections.length > 1;
    }

    private byte[] applySubVolume()
    {
        if(data == null) return null;
        if(!fractureSetClass.getDisplayOnSubVolumes()) return data;
        StsSubVolumeClass subVolumeClass = (StsSubVolumeClass)StsModel.getCurrentModel().getStsClass(StsSubVolume.class);
        if(subVolumeClass == null || subVolumeClass.getSubVolumes().length == 0) return data;

        if(subVolumeData == null)
        {
            if(debug) System.out.println("recomputing subVolume for dir: " + dirNo + " dirCoor: " + dirCoordinate);
            if(!computeSubVolumePlane()) return data;
        }

        if(debug) System.out.println("applying subVolume for dir: " + dirNo + " dirCoor: " + dirCoordinate);
        for(int n = 0; n < data.length; n++)
        {
            if(subVolumeData[n] != 0)
                subVolumeData[n] = data[n];
        }
        return subVolumeData;
    }

    public void drawTextureTileSurface(StsTextureTile tile, GL gl, boolean is3d)
    {
		StsClassFractureDisplayable[] fractureCursorSectionClasses = getClassFractureDisplayables();
		int nCursorSectionClasses = fractureCursorSectionClasses.length;
		for(int n = 0; n < nCursorSectionClasses; n++)
		{
			StsRotatedGridBoundingBox fractureBoundingBox = fractureCursorSectionClasses[n].getFractureBoundingBox();
			drawTextureTileSurface(fractureBoundingBox, tile, gl, is3d, zDomain);
		}
    }

	public StsClassFractureDisplayable[] getClassFractureDisplayables()
	{
		StsArrayList fractureCursorSectionClassList = model.getFractureDisplayableClasses();
		int nCursorSectionClasses = fractureCursorSectionClassList.size();
		StsClassFractureDisplayable[] fractureCursorSectionClasses = new StsClassFractureDisplayable[nCursorSectionClasses];
		for(int n = 0; n < nCursorSectionClasses; n++)
			fractureCursorSectionClasses[n] = (StsClassFractureDisplayable)fractureCursorSectionClassList.get(n);
		return fractureCursorSectionClasses;
	}
/*
    public void drawTextureTileSurface(StsTextureTile tile, GL gl, boolean is3d, int n)
    {
		if(is3d)
			tile.drawQuadSurface3d(gl);
		else
			tile.drawQuadSurface2d(gl);
    }
*/
/*
	public void drawTextureTileSurface2d(StsTextureTile tile, GL gl)
	{
		tile.drawQuadSurface2d(gl);
	}
*/
    // this could be put in a displayList, but is so small in terms of the number of colors
    // that it is easy to recompute on the fly
 
    public void display(StsGLPanel3d glPanel3d, boolean is3d)
    {
    }

    public boolean hasData(int nPlane) { return nPlane == this.nPlane; }


    public void cropChanged()
    {
        textureTiles = null;
        textureChanged();
    }

    public void subVolumeChanged()
    {
        if(debug) StsException.systemDebug(this, "subVolumeChanged", "subVolumePlane cleared for dir: " + dirNo + " dirCoor: " + dirCoordinate);
        subVolumeData = null;
    }

    private boolean computeSubVolumePlane()
    {
        StsModel model = StsModel.getCurrentModel();
        StsSubVolumeClass subVolumeClass = (StsSubVolumeClass)model.getStsClass(StsSubVolume.class);
        if(subVolumeClass == null || subVolumeClass.getSubVolumes().length == 0)return false;

		StsClassFractureDisplayable[] fractureClasses = getClassFractureDisplayables();
		StsRotatedGridBoundingBox fracturesBoundingBox = new StsRotatedGridBoundingBox();
		for(StsClassFractureDisplayable fractureClass : fractureClasses)
		{
        	StsRotatedGridBoundingBox fractureBoundingBox = fractureClass.getFractureBoundingBox();
			fracturesBoundingBox.addBoundingBox(fractureBoundingBox);
		}
        byte zDomainProject = StsModel.getCurrentZDomain();
        subVolumeData = subVolumeClass.getSubVolumePlane(dirNo, dirCoordinate, fracturesBoundingBox, zDomainProject);
        return subVolumeData != null;
    }

    public String propertyReadout(StsPoint point)
    {
        StsRotatedGridBoundingBox boundingBox = model.getProject().getRotatedBoundingBox();
        if(boundingBox == null) return "";

        if(dirNo == XDIR)
        {
            int row = Math.round(boundingBox.getRowCoor(point.v[1]));
            int col = Math.round(boundingBox.getSliceCoor(point.v[2]));
            return getTypeString(row, col);
        }
        else if(dirNo == YDIR)
        {
            int row = Math.round(boundingBox.getColCoor(point.v[0]));
            int col = Math.round(boundingBox.getSliceCoor(point.v[2]));
            return getTypeString(row, col);
        }
        else
        {
            int row = Math.round(boundingBox.getRowCoor(point.v[1]));
            int col = Math.round(boundingBox.getColCoor(point.v[0]));
            return getTypeString(row, col);
        }
    }

    private String getTypeString(int row, int col)
    {
        if(data == null || !isInRange(row, col)) return " ";
        int nFracture = data[row* nTextureCols + col] - 1;
        StsFracture fracture = (StsFracture)fractureSetClass.getElement(nFracture);
        return fracture.getName();
    }     

    public String getName()
    {
        return "Cursor view[" + dirNo + "] of: " + "fractureSetClass";
    }
}