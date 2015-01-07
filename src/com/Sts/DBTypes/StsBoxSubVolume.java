package com.Sts.DBTypes;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;

public class StsBoxSubVolume extends StsRotatedGridBoundingBox
{
    StsGridPoint boxCenter = null;
    StsSeismicVolume volume;
	/** Display lists should be used (controlled by View:Display Options) */
	boolean useDisplayLists;
	/** Display lists currently being used for surface geometry */
	boolean usingDisplayLists = true;	/** indicates whether cursorSection is being drawn in time or depth. */
	boolean textureChanged = true;

    public StsBoxSubVolume()
    {
    }

    public StsBoxSubVolume(boolean persistent)
    {
        super(persistent);
    }

    public StsBoxSubVolume(StsCursorPoint cursorPoint, StsSeismicVolume volume, int nRows, int nCols, int nSlices)
    {
        StsToolkit.copySubToSuperclass(volume, this, StsRotatedGridBoundingBox.class, StsBoundingBox.class, false);
        this.volume = volume;
        this.nRows = nRows;
        this.nCols = nCols;
        this.nSlices = nSlices;
        float[] xyz = cursorPoint.point.v;
        boxCenter = new StsGridPoint(cursorPoint.point, volume);

        int centerRow = volume.getNearestBoundedRowCoor(xyz[1]);
        int rowMin = centerRow - nRows/2;
        int rowMax = centerRow + nRows/2;
        nRows = rowMax - rowMin + 1;
        yMin = volume.getYCoor(rowMin);
        yMax = volume.getYCoor(rowMax);
        rowNumMin = volume.getNearestBoundedRowNumFromY(yMin);
        rowNumMax = volume.getNearestBoundedRowNumFromY(yMax);

        int centerCol = volume.getNearestBoundedColCoor(xyz[0]);
        int colMin = centerCol - nCols/2;
        int colMax = centerCol + nCols/2;
        nCols = colMax - colMin + 1;
        xMin = volume.getXCoor(colMin);
        xMax = volume.getXCoor(colMax);
        colNumMin = volume.getNearestBoundedColNumFromX(xMin);
        colNumMax = volume.getNearestBoundedColNumFromX(xMax);

        int centerSlice = volume.getNearestBoundedSliceCoor(xyz[2]);
        int sliceMin = centerSlice - nSlices/2;
        int sliceMax = centerSlice + nSlices/2;
        nSlices = sliceMax - sliceMin + 1;
        zMin = volume.getZCoor(sliceMin);
        zMax = volume.getZCoor(sliceMax);
    }

    public StsBoxSubVolume(StsSeismicVolume volume, float rowNumMin, float rowNumMax, float colNumMin, float colNumMax, float sliceMin, float sliceMax)
    {
        StsToolkit.copySubToSuperclass(volume, this, StsRotatedGridBoundingBox.class, StsBoundingBox.class, false);
        this.volume = volume;
        volume.setStsDirectory(currentModel.getProject().getRootDirString());

        setRowNumMin(rowNumMin);
        setRowNumMax(rowNumMax);
        setColNumMin(colNumMin);
        setColNumMax(colNumMax);
        setSliceMin(sliceMin);
        setSliceMax(sliceMax);

        zMin = volume.getZCoor(sliceMin);
        zMax = volume.getZCoor(sliceMax);
        float[] xyz = new float[3];
        xyz[0] = (xMax + xMin)/2;
        xyz[1] = (yMax + yMin)/2;
        xyz[2] = (zMax + zMin)/2;

        StsCursorPoint cursorPoint = new StsCursorPoint(0, new StsPoint(xyz));
        boxCenter = new StsGridPoint(cursorPoint.point, volume);
    }

	public boolean initialize(StsModel model)
	{
		isVisible = false;
		return true;
	}

    public void setXMin(float xMin)
    {
        xMin = Math.max(xMin, volume.xMin);
        this.xMin = xMin;
        colNumMin = volume.getNearestBoundedColNumFromX(xMin);
        if(boxCenter != null) resetCenter();
        dbFieldChanged("xMin", xMin);
    }

    public void setXMax(float xMax)
    {
        xMax = Math.min(xMax, volume.xMax);
        this.xMax = xMax;
        colNumMax = volume.getNearestBoundedColNumFromX(xMax);
        if(boxCenter != null) resetCenter();
        dbFieldChanged("xMax", xMax);
    }

    public void setYMin(float yMin)
    {
        yMin = Math.max(yMin, volume.yMin);
        this.yMin = yMin;
        rowNumMin = volume.getNearestBoundedRowNumFromY(yMin);
        if(boxCenter != null) resetCenter();
        dbFieldChanged("yMin", yMin);
    }

    public void setYMax(float yMax)
    {
        yMax = Math.min(yMax, volume.yMax);
        this.yMax = yMax;
        rowNumMax = volume.getNearestBoundedRowNumFromY(yMax);
        if(boxCenter != null) resetCenter();
        dbFieldChanged("yMax", yMax);
    }

    public void setRowNumMin(float rowNumMin)
    {
        yMin = volume.getYFromRowNum(rowNumMin);
        this.rowNumMin = rowNumMin;
        if(boxCenter != null) resetCenter();
        dbFieldChanged("yMin", yMin);
    }

    public void setRowNumMax(float rowNumMax)
    {
        yMax = volume.getYFromRowNum(rowNumMax);
        this.rowNumMax = rowNumMax;
        if(boxCenter != null) resetCenter();
        dbFieldChanged("yMax", yMax);
    }

    public void setColNumMin(float colNumMin)
    {
        xMin = volume.getXFromColNum(colNumMin);
        this.colNumMin = colNumMin;
        if(boxCenter != null) resetCenter();
        dbFieldChanged("xMin", xMin);
    }

    public void setColNumMax(float colNumMax)
    {
        xMax = volume.getXFromColNum(colNumMax);
        this.colNumMax = colNumMax;
        if(boxCenter != null) resetCenter();
        dbFieldChanged("xMax", xMax);
    }

    public void setSliceMin(float sliceMin)
    {
        zMin = volume.getZCoor(sliceMin);
        if(boxCenter != null) resetCenter();
        dbFieldChanged("zMin", zMin);
    }

    public void setSliceMax(float sliceMax)
    {
        zMax = volume.getZCoor(sliceMax);
        if(boxCenter != null) resetCenter();
        dbFieldChanged("zMax", zMax);
    }

    public int getSliceMin()
    {
        return volume.getNearestSliceCoor(zMin);
    }

    public int getSliceMax()
    {
        return volume.getNearestSliceCoor(zMax);
    }

	public void adjustXMin(float xMin)
	{
		setXMin(xMin);
		currentModel.addMethodCmd(this, "setXMin", new Object[] { new Float(xMin) }, "boxSubVolume.setXMin");
	}

	public void adjustXMax(float xMax)
	{
		setXMax(xMax);
		currentModel.addMethodCmd(this, "setXMax", new Object[] { new Float(xMax) }, "boxSubVolume.setXMax");
	}

	public void adjustYMin(float yMin)
	{
		setYMin(yMin);
		currentModel.addMethodCmd(this, "setYMin", new Object[] { new Float(yMin) }, "boxSubVolume.setYMin");
	}

	public void adjustYMax(float yMax)
	{
		setYMax(yMax);
		currentModel.addMethodCmd(this, "setYMax", new Object[] { new Float(yMax) }, "boxSubVolume.setYMax");
	}

	public void adjustZMin(float zMin)
	{
		setZMin(zMin);
		currentModel.addMethodCmd(this, "setZMin", new Object[] { new Float(zMin) }, "boxSubVolume.setZMin");
	}

	public void adjustZMax(float zMax)
	{
		setZMin(zMax);
		currentModel.addMethodCmd(this, "setZMax", new Object[] { new Float(zMax) }, "boxSubVolume.setZMax");
	}

    public StsSeismicVolume getVolume()
    {
        return volume;
    }
    public void setVolume(StsSeismicVolume volume)
    {
        this.volume = volume;
    }

    public int getRowCenter() { return boxCenter.row; }
    public int getColCenter() { return boxCenter.col; }
    public int getSliceCenter() { return volume.getNearestSliceCoor(boxCenter.getZorT()); }

    public void display(StsGLPanel glPanel, StsColor stsColor, byte action, boolean editing)
    {
		GL gl = glPanel.getGL();
        displayBoundingBox(gl, stsColor, StsGraphicParameters.gridLineWidth);

        if (boxCenter == null)
        {
            return;
        }

        int boxSize = 4;
        if (editing && action == StsBoxSetSubVolume.ACTION_MOVE_BOX)
            boxSize = 8;

        float[] xyz = boxCenter.getXYZorT(currentModel.isDepth);
        StsGLDraw.drawPoint(xyz, StsColor.BLACK, glPanel, boxSize + 4, 1.0);
        StsGLDraw.drawPoint(xyz, stsColor, glPanel, boxSize, 2.0);
        displayBoundingBox(gl, stsColor, StsGraphicParameters.gridLineWidth);
        if (editing && action == StsBoxSetSubVolume.ACTION_MOVE_POINT)
            drawFacePoints(gl, false);
    }

    public void pickCenterPoint(GL gl, boolean editing, int boxID)
    {
        if (boxCenter == null)
        {
            return;
        }

        int boxSize = 4;
        if (editing) boxSize = 8;

        gl.glInitNames();
        gl.glPushName(boxID);
        gl.glPushName( -1);
        StsGLDraw.drawPoint(boxCenter.getXYZorT(), gl, boxSize + 4);
        gl.glPopName();
        gl.glPopName();
    }

    public void pickFacePoint(GL gl, int boxID)
    {
        gl.glInitNames();
        gl.glPushName(boxID);
        drawFacePoints(gl, true);
    }

    public void moveBox(StsMouse mouse, StsGLPanel3d glPanel3d)
    {
        StsPoint newCenterPoint;

        StsView currentView = glPanel3d.getView();
        if (currentView instanceof StsView3d)
        {
            int planeDir;

            double[] lineVector = glPanel3d.getViewVectorAtMouse(mouse);
            if (Math.abs(lineVector[0]) > Math.abs(lineVector[1]))
            {
                planeDir = 0;
                if (Math.abs(lineVector[2]) >= Math.abs(lineVector[0]))
                    planeDir = 2;
            }
            else
            {
                planeDir = 1;
                if (Math.abs(lineVector[2]) >= Math.abs(lineVector[1]))
                    planeDir = 2;
            }
            StsCursor3d cursor3d = currentModel.win3d.getCursor3d();
            newCenterPoint = cursor3d.getPointInPlaneAtMouse(glPanel3d, planeDir, boxCenter.getPoint(), mouse);
            if (newCenterPoint == null)
            {
                return;
            }
        }
        else if (currentView instanceof StsViewCursor)
        {
            StsCursorPoint cubePoint = ( (StsViewCursor) currentView).getCursorPoint(mouse);
            if (cubePoint == null)
            {
                return;
            }
            newCenterPoint = cubePoint.point;
        }
        else
        {
            return;
        }

        StsPoint dCenter = StsPoint.subPointsStatic(newCenterPoint, boxCenter.getPoint());
        adjustXYZPosition(dCenter);
        boxCenter = new StsGridPoint(dCenter, volume);
    }

    public void movePoint(StsMouse mouse, StsGLPanel3d glPanel3d, int pickedPointIndex)
    {
        float[] faceCenter;
        StsPoint newFaceCenterPoint;

        StsView currentView = glPanel3d.getView();
        if (! (currentView instanceof StsView3d))
        {
            return;
        }

        // X and Y face points will be moved in the Z-plane thru the face point
        // Z face points will be moved in the X or Y planes depending on orientation
        int planeDir = 2;
        if (pickedPointIndex > 3)
        {
            double[] lineVector = glPanel3d.getViewVectorAtMouse(mouse);
            if (Math.abs(lineVector[0]) > Math.abs(lineVector[1]))
            {
                planeDir = 0;
            }
            else
            {
                planeDir = 1;
            }
        }
        int moveDir = pickedPointIndex / 2;
        faceCenter = getFaceCenter(pickedPointIndex);
        if(faceCenter == null) return;
        StsCursor3d cursor3d = currentModel.win3d.getCursor3d();
        newFaceCenterPoint = cursor3d.getPointInPlaneAtMouse(glPanel3d, planeDir, faceCenter, mouse);
        if (newFaceCenterPoint == null) return;
        float[] dCenter = StsMath.subtract(newFaceCenterPoint.v, faceCenter);
        adjustRange(pickedPointIndex, dCenter[moveDir]);
        resetCenter();
    }

    private void resetCenter()
    {
        boxCenter.setX((xMin + xMax)/2);
        boxCenter.setY((yMin + yMax)/2);
        boxCenter.setZorT(currentModel.isDepth,(zMin + zMax)/2);
    }

    public StsGridBoundingBox getGridBoundingBox()
    {
        StsRotatedGridBoundingBox cursor3dBoundingBox = currentModel.getProject().getRotatedBoundingBox();
        StsGridBoundingBox boundingBox = new StsGridBoundingBox(false);
        boundingBox.initialize(cursor3dBoundingBox, this);
        return boundingBox;
    }

	public void displayVoxels(StsGLPanel glPanel)
	{
/*
		GL gl = glPanel.getGL();

		if (seismicVolumeChanged())
		  {
			  textureChanged = true;
			  if (texture != 0)
			  {
				  deleteTexture(gl);
			  }
		  }

		if((textureChanged || texture == 0) )
		{
			initializeTexture();
		}

		useDisplayLists = StsObject.getCurrentModel().getBooleanProperty("Use Display Lists");

		if(textureChanged || useDisplayLists && !usingDisplayLists)
		{
			if(textureChanged)       textureChanged = false;
			else if(useDisplayLists) usingDisplayLists = true;
			deleteDisplayLists(gl);
			usingDisplayLists = false;
		}
		if (nPlane == -1)
		{
			return;
		}
		if (textureTiles.isDirCoordinateCropped())
		{
			return;
		}

		gl.glDisable(GL.GL_LIGHTING);
		gl.glEnable(GL.GL_TEXTURE_2D);
		gl.glEnable(GL.GL_BLEND);
		gl.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE);
		gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
		gl.glShadeModel(GL.GL_FLAT);

		int colorListNum = line2d.getColorDisplayListNum(gl, false);
		if(mainDebug) System.out.println("Using colorListNum " + colorListNum + " for seismicCursorSection in dirNo " + dirNo);
		gl.glCallList(colorListNum);

		byte[] isDisplayData = null;
		if (textureChanged || useDisplayLists && !usingDisplayLists)
		{
			planeData = (byte[]) line2d.readPlaneData(dirNo, dirCoordinate);
			if(planeData != null)
			{
				textureChanged = false;
				if (subVolumePlane != null)
				{
					planeData = applySubVolume(planeData, subVolumePlane);
				}
			}
			else
			{
				planeData = transparentData();
				if (subVolumePlane != null)
				{
					planeData = applySubVolume(planeData, subVolumePlane);
				}

			}
			isDisplayData = planeData;
		}

		if(useDisplayLists && !usingDisplayLists)
		 {
			 if(textureTiles == null) StsException.systemError("StsSurface.displaySurfaceFill() failed. textureTiles should not be null.");
			 textureTiles.constructSurface(this, gl, useDisplayLists);
			 usingDisplayLists = true;
	 }
		if (is3d)
		{
			textureTiles.displayTiles(this, gl, isPixelMode, isDisplayData);
			gl.glDisable(GL.GL_TEXTURE_2D);
			gl.glDisable(GL.GL_BLEND);
			gl.glEnable(GL.GL_LIGHTING);
		}
		else
		{
			StsViewCursor viewCursor = (StsViewCursor) glPanel3d.getCurrentView();
			textureTiles.displayTiles2d(gl, isPixelMode, viewCursor.axesFlipped, planeData);
			gl.glDisable(GL.GL_TEXTURE_2D);
			gl.glDisable(GL.GL_BLEND);
			gl.glEnable(GL.GL_LIGHTING);
		}
*/
    }
/*
	private boolean seismicVolumeChanged()
	{
		StsSeismicVolumeClass lineSet = (StsSeismicVolumeClass) glPanel3d.model.getStsClass(line2d.getClass());
		StsSeismicVolume currentSeismicVolume = lineSet.getCurrentSeismicVolume();
		if(volume == currentSeismicVolume) return false;
		if(currentSeismicVolume == null) return true;
		return true;
    }
*/
}
