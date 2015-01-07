package com.Sts.DBTypes;

/**
 * <p>Title: PreStack3d Data Loading</p>
 * <p>Copyright: Copyright (c) 2006</p>
 * <p>Company: S2SSystems LLC</p>
 * @author T.J.Lasseter
 * @version c75k
 */

import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.PreStack.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.io.*;
import java.nio.*;
import java.nio.channels.*;


/** class instance for each 2d prestack line which is loaded */
public class StsPreStackLine2d extends StsPreStackLine implements StsTextureSurfaceFace
{
    transient boolean planeOK = false;
    transient boolean lighting = true;
	transient float tMin, tMax, tInc;
	transient boolean textureChanged = true;
	transient boolean deleteTexture = false;
	transient StsTextureTiles textureTiles = null;
	/** index of current gather on this line */
	transient int currentIndex = 0;
	/** Display lists should be used (controlled by View:Display Options) */
	transient boolean useDisplayLists;
	/** Display lists currently being used for surface geometry */
	transient boolean usingDisplayLists = true;
	/** indicates velocity is being isVisible on this line; otherwise stack is isVisible */
//	transient boolean displayVelocity = true;
    /** current stacking option for this line */
    transient byte stackOption = StsPreStackLineSetClass.STACK_LINES;

    transient boolean lineStackOk = false;
    transient boolean drawBlank = true;

    transient ByteBuffer velocityByteBuffer = null;

    transient static public final boolean debugTimer = false;
    transient static final boolean debug = false;
	public final static byte nullByte = StsParameters.nullByte;
	final static double roundOff = StsParameters.roundOff;
	
    public float analysisColStart;
    public int analysisColInc = 25;

    public StsPreStackLine2d()
	{
	}

	public StsPreStackLine2d(boolean persistent)
	{
		super(persistent);
	}

	/*
	 public StsPreStackLine2d(StsFile file, StsModel model) throws FileNotFoundException
	 {
	  this(file.getDirectory(), file.createFilename(), model);
	 }
	 */
	private StsPreStackLine2d(String directory, String filename, int lineIndex, StsModel model, StsPreStackLineSet2d lineSet) throws FileNotFoundException, StsException
	{
        super(directory, filename, lineIndex, model, lineSet);
		initializeGatherIndexes();
	}

	static public StsPreStackLine2d constructor(StsFile file, int lineIndex, StsModel model, StsPreStackLineSet2d lineSet)
	{
		return constructor(file.getDirectory(), file.getFilename(), lineIndex, model, lineSet);
	}

	static public StsPreStackLine2d constructor(String directory, String filename, int lineIndex, StsModel model, StsPreStackLineSet2d lineSet)
	{
		try
		{
			return new StsPreStackLine2d(directory, filename, lineIndex, model, lineSet);
		}
		catch(FileNotFoundException e)
		{
			return null;
		}
		catch(Exception e)
		{
			StsException.outputException("StsPreStackLine2d.constructor() failed.", e, StsException.WARNING);
			return null;
		}
	}

    public boolean initialize(StsModel model)
    {
        super.initialize(model);
        computeNormals();
        return true;
    }

    protected void initializeGatherIndexes()
	{
		minGatherIndex = 0;
		maxGatherIndex = nCols - 1;
	}
/*
	public boolean classInitialize(StsModel model)
	{
        return super.classInitialize(model);
	}
*/
    public void setPlaneOK(boolean ok)
    {
        planeOK = ok;
    }

    public int getDefaultShader() { return StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS; }
    public boolean getUseShader() { return lineSet.lineSetClass.getContourColors(); }

	public double[] getAttributeArray(String attributeName, int index) throws FileNotFoundException
	{
		int start = 0;
        int num = 0;

		StsMappedDoubleBuffer doubleBuffer = getAttributeArrayBuffer(attributeName, nLineTraces);
		if(doubleBuffer == null) return null;

		start = getNFirstGatherTrace(index);
		num = getNTracesInGather(index);
		if(num == 0)
			return null;

		double[] newArray = new double[num];
		doubleBuffer.position(start);
		doubleBuffer.get(newArray);
		doubleBuffer.close();
		return newArray;
	}

	public long getGatherFileOffset(int row, int col)
	{
	    return getGatherFileOffset(col);
	}

	public int getGatherFileLength(int row, int col)
	{
		return getGatherFileLength(col);
	}

    public String getGatherDescription(int row, int col)
    {
        return new String(stemname + " line " + row + " gather " + cdp[col]);
    }

    public boolean setLineXYs()
	{
        if(!setLineXYsFromCDPs())
            return setLineXYsFromShtRecs();
        else
            return true;
    }

    public boolean setLineXYsFromShtRecs()
    {
        try
		{
			double[] lineShtX = getAttributeArray(StsSEGYFormat.SHT_X);
			double[] lineShtY = getAttributeArray(StsSEGYFormat.SHT_Y);
			double[] lineRcvX = getAttributeArray(StsSEGYFormat.REC_X);
			double[] lineRcvY = getAttributeArray(StsSEGYFormat.REC_Y);
            double[] lineCdp = getAttributeArray(StsSEGYFormat.CDP);
			cdpX = new float[nCols];
			cdpY = new float[nCols];
			cdp = new int[nCols];
            int nLastTrace = -1;
            for(int n = 0; n < nCols; n++)
			{
                int nFirstTrace = nLastTrace + 1;
                nLastTrace = nLastTraceInGathers[n];
                double x = 0;
                double y = 0;
                for(int t = nFirstTrace; t <= nLastTrace; t++)
                {
                    x += lineShtX[t] + lineRcvX[t];
                    y += lineShtY[t] + lineRcvY[t];
                }
                int nTraces = nLastTrace - nFirstTrace + 1;
                x = x/nTraces/2;
                y = y/nTraces/2;
                float[] xy = currentModel.getProject().getRotatedRelativeXYFromUnrotatedAbsoluteXY(x, y);
				cdpX[n] = xy[0];
				cdpY[n] = xy[1];
				cdp[n] = (int)lineCdp[nLastTrace];

                //  Must reset origin since it is likely wrong because the cdp XY were wrong.
                if(n == 0)
                {
                    this.originSet = false;
                    this.checkSetOriginAndAngle(cdpX[n], cdpY[n], angle);
                }
			}
            computeNormals();
            return true;
		}
		catch(FileNotFoundException fnfe)
		{
			new StsMessage(currentModel.win3d, StsMessage.WARNING, "Failed to find attribute file " + attributeFilename);
			return false;
		}
	}


    /** Given a line into the screen at the mousePoint clipped to the project box,
	 *  find its XY intersection with points along the 2d line.
	 *  Compute the intersection point (XYZ) and interpolation factor along line (F)
	 *  and returned the sequence index of the nearest gather.
     *  If there is no intersection, set the interfaction value XYZF[3] to large number.
	 */

	public int getMouseSeismicLine2dIntersect(double[][] mouseLine, double[] intersectionXYZF)
	{
		double mouseLineLenSq;
		double[] lineVector;
		double[] pointVector0, pointVector1;
		double x0, y0, x1, y1;

        if(nCols < 1)
        {
            return -1;
        }
		try
		{
			mouseLineLenSq = StsMath.distanceSq(mouseLine[0], mouseLine[1], 2);
			lineVector = StsMath.subtract(mouseLine[1], mouseLine[0]);
			pointVector1 = StsMath.vector2(mouseLine[0][0], mouseLine[0][1], cdpX[0], cdpY[0]);
			y1 = computePointCoordinatesY(pointVector1, lineVector, mouseLineLenSq);
			int nPointNearest = -1;
            intersectionXYZF[3] = StsParameters.largeDouble;
            double fLineNearest = StsParameters.largeDouble;
			for (int n = 1; n < nCols; n++)
			{
				pointVector0 = pointVector1;
				pointVector1 = StsMath.vector2(mouseLine[0][0], mouseLine[0][1], cdpX[n], cdpY[n]);
				y0 = y1;
				y1 = computePointCoordinatesY(pointVector1, lineVector, mouseLineLenSq);

				if (y1 * y0 < 0.0)
				{
					x0 = computePointCoordinatesX(pointVector0, lineVector, mouseLineLenSq);
					x1 = computePointCoordinatesX(pointVector1, lineVector, mouseLineLenSq);

					double fPoint = -y0 / (y1 - y0);
					double fLine = x0 + fPoint*(x1 - x0);
                    double z = mouseLine[0][2] + fLine*(mouseLine[1][2] - mouseLine[0][2]);
                    if(fLine < fLineNearest && StsMath.betweenInclusive(z, zMin, zMax))
					{
						fLineNearest = fLine;

						if (fPoint <= 0.5)
							nPointNearest = n - 1;
						else
							nPointNearest = n;
					}
				}
			}
            if(nPointNearest == -1) return -1;
            if (!StsMath.interpolate(mouseLine[0], mouseLine[1], fLineNearest, 3, intersectionXYZF))
				return -1;
			intersectionXYZF[3] = fLineNearest;
			return nPointNearest;
		}
		catch(Exception e)
		{
			return -1;
		}
	}

    private double computePointCoordinatesX(double[] pointVector, double[] lineVector, double lineLenSq)
	{
		return StsMath.dot(lineVector, pointVector)/lineLenSq;
	}

	private double computePointCoordinatesY(double[] pointVector, double[] lineVector, double lineLenSq)
	{
		return StsMath.cross2(lineVector, pointVector)/lineLenSq;
	}

    public boolean textureChanged()
	{
	    textureChanged = true;
        return true;
    }

    public boolean dataChanged()
    {
        lineStackOk = false;
        textureChanged();
        return true;
    }

    public void geometryChanged()
    {
    }

    /** Called to actually delete the displayables on the delete list. */
	public void deleteTexturesAndDisplayLists(GL gl)
	{
		if(textureTiles == null) return;
		textureTiles.deleteTextures(gl);
		textureTiles.deleteDisplayLists(gl);
        textureTiles = null;
        textureChanged = true;
	}
/*
	public void addTextureToDeleteList()
	{
	    if (textureTiles != null)
		{
			StsTextureList.addTextureToDeleteList(this);
		}
	    textureChanged = true;
        lineStackOk = false;
    }
*/
    public void drawTextureTileSurface(StsTextureTile tile, GL gl, boolean is3d)
	{
		byte volumeZDomain = lineSet.getZDomain();
		if(volumeZDomain == StsParameters.TD_TIME)
			drawTextureTileTimeSurface(tile, gl);
		else
		{
			drawTextureTileDepthSurface(tile, gl);
		}
	}

	private void drawTextureTileTimeSurface(StsTextureTile tile, GL gl)
	{
		StsGridCrossingPoint gridCrossingPoint;

			gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
//            gl.glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
		 double rowTexCoor = 0;
		 double dRowTexCoor = tile.dRowTexCoor;
         float tileZMin, tileZMax;
         gl.glGetError();

//         if(drawBlank)         //this was incorrectly initializing the texture tiles, only first ~1 second of velocities isVisible
//         {
//            tileZMin = zMin;
//            tileZMax = zMax;
//         }
//         else
//         {
		    tileZMin = zMin + tile.colMin * zInc;
		    tileZMax = zMin + tile.colMax * zInc;
//         }
         gl.glBegin(GL.GL_QUAD_STRIP);
		 for (int row = tile.rowMin; row <= tile.rowMax; row++, rowTexCoor += dRowTexCoor)
		 {
             gl.glNormal3fv(normals[row], 0);
             float x = cdpX[row];
			 float y = cdpY[row];
			 gl.glTexCoord2d(tile.minColTexCoor, rowTexCoor);
			 gl.glVertex3f(x, y, tileZMin);
			 gl.glTexCoord2d(tile.maxColTexCoor, rowTexCoor);
			 gl.glVertex3f(x, y, tileZMax);
		 }
		 gl.glEnd();
	}

	public void drawTextureTileDepthSurface(StsTextureTile tile, GL gl)
	{
		StsGridCrossingPoint gridCrossingPoint;
		/*
		 gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
//            gl.glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
		 double rowTexCoor = 0;
//		double rowTexCoor = tile.minRowTexCoor;
		 double dRowTexCoor = tile.dRowTexCoor;
		 double dColTexCoor = tile.dColTexCoor;
		 float[] xyz = cellGridCrossingPoints[tile.rowMin].getXYZ();
		 float x1 = xyz[0];
		 float y1 = xyz[1];
		 StsSeismicVolume velocityVolume = velocityModel.getVelocityVolume();
		 int r1 = velocityVolume.getNearestBoundedRowCoor(y1);
		 int c1 = velocityVolume.getNearestBoundedColCoor(x1);
		 float[] velocityTrace1 = velocityVolume.getTraceValues(r1, c1);
		 float depthMin = velocityModel.depthDatum;

		 for (int row = tile.rowMin + 1; row <= tile.rowMax; row++, rowTexCoor += dRowTexCoor)
		 {
		  float x0 = x1;
		  float y0 = y1;
		  float[] velocityTrace0 = velocityTrace1;

		  xyz = cellGridCrossingPoints[row].getXYZ();
		  x1 = xyz[0];
		  y1 = xyz[1];
		  r1 = velocityVolume.getNearestBoundedRowCoor(y1);
		  c1 = velocityVolume.getNearestBoundedColCoor(x1);
		  velocityTrace1 = velocityVolume.getTraceValues(r1, c1);

		  gl.glBegin(GL.GL_QUAD_STRIP);

		  double colTexCoor = tile.minColTexCoor;
		  float t = tMin + tile.colMin*tInc;
		  if (lighting)
		  {
		   for (int col = tile.colMin; col <= tile.colMax; col++, t += tInc, colTexCoor += dColTexCoor)
		   {
		 float v0 = velocityTrace0[col];
		 float z0 = v0*t + depthMin;
		 gl.glTexCoord2d(colTexCoor, rowTexCoor);
		 gl.glVertex3f(x0, y0, z0);
		 float v1 = velocityTrace1[col];
		 float z1 = v1*t + depthMin;
		 gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
		 gl.glVertex3f(x1, y1, z1);
		   }
		  }
		  else
		  {
		   for (int col = tile.colMin; col <= tile.colMax; col++, t += tInc, colTexCoor += dColTexCoor)
		   {
		 float v0 = velocityTrace0[col];
		 float z0 = v0 * t + depthMin;
		 gl.glTexCoord2d(colTexCoor, rowTexCoor);
		 gl.glVertex3f(x0, y0, z0);
		 float v1 = velocityTrace1[col];
		 float z1 = v1 * t + depthMin;
		 gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
		 gl.glVertex3f(x1, y1, z1);
		   }
		  }
		 gl.glEnd();
		 }
		 */
	}

	public boolean displayTexture(StsGLPanel3d glPanel3d)
	{
		GL gl = glPanel3d.getGL();
		try
		{
			StsModel model = StsObject.getCurrentModel();

			if (lineSet.isPixelModeChanged())
			{
				// glPanel3d.deleteTexture(textureTiles, gl);
				textureChanged = true;
			}
			textureChanged = true; //forcing updates in velocity model to be drawn
			//byte projectZDomain = model.getProject().getZDomain();

			checkStackOptionChanged();
			//boolean displayVelocity = lineSet.lineSetClass.getDisplayVelocity();
			byte displayType = lineSet.lineSetClass.getDisplayType();

			checkSeismicVolumeChanged();
			if (textureTiles == null || !textureTiles.isSameSize(nCols, nSlices)) //Tom???? had to add this to avoid NPE while loading 2D
				textureTiles = StsTextureTiles.constructor(currentModel, this, nCols, nSlices, lineSet.isPixelMode);
			if(textureTiles.shaderChanged()) textureChanged = true;
/*
			useDisplayLists = model.useDisplayLists;


            if (!useDisplayLists && usingDisplayLists)
			{
				deleteDisplayLists(gl);
				usingDisplayLists = false;
			}
*/
			ByteBuffer planeData = null;
			if (textureChanged)
			{
				if (lineSet.velocityModel != null)
				{
					if (displayType == StsPreStackLineSetClass.DISPLAY_VELOCITY)
						planeData = ((StsPreStackVelocityModel2d)lineSet.velocityModel).computeVelocity2dBytes(this);
					else if (displayType == StsPreStackLineSetClass.DISPLAY_STACKED)
						planeData = getStackedByteData2d();
					else if (displayType == StsPreStackLineSetClass.DISPLAY_SEMBLANCE)
						planeData = getSemblanceByteData2d();
					else
						return false;
				}
				if (planeData != null)
				{
					drawBlank = false;
				}
				else
				{
                    drawBlank = true;
					planeData = ByteBuffer.wrap(new byte[nSlices * nCols ]); //needs to be number of cdps x number of samples
				}

//				if (textureTiles == null || !textureTiles.isSameSize(this.nCols, nCols))
//					textureTiles = StsTextureTiles.constructor(currentModel, this, nCols, nCols, lineSet.isPixelMode);
			}

            textureTiles.checkBuildDisplayLists(gl, true);
            /*
            if (textureChanged || useDisplayLists && !usingDisplayLists)
			{
				if (useDisplayLists)
					usingDisplayLists = true;
				if (textureTiles == null)
					StsException.systemError("StsSurface.displaySurfaceFill() failed. textureTiles should not be null.");
				textureTiles.constructSurface(this, gl, useDisplayLists, true);
			}
			*/

			gl.glEnable(GL.GL_TEXTURE_2D);
			gl.glEnable(GL.GL_BLEND);
			gl.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE);
			gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
			gl.glShadeModel(GL.GL_FLAT);

//            gl.glLightModeli(GL.GL_LIGHT_MODEL_COLOR_CONTROL, GL.GL_SEPARATE_SPECULAR_COLOR);

			if (drawBlank)
			{
				if (textureTiles.shader != StsJOGLShader.NONE)
				{
					StsJOGLShader.loadEnableARBColormap(gl, StsColor.GREY.getRGBAArray(), textureTiles.shader);
				}
			}
			else
			{
//                gl.glDisable(GL.GL_LIGHTING);
				if (!lineSet.setGLColorList(gl, true, displayType, textureTiles.shader))
					return false;
			}

			textureTiles.displayTiles(this, gl, lineSet.isPixelMode, planeData, nullByte);
			textureChanged = false;
			return true;
		}
		catch (Exception e)
		{
			StsException.outputWarningException(this, "displayTexture(GLPanel3d)", "", e);
			return false;
		}
		finally
		{
			gl.glDisable(GL.GL_TEXTURE_2D);
			gl.glDisable(GL.GL_BLEND);
//			gl.glEnable(GL.GL_LIGHTING);
			if (textureTiles.shader != StsJOGLShader.NONE)
				StsJOGLShader.disableARBShader(gl);
		}
	}

    /*
    private boolean checkDisplayVelocityToggleChanged()
    {
        boolean currentDisplayVelocity = lineSet.lineSetClass.displayVelocity;
        if (currentDisplayVelocity == displayVelocity) return false;
		this.displayVelocity = currentDisplayVelocity;
        textureChanged = true;
        return true;
    }
*/
    private boolean checkStackOptionChanged()
    {
		byte currentStackOption = lineSet.lineSetClass.getStackOption();
		if (currentStackOption == stackOption) return false;
		stackOption = currentStackOption;
        textureChanged = true;
        return true;
	}

    private boolean checkSeismicVolumeChanged()
	{
		StsPreStackLineSet currentLineSet = lineSet.lineSetClass.getCurrentProjectLineSet();
		if(lineSet == currentLineSet)return false;
        textureChanged = true;
        if(currentLineSet == null)return true;
        lineSet = currentLineSet;
        return true;
	}

    public ByteBuffer getStackedByteData2d()
	{
        StsPreStackVelocityModel2d velocityModel = (StsPreStackVelocityModel2d)lineSet.velocityModel;
        if(velocityModel == null)return null;
		if(!velocityModel.hasProfiles(this))return null;
		switch(stackOption)
		{
			case StsPreStackLineSet.STACK_NONE:
				return null;
			case StsPreStackLineSet.STACK_NEIGHBORS:
				ByteBuffer lineBuffer = getStackedLine2d(true);
				System.out.println("Neighbor 2d not up ");
				//return getNeighborBuffer(dirNo, nPlane, lineBuffer);
			case StsPreStackLineSet.STACK_LINES:
				return getStackedLine2d(false);
			default:
				return null;
		}
	}

	private ByteBuffer getStackedLine2d(boolean neighborsOnly)
	{
		MappedByteBuffer byteBlockBuffer = null;
		ByteBuffer lineBuffer = null; // buffer contains full line of traces
		int row = -1, col = -1;
        StsPreStackVelocityModel2d velocityModel = (StsPreStackVelocityModel2d)lineSet.velocityModel;
//		if(!velocityInterpolationComplete())return null;
        System.out.println("Getting stacked line "+ getName());

		try
		{
			if(debugTimer)
			{
				StsSeismicTimer.clear();
				StsSeismicTimer.overallTimer.start();
				StsSeismicTimer.getOutputBlockTimer.start();
			}

 //           FloatBuffer floatBlockBuffer = lineSet.fileMapRowFloatBlocks.getByteBufferPlane(lineIndex, FileChannel.MapMode.READ_WRITE).asFloatBuffer();
			if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
            lineBuffer = ((StsPreStackLineSet2d)lineSet).fileMapBlocks2d.getByteBufferPlane(lineIndex, FileChannel.MapMode.READ_WRITE);
			if(lineStackOk)return lineBuffer;

			if(debug)System.out.println("  ByteBuffer capacity: " + lineBuffer.capacity());
			byte[] data = new byte[nSlices];
			StsSuperGather superGather = new StsSuperGather(currentModel, lineSet); // tempory superGather used for stacking operation

			StsPreStackLineSet.checkTransparentTrace(nSlices);

			float maxZToStack = lineSet.lineSetClass.getMaxZToStack();
			int nMaxDisplaySlices = Math.min(nSlices, getNearestBoundedSliceCoor(maxZToStack)) + 1;
			if(nMaxDisplaySlices < nSlices)
				System.arraycopy(StsPreStackLineSet.byteTransparentTrace, nMaxDisplaySlices, data, nMaxDisplaySlices, nSlices - nMaxDisplaySlices);

			int position = 0;
			StsStatusArea statusArea = currentModel.win3d.statusArea;
            int nUpdatedCols = 0;

            statusArea.setStatus("2d line stacking...");
            statusArea.addProgress();

            int nCols = this.nCols;
            statusArea.setMaximum(nCols);
            row = lineIndex;
            for(col = 0; col < nCols; col++)
            {
                if(velocityModel.interpolation.isNeighbor[row][col] || !neighborsOnly)
                {
                    if(!superGather.initializeSuperGather(row, col))
                    {
                        lineBuffer.put(StsPreStackLineSet.byteTransparentTrace);
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
                        continue;
                    }
                    // already filtered: shouldn't do it again. TJL 7/15/06
                    float[] gatherData = superGather.centerGather.computeStackedTrace(nMaxDisplaySlices);
                    if(gatherData == null)
                    {
//							//if(debug)
                        System.out.println("NULL superGather! lineBuffer.put transparent 2d trace for row " + row + " col " + col);
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.start();
                        lineBuffer.put(StsPreStackLineSet.byteTransparentTrace);
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
                    }
                    else
                    {
                        StsSeismicFilter.normalizeAmplitude(gatherData, nMaxDisplaySlices);   
                        nUpdatedCols++;

                        for(int s = 0; s < nMaxDisplaySlices; s++)
                            data[s] = (byte)(127 + 127 * gatherData[s]);
//								if(debug) System.out.println("lineBuffer.put data trace for row " + row + " col " + col);
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.start();
                        lineBuffer.put(data);
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
                    }
                }
                else // PLANE_STATUS_DIRTY && gridChanged && !isNeighbor && neighborsOnly
                {
                    if(debugTimer) StsSeismicTimer.getOutputBlockTimer.start();
                    lineBuffer.put(StsPreStackLineSet.byteTransparentTrace, 0, nSlices);
                    if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
                    //					System.arraycopy(transparentTrace, 0, data, n, nCroppedSlices);
                    //					n += nCroppedSlices;
                }
                statusArea.setProgress(col);
            }
            if(debug)System.out.println("Updated " + nUpdatedCols + " out of " + nCols + " cols for row " + row);
            statusArea.removeProgress();
            statusArea.clearStatus();
			lineStackOk = true;

			//        progressBarDialog.hide();
			//		System.out.println("Finished stacking line....");
//			blockBuffer.force();
			lineBuffer.rewind();
			if(debugTimer)
			{
				StsSeismicTimer.printTimers(" Line index " + getName());
				StsSeismicTimer.overallTimer.stopPrint("Total time: ");
			}
			return lineBuffer;
		}
		catch(Exception e)
		{
			e.printStackTrace();
			if(lineBuffer != null)
				StsException.systemError("StsPreStackSeismicVolume.getStackedByteData() failed. \n" +
										 "byteBuffer capacity: " + lineBuffer.capacity() + " position: " + lineBuffer.position());
			else
				StsException.outputException("StsPreStackSeismicVolume.getStackedByteData() failed: null lineBuffer.", e, StsException.WARNING);
			return null;
		}
	}
    public ByteBuffer getSemblanceByteData2d()
	{
        StsPreStackVelocityModel2d velocityModel = (StsPreStackVelocityModel2d)lineSet.velocityModel;
        if(velocityModel == null)return null;
		if(!velocityModel.hasProfiles(this))return null;
		switch(stackOption)
		{
			case StsPreStackLineSet.STACK_NONE:
				return null;
			case StsPreStackLineSet.STACK_NEIGHBORS:
				ByteBuffer lineBuffer = getStackedLine2d(true);
				System.out.println("Neighbor 2d not up ");
				//return getNeighborBuffer(dirNo, nPlane, lineBuffer);
			case StsPreStackLineSet.STACK_LINES:
				return getSemblanceLine2d(false);
			default:
				return null;
		}
	}

	private ByteBuffer getSemblanceLine2d(boolean neighborsOnly)
	{
		MappedByteBuffer byteBlockBuffer = null;
		ByteBuffer lineBuffer = null; // buffer contains full line of traces
		int row = -1, col = -1;
        StsPreStackVelocityModel2d velocityModel = (StsPreStackVelocityModel2d)lineSet.velocityModel;
//		if(!velocityInterpolationComplete())return null;
        System.out.println("Getting stacked line "+ getName());

		try
		{
			if(debugTimer)
			{
				StsSeismicTimer.clear();
				StsSeismicTimer.overallTimer.start();
				StsSeismicTimer.getOutputBlockTimer.start();
			}

 //           FloatBuffer floatBlockBuffer = lineSet.fileMapRowFloatBlocks.getByteBufferPlane(lineIndex, FileChannel.MapMode.READ_WRITE).asFloatBuffer();
			if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
            lineBuffer = ((StsPreStackLineSet2d)lineSet).fileMapBlocks2d.getByteBufferPlane(lineIndex, FileChannel.MapMode.READ_WRITE);
			if(lineStackOk)return lineBuffer;

			if(debug)System.out.println("  ByteBuffer capacity: " + lineBuffer.capacity());
			byte[] data = new byte[nSlices];
			StsSuperGather superGather = new StsSuperGather(currentModel, lineSet); // tempory superGather used for stacking operation

			StsPreStackLineSet.checkTransparentTrace(nSlices);

			float maxZToStack = lineSet.lineSetClass.getMaxZToStack();
			int nMaxDisplaySlices = Math.min(nSlices, getNearestBoundedSliceCoor(maxZToStack)) + 1;
			if(nMaxDisplaySlices < nSlices)
				System.arraycopy(StsPreStackLineSet.byteTransparentTrace, nMaxDisplaySlices, data, nMaxDisplaySlices, nSlices - nMaxDisplaySlices);

			int position = 0;
			StsStatusArea statusArea = currentModel.win3d.statusArea;
            int nUpdatedCols = 0;

            statusArea.setStatus("2d line stacking...");
            statusArea.addProgress();

            int nCols = this.nCols;
            statusArea.setMaximum(nCols);
            row = lineIndex;
            for(col = 0; col < nCols; col++)
            {
                if(velocityModel.interpolation.isNeighbor[row][col] || !neighborsOnly)
                {
                    if(!superGather.initializeSuperGather(row, col))
                    {
                        lineBuffer.put(StsPreStackLineSet.byteTransparentTrace);
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
                        continue;
                    }
                    float[] gatherData = superGather.computeSemblanceTrace(nMaxDisplaySlices);
                    if(gatherData == null)
                    {
//							//if(debug)
                        System.out.println("NULL superGather! lineBuffer.put transparent 2d trace for row " + row + " col " + col);
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.start();
                        lineBuffer.put(StsPreStackLineSet.byteTransparentTrace);
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
                    }
                    else
                    {
                        nUpdatedCols++;

                        for(int s = 0; s < nMaxDisplaySlices; s++)
                        {
                            int iValue = (int)Math.round(254 * StsMath.minMax(gatherData[s], 0.0, 1.0));
                            data[s] = StsMath.unsignedIntToUnsignedByte(iValue);
                        }
//								if(debug) System.out.println("lineBuffer.put data trace for row " + row + " col " + col);
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.start();
                        lineBuffer.put(data);
                        if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
                    }
                }
                else // PLANE_STATUS_DIRTY && gridChanged && !isNeighbor && neighborsOnly
                {
                    if(debugTimer) StsSeismicTimer.getOutputBlockTimer.start();
                    lineBuffer.put(StsPreStackLineSet.byteTransparentTrace, 0, nSlices);
                    if(debugTimer) StsSeismicTimer.getOutputBlockTimer.stopAccumulateIncrementCount();
                    //					System.arraycopy(transparentTrace, 0, data, n, nCroppedSlices);
                    //					n += nCroppedSlices;
                }
                statusArea.setProgress(col);
            }
            if(debug)System.out.println("Updated " + nUpdatedCols + " out of " + nCols + " cols for row " + row);
            statusArea.removeProgress();
            statusArea.clearStatus();
			lineStackOk = true;

			//        progressBarDialog.hide();
			//		System.out.println("Finished stacking line....");
//			blockBuffer.force();
			lineBuffer.rewind();
			if(debugTimer)
			{
				StsSeismicTimer.printTimers(" Line index " + getName());
				StsSeismicTimer.overallTimer.stopPrint("Total time: ");
			}
			return lineBuffer;
		}
		catch(Exception e)
		{
			e.printStackTrace();
			if(lineBuffer != null)
				StsException.systemError("StsPreStackSeismicVolume.getSemblanceByteData() failed. \n" +
										 "byteBuffer capacity: " + lineBuffer.capacity() + " position: " + lineBuffer.position());
			else
				StsException.outputException("StsPreStackSeismicVolume.getSemblanceByteData() failed: null lineBuffer.", e, StsException.WARNING);
			return null;
		}
	}

    public void display2dLine(StsGLPanel3d glPanel3d, GL gl)
	{
		if(cdpX == null || cdpY == null) return;

        gl.glDisable(GL.GL_LIGHTING);

        gl.glColor3f(0.8f, 0.f, 0.f);
        gl.glLineWidth(2.f);
        gl.glBegin(GL.GL_LINE_STRIP);
        for(int n = 0; n < nCols; n++)
        {
            gl.glVertex3f(cdpX[n], cdpY[n], zMin);
        }
        for(int n = 0; n < nCols; n++)
        {
            gl.glVertex3f(cdpX[nCols - n - 1], cdpY[nCols - n - 1], zMax);
        }
        gl.glVertex3f(cdpX[0], cdpY[0], zMin);
        gl.glEnd();

        gl.glColor3f(1.0f, 0.0f, 0.0f);
        gl.glPointSize(1.f);
        gl.glBegin(GL.GL_POINTS);
        for(int n = 0; n < nCols; n++)
        {
            gl.glVertex3f(cdpX[n], cdpY[n], zMin);
        }
        for(int n = 0; n < nCols; n++)
        {
            gl.glVertex3f(cdpX[n], cdpY[n], zMax);
        }
        gl.glEnd();
        gl.glEnable(GL.GL_LIGHTING);

        if(displayTexture(glPanel3d))
            return;
        else
        {
            gl.glDisable(GL.GL_LIGHTING);

            gl.glColor3f(0.8f, 0.8f, 0.8f);
            gl.glBegin(GL.GL_QUAD_STRIP);

            for(int n = 0; n < nCols; n++)
            {
                gl.glNormal3fv(normals[n], 0);
                gl.glVertex3f(cdpX[n], cdpY[n], zMin);
                gl.glVertex3f(cdpX[n], cdpY[n], zMax);
            }
            gl.glEnd();

            gl.glEnable(GL.GL_LIGHTING);
        }
	}

	public void drawVerticalLineAtCDP(GL gl, int nCDP)
	{
        // hack until we get bogus value fixed
        nCDP = StsMath.minMax(nCDP, 0, nCols-1);
        StsPoint[] points = new StsPoint[2];
	    points[0] = new StsPoint(cdpX[nCDP], cdpY[nCDP], zMin);
		points[1] = new StsPoint(cdpX[nCDP], cdpY[nCDP], zMax);
		StsGLDraw.drawDottedLine(gl, StsColor.WHITE, 2, points, StsColor.BLACK);
	}

    public String getGroupname()
    {
        return StsSeismicBoundingBox.group2dPrestack;
    }

    public Class getDisplayableClass() { return StsPreStackLineSet2d.class; }
    
    public float getAnalysisColStart() {
        if (analysisColStart < colNumMin) setAnalysisColStart(colNumMin);
        return analysisColStart;
    }

    public void setAnalysisColStart(float analysisColStart) 
    {
    	this.analysisColStart = analysisColStart;
    	if (analysisColStart < colNumMin) analysisColStart = colNumMin;
    	if (analysisColStart > colNumMax) analysisColStart = colNumMax;
    	dbFieldChanged("analysisColStart", analysisColStart);
    }

    public int getAnalysisColInc() {return analysisColInc;}

    public void setAnalysisColInc(int analysisColInc) 
    {
    	this.analysisColInc = analysisColInc;
    	dbFieldChanged("analysisColInc", analysisColInc);
    }

}
