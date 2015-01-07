package com.Sts.DBTypes;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.PreStack.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.io.*;

public class StsPreStackLine3d extends StsPreStackLine implements Comparable
{
    /** number label for this line; used by 3d only */
    public float lineNum = -StsParameters.largeFloat;

    public String getGroupname()
    {
        return StsSeismicBoundingBox.group3dPrestack;
    }

    public StsPreStackLine3d()
	{
	}

	public StsPreStackLine3d(boolean persistent)
	{
		super(persistent);
	}

	/*
	 public StsPreStackLine3d(StsFile file, StsModel model) throws FileNotFoundException
	 {
	  this(file.getDirectory(), file.createFilename(), model);
	 }
	 */
	private StsPreStackLine3d(StsModel model, StsFile file, StsPreStackLineSet3d lineSet3d) throws FileNotFoundException, StsException
	{
		super(false);
		if (!file.exists())
		{
			throw new FileNotFoundException();
		}
		if(!StsParameterFile.initialReadObjectFields(file.getPathname(), this, StsPreStackLine.class, StsBoundingBox.class))
        {
//            throw new StsException("Failed reading file: " + file.getFilename());
        }

        if(lineSet3d.stsDirectory == null) lineSet3d.stsDirectory = stsDirectory;
        this.lineSet = lineSet3d;
		initializeScale();
		addToModel();
		if (!initialize(model))
		{
			throw new StsException("Failed initializing: " + getName());
		}
        stsDirectory = lineSet3d.stsDirectory;
        isVisible = true;
	}

	static public StsPreStackLine3d constructor(StsModel model, StsFile file, StsPreStackLineSet3d lineSet3d)
	{
		try
		{
			return new StsPreStackLine3d(model, file, lineSet3d);
		}
		catch (FileNotFoundException e)
		{
			return null;
		}
		catch (Exception e)
		{
			StsException.outputException("StsPreStackLine3d.constructor() failed.", e, StsException.WARNING);
			return null;
		}
	}
/*
	public boolean classInitialize(StsModel model)
	{
		try
		{
//			traceOffsets = getAttributeArray(stsDirectory, getAttributeFilename(), "OFFSET");
			StsBlocksMemoryManager memoryManager = currentModel.project.getBlocksMemoryManager();
			fileBlocks = new StsPreStackLineFileBlocks(nGathers, nOffsetsMax, nLastTraceInGathers, nCroppedSlices, stsDirectory, gatherFilename, 4, memoryManager);
			return true;
		}
		catch(Exception e)
		{
			StsException.systemError("StsPreStackLine3d.classInitialize(model) failed.");
			return false;
		}
	}
*/
    public double[] getAttributeArray(String attributeName, int row, int col) throws FileNotFoundException
    {
        int start = 0;
        int idx = 0, num = 0;

        double[] attArray = getAttributeArray(attributeName);
        if(attArray == null)
            return null;

        if(isInline)
            idx = col;
        else
            idx = row;

        start = getNFirstGatherTrace(idx);
        num = getNTracesInGather(idx);
        if(num == 0)
            return null;

        double[] newArray = new double[num];
        for(int i=0; i<num; i++)
            newArray[i] = attArray[start+i];

        return newArray;
	}

	public long getGatherFileOffset(int row, int col)
	{
		if(isInline)
			return getGatherFileOffset(col);
		else
			return getGatherFileOffset(row);
	}

	public int getGatherFileLength(int row, int col)
	{
		if(isInline)
			return getGatherFileLength(col);
		else
			return getGatherFileLength(row);
	}
/*
	public int getPlaneValue(float[] xyz)
	{
		return -1;
	}
*/
	public int compareTo(Object other)
	{
		StsPreStackLine3d otherLine3d = (StsPreStackLine3d) other;

		if(lineIndex < otherLine3d.lineIndex) return -1;
		if(lineIndex > otherLine3d.lineIndex) return 1;
		if(maxGatherIndex < otherLine3d.minGatherIndex) return -1;
		if(minGatherIndex > otherLine3d.maxGatherIndex) return 1;
		StsException.systemError("StsPreStackLine3d.compareTo() failed.  Lines overlap for lineIndex " + lineIndex + " firstLineGatherIndices " + minGatherIndex + "-" + maxGatherIndex +
			 " otherLineIndex " + otherLine3d.lineIndex + " secondLineGatherIndices " + otherLine3d.minGatherIndex + "-" + otherLine3d.maxGatherIndex);
	    return 0;
	}

	public int compareTo(int otherLineIndex, int otherGatherIndex)
	{
		if(lineIndex < otherLineIndex) return -1;
		if(lineIndex > otherLineIndex) return 1;
		if(maxGatherIndex < otherGatherIndex) return -1;
		if(minGatherIndex > otherGatherIndex) return 1;
		return 0;
	}

	public void initializeLineIndex(StsPreStackLineSet3d seismicVolume)
	{
		if (isInline)
		{
			lineNum = rowNumMin;
			lineIndex = seismicVolume.getNearestRowCoor(yMin);
 //           System.out.println("yMin=" + yMin + " Linenum=" + lineNum + " LineIndex=" + lineIndex);
		}
		else
		{
			lineNum = colNumMin;
			lineIndex = seismicVolume.getNearestColCoor(xMin);
		}
	}

	public void initializeGatherIndexes()
	{
		if (isInline)
		{
			minGatherIndex = lineSet.getColFromColNum(colNumMin);
			maxGatherIndex = lineSet.getColFromColNum(colNumMax);
		}
		else
		{
			minGatherIndex = lineSet.getRowFromRowNum(rowNumMin);
			maxGatherIndex = lineSet.getRowFromRowNum(rowNumMax);
		}
	}

    public String toString()
	{
        if (this.isInline)
			return new String(getName() + " inline " + lineNum + " crosslines " + colNumMin + " to " + colNumMax);
		else
			return new String(getName() + " crossline " + lineNum + " inlines " + rowNumMin + " to " + rowNumMax);
	}

	/*
	  public void displayWiggleTraces(GL gl, float[][] axisRanges, byte[] planeData)
	  {
	 float horizMin = axisRanges[0][0];
	 float height = axisRanges[1][0] - axisRanges[1][1];
	 float verticalMin = axisRanges[1][1];
	 int nTraces = nLastGatherTrace - nFirstGatherTrace + 1;
	 float horizInc = (axisRanges[0][1] - axisRanges[0][0])/nTraces;

	 // If density is less than 1:4 traces to pixels, display wiggles
//		if(((axisRanges[0][1] - axisRanges[0][0])/horizInc) > glPanel3d.winRectGL.width/getWiggleToPixelRatio())
//			return;

	 gl.pushMatrix();
	 gl.glLoadIdentity();
	 float horizScale = 3*horizInc / 254;
	 float verticalScale = height/(nCroppedSlices - 1);
	 horizInc /= horizScale;
	 gl.glScalef(horizScale, verticalScale, 1.0f);
	 gl.glTranslatef(horizMin/horizScale, verticalMin/verticalScale, 0.0f);
	 StsColor.setGLColor(gl, StsColor.BLACK);
	 gl.glDisable(GL.GL_LIGHTING);
	 gl.glLineWidth(0.5f);
	 boolean isDrawing = false;
	 int n = 0;
	 for (int t = 0; t < nTraces; t++)
	 {
	  for (int s = 0; s < nCroppedSlices; s++, n++)
	  {
	   if (planeData[n] == -1)
	   {
	 if (isDrawing)
	 {
	  gl.glEnd();
	  isDrawing = false;
	 }
	 continue;
	   }
	   else if (!isDrawing)
	   {
	 isDrawing = true;
	 gl.glBegin(GL.GL_LINE_STRIP);
	   }
	   gl.glVertex2i(StsMath.unsignedByteToSignedInt(planeData[n]), s);
	  }
	  if (isDrawing) gl.glEnd();
	  isDrawing = false;
	  gl.glTranslatef(horizInc, 0.0f, 0.0f);
	 }
	 gl.popMatrix();
	 gl.glEnable(GL.GL_LIGHTING);
	  }
	 */

	public void drawTextureTileSurface(StsTextureTile tile, GL gl, boolean is3d)
	{
		byte projectZDomain = StsObject.getCurrentModel().getProject().getZDomain();
//		byte volumeZDomain = getZDomain();
		byte volumeZDomain = StsParameters.TD_TIME;

		if (projectZDomain == StsParameters.TD_TIME)
		{
			if (volumeZDomain == StsParameters.TD_TIME)
				drawTextureTileTimeSurface(tile, gl, is3d);
		}
		else if (projectZDomain == StsParameters.TD_DEPTH)
		{
			// volumeZDomain is TD_TIME
			StsModel model = StsObject.getCurrentModel();
			if (volumeZDomain == StsParameters.TD_DEPTH) // seismic already in depth, don't need to convert so draw as if in time
				drawTextureTileTimeSurface(tile, gl, is3d);
			else
			{
				StsPreStackVelocityModel velocityVolume = lineSet.velocityModel;
				if (velocityVolume == null)return;
				drawTextureTileDepthSurface(velocityVolume, tile, gl);
			}
		}
	}

	public void drawTextureTileTimeSurface(StsTextureTile tile, GL gl, boolean is3d)
	{
		if (is3d)
			tile.drawQuadSurface3d(gl, StsProject.TD_TIME);
		else
			tile.drawQuadSurface2d(gl);
	}

	public void drawTextureTileDepthSurface(StsPreStackVelocityModel velocityModel, StsTextureTile tile, GL gl)
	{
	/*
		float cursorXInc, cursorYInc;
		if (dirNo == StsCursor3d.ZDIR)
		{
			return;
		}

		gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
		double rowTexCoor = tile.minRowTexCoor;
		double dRowTexCoor = tile.dRowTexCoor;
		double dColTexCoor = tile.dColTexCoor;
		double[] xyz = tile.xyzPlane[0];
		double x1 = xyz[0];
		double y1 = xyz[1];
		double t1 = xyz[2];
		int volumeRow = this.getNearestBoundedRowCoor( (float)y1);
		int volumeCol = getNearestBoundedColCoor( (float)x1);

//		StsSeismicVolume velocityVolume = velocityModel.getVelocityVolume();
//		if (velocityVolume == null)return;
//		float depthMin = velocityModel.depthDatum;
		float depthMin = 0.0f;
		int volumeRowInc = 0;
		int volumeColInc = 0;
		if (dirNo == XDIR)
		{
			cursorXInc = 0;
			cursorYInc = getYInc();
			volumeRowInc = 1;
		}
		else // dirNo == StsCursor3d.YDIR
		{
			cursorXInc = getXInc();
			cursorYInc = 0;
			volumeColInc = 1;
		}
		double tInc = getZInc();
		for (int row = tile.rowMin + 1; row <= tile.rowMax; row++, rowTexCoor += dRowTexCoor)
		{
			double x0 = x1;
			double y0 = y1;
			x1 += cursorXInc;
			y1 += cursorYInc;

			gl.glBegin(GL.GL_QUAD_STRIP);

			double colTexCoor = tile.minColTexCoor;
			double t = t1 + tile.colMin * tInc;

			for (int col = tile.colMin; col <= tile.colMax; col++, t += tInc, colTexCoor += dColTexCoor)
			{
				float v0 = velocityVolume.getValue(volumeRow, volumeCol, col);
				float z0 = (float) (v0 * t + depthMin);
				gl.glTexCoord2d(colTexCoor, rowTexCoor);
				gl.glVertex3d(x0, y0, z0);
				float v1 = velocityVolume.getValue(volumeRow + volumeRowInc, volumeCol + volumeColInc, col);
				float z1 = (float) (v1 * t + depthMin);
				gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
				gl.glVertex3d(x1, y1, z1);
			}
			gl.glEnd();
			volumeRow += volumeRowInc;
			volumeCol += volumeColInc;
		}
	*/
	}
}
