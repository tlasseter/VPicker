package com.Sts.Types;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Utilities.*;
import com.magician.fonts.*;

import javax.media.opengl.*;
import java.io.*;
import java.text.*;

/** The displayed bounding box for the project.  Constructed by taking the boundingBox
 *  around the project and expanding it to "rounded" values of xMin, xMax, yMin, and yMax
 */

public class StsDisplayBoundingBox extends StsBoundingBox implements Cloneable, Serializable
{
	/** displayed grid spacing in x */
    public float gridDX;
	/** displayed grid spacing in y */
    public float gridDY;
	/** z value at which grid is displayed */
	public float gridZ = StsParameters.nullValue;
	/** Approx number of XY Grid increments */
	static public int approxNumXYGridIncrements = 20;

//    public StsVerticalFont verticalFont = new StsVerticalFont(new GLHelvetica12BitmapFont());
//    public DecimalFormat labelFormat = new DecimalFormat("#,##0.0#");

    public StsDisplayBoundingBox()
    {
    }

	public StsDisplayBoundingBox(boolean persistent)
	{
		super(persistent);
	}

	public StsDisplayBoundingBox(boolean persistent, String name)
	{
		super(persistent, name);
	}

    public boolean initialize(StsModel model)
    {
        return true;
    }

	public void setZRange(float zMin, float zMax)
	{
		this.zMin = zMin;
		this.zMax = zMax;
		gridZ = zMax;
	}

	/* Starting with current boundingBox, adjust x, y, and z ranges so they are "nice".
	 * @parameters model the current StsModel object
	 * @see setValuesFromBoundingBox
	 * @see StsMath.niceScale2
	 *  */
	public void adjustToNiceSize()
	{
		float dSize = 0.0f;
		// if bounding box is tall and skinny, extend size of box laterally
		float xSize = xMax - xMin;
		float ySize = yMax - yMin;
		float zSize = zMax - zMin;
		if (zSize > 2.0f * xSize || zSize > 2.0f * ySize) dSize = 0.5f*zSize;
		xMin -= dSize;
		xMax += dSize;
		yMin -= dSize;
		yMax += dSize;

		// find nice integral range for x and y
		double[] xScale = StsMath.niceScale(xMin + xOrigin, xMax + xOrigin, approxNumXYGridIncrements, true);
		xMin = (float) (xScale[0] - xOrigin);
		xMax = (float) (xScale[1] - xOrigin);
		gridDX = (float) xScale[2];
		double[] yScale = StsMath.niceScale(yMin + yOrigin, yMax + yOrigin, approxNumXYGridIncrements, true);
		yMin = (float) (yScale[0] - yOrigin);
		yMax = (float) (yScale[1] - yOrigin);
		gridDY = (float) yScale[2];
	}

	/* Starting with current boundingBox, adjust x, y, and z ranges so they are "nice".
	 * @parameters model the current StsModel object
	 * @see setValuesFromBoundingBox
	 * @see StsMath.niceScale2
	 *  */
	public void adjustBoundingBoxXYRange(StsBoundingBox boundingBox)
	{
		xOrigin = boundingBox.xOrigin;
		yOrigin = boundingBox.yOrigin;
		xMin = boundingBox.xMin;
		xMax = boundingBox.xMax;
		yMin = boundingBox.yMin;
		yMax = boundingBox.yMax;
		zMin = boundingBox.zMin;
		zMax = boundingBox.zMax;
		float dSize = 0.0f;

		// if display bounding box is tall and skinny, extend size of box laterally
		float xSize = xMax - xMin;
		float ySize = yMax - yMin;
		float zSize = zMax - zMin;
		if (zSize > 2.0f * xSize || zSize > 2.0f * ySize) dSize = 0.5f*zSize;
		xMin -= dSize;
		xMax += dSize;
		yMin -= dSize;
		yMax += dSize;

		// find nice integral range for x and y
		double[] xScale = StsMath.niceScale(xMin + xOrigin, xMax + xOrigin, approxNumXYGridIncrements, true);
		xMin = (float) (xScale[0] - xOrigin);
		xMax = (float) (xScale[1] - xOrigin);
		gridDX = (float) xScale[2];
		double[] yScale = StsMath.niceScale(yMin + yOrigin, yMax + yOrigin, approxNumXYGridIncrements, true);
		yMin = (float) (yScale[0] - yOrigin);
		yMax = (float) (yScale[1] - yOrigin);
		gridDY = (float) yScale[2];
	}

	/* Starting with current boundingBox, adjust x, y, and z ranges so they are "nice".
	 * @parameters model the current StsModel object
	 * @see setValuesFromBoundingBox
	 * @see StsMath.niceScale2
	 *  */
/*
	public void adjustBoundingBoxToNiceSize(StsRotatedGridBoundingBox boundingBox)
	{
		xOrigin = boundingBox.xOrigin;
		yOrigin = boundingBox.yOrigin;
		computeUnrotatedBoundingBox(boundingBox);
		zMin = boundingBox.zMin;
		zMax = boundingBox.zMax;

		float dSize = 0.0f;

		// if bounding box is tall and skinny, extend size of box laterally
		float xSize = xMax - xMin;
		float ySize = yMax - yMin;
		float zSize = zMax - zMin;
		if (zSize > 2.0f * xSize || zSize > 2.0f * ySize) dSize = 0.5f*zSize;
		xMin -= dSize;
		xMax += dSize;
		yMin -= dSize;
		yMax += dSize;

		// find nice integral range for x and y
		double[] xScale = StsMath.niceScale(xMin + xOrigin, xMax + xOrigin, approxNumXYGridIncrements, true);
		xMin = (float) (xScale[0] - xOrigin);
		xMax = (float) (xScale[1] - xOrigin);
		gridDX = (float) xScale[2];
		double[] yScale = StsMath.niceScale(yMin + yOrigin, yMax + yOrigin, approxNumXYGridIncrements, true);
		yMin = (float) (yScale[0] - yOrigin);
		yMax = (float) (yScale[1] - yOrigin);
		gridDY = (float) yScale[2];

		double[] zScale = StsMath.niceScale(zMin, zMax, approxNumZGridIncrements, true);
		zMin = (float)(zScale[0]);
		zMax = (float)(zScale[1]);
		gridZ = zMax;
		boundingBox.zInc = (float)zScale[2];
	}
*/
	public void addUnrotatedBoundingBox(StsBoundingBox boundingBox)
	{
		super.addUnrotatedBoundingBox(boundingBox);
		gridZ = zMax;
	}
	/** Starting with current boundingBox, adjust x, y, and z grid lines so they are "nice".
	 * @parameters model the current StsModel object
	 * @see setValuesFromBoundingBox
	 * @see StsMath.niceScale2
	 */
/*
	public void adjustBoundingBoxGridLines(StsRotatedGridBoundingBox boundingBox)
	{
		// find nice integral range for x and y
		double xOrigin = boundingBox.xOrigin;
		double[] xScale = StsMath.niceScale2(xMin+xOrigin, xMax+xOrigin, approxNumXYGridIncrements, true);
		gridDX = (float)xScale[2];
		double yOrigin = boundingBox.yOrigin;
		double[] yScale = StsMath.niceScale2(yMin+yOrigin, yMax+yOrigin, approxNumXYGridIncrements, true);
		gridDY = (float)yScale[2];
		double[] zScale = StsMath.niceScale2(zMin, zMax, approxNumZGridIncrements, true);
		zMin = (float)(zScale[0]);
		zMax = (float)(zScale[1]);
		gridZ = zMax;
		boundingBox.zInc = (float)zScale[2];
	}
*/
	/** Display method for Project */
	public void display(StsGLPanel3d glPanel3d, StsModel model, float angle)
	{
		GL gl = glPanel3d.getGL();
		if(gl == null)return;

		try
		{
			if(model == null)return;
			StsProject project = model.getProject();

			if(angle != 0.0f)
			{
                gl.glPushMatrix();
				gl.glRotatef(angle, 0.0f, 0.0f, -1.0f);
			}

			float x, y;
			float vec[] = new float[3];

			/** IMPORTANT: turn LIGHTING off and then back on for lines */
			gl.glDisable(GL.GL_LIGHTING);

			gl.glLineWidth(StsGraphicParameters.gridLineWidth);

			gl.glLineStipple(1, StsGraphicParameters.dottedLine);
			gl.glEnable(GL.GL_LINE_STIPPLE);

			/* All lines are at grid_z */
			if(project.getShowGrid())
			{
				/** draw lines with gray color */
				project.getGridColor().setGLColor(gl);

				vec[2] = gridZ;

				gl.glBegin(GL.GL_LINES);
				{

					/* Draw the grid lines in the x-direction */
					for(x = xMin + gridDX; x < xMax; x += gridDX)
					{
						vec[0] = x;
						vec[1] = yMin;
						gl.glVertex3fv(vec, 0);
						vec[0] = x;
						vec[1] = yMax;
						gl.glVertex3fv(vec, 0);
					}

					/* Draw the grid lines in the y-direction */
					for(y = yMin + gridDY; y < yMax; y += gridDY)
					{
						vec[0] = xMin;
						vec[1] = y;
						gl.glVertex3fv(vec, 0);
						vec[0] = xMax;
						vec[1] = y;
						gl.glVertex3fv(vec, 0);
					}
				}
				gl.glEnd();
			} // showGrid
			gl.glDisable(GL.GL_LINE_STIPPLE);
			if(!project.getIsVisible())return;

            // if(StsGLPanel.debugProjectionMatrix) glPanel3d.debugPrintProjectionMatrix("StsDisplayBoundingBox.display(). proj matrix before displayBoundingBox call: ");
            displayBoundingBox(gl, project.getStsGridColor(), StsGraphicParameters.gridLineWidth);
            // if(StsGLPanel.debugProjectionMatrix) glPanel3d.debugPrintProjectionMatrix("StsDisplayBoundingBox.display(). proj matrix after displayBoundingBox call: ");
			/* Draw
			  //	set_line_smoothing(current_window, FALSE);

			   /* draw dotted line from grid to lookat_point 	*/
			 /* and a square at lookat_point				 	*/
			 /*
			   // draw line from grid to centerViewPoint with cyan
			   StsColor.setGLColor(gl, StsColor.CYAN);

			   gl.glLineStipple(1, DOTTED_LINE);
			   gl.glEnable(GL.GL_LINE_STIPPLE);

			   StsPoint centerViewPoint = glPanel3d.getCenterViewPoint();

			   gl.glLineWidth(4.0f);
			   gl.glBegin(GL.GL_LINE_STRIP);
			   {
				vec[0] = centerViewPoint.v[0]; vec[1] = centerViewPoint.v[1]; vec[2] = gridZ;
				gl.glVertex3fv(vec);
				vec[0] = centerViewPoint.v[0]; vec[1] = centerViewPoint.v[1]; vec[2] = centerViewPoint.v[2];
				gl.glVertex3fv(vec);
			   }
			   gl.glEnd();

			   gl.glDisable(GL.GL_LINE_STIPPLE);

			   StsGLDraw.drawPoint(vec, StsColor.BLACK, glPanel3d, 6, -1.0);
			   StsGLDraw.drawPoint(vec, StsColor.CYAN, glPanel3d, 4, 0.0);
			  */
			 if(project.getShowLabels())
			 {
				 String label1, label2, label3, label4;
				 GLBitmapFont horizontalFont = GLHelvetica12BitmapFont.getInstance(gl);
				 DecimalFormat labelFormat = new DecimalFormat("###0");

				 labelFormat = model.getProject().getLabelFormat();

				 // Draw Edge Labels
		//			if(!glPanel3d.getCursor3d().isGridCoordinates())  // Actual Coordinates
				 {
					 label1 = new String("Xmin=" + labelFormat.format(xMin + xOrigin) + " Ymin=" + labelFormat.format(yMin + yOrigin));
					 label2 = new String("Xmax=" + labelFormat.format(xMax + xOrigin));
					 label3 = new String("Ymax=" + labelFormat.format(yMax + yOrigin));
		//				label4 = new String("Z=" + labelFormat.format(zMax));
				 }
				 /*
					else  // Relative Coordinates
					{
					 label1 = new String("XL=" + labelFormat.format(glPanel3d.cursor3d.currentSeismicVolume.getLabelFromCoor(StsCursor3d.XDIR,xMin)) + " IL=" + labelFormat.format(yMin));
					 label2 = new String("XL=" + labelFormat.format(glPanel3d.cursor3d.currentSeismicVolume.getLabelFromCoor(StsCursor3d.XDIR,xMax)));
					 label3 = new String("YL=" + labelFormat.format(glPanel3d.cursor3d.currentSeismicVolume.getLabelFromCoor(StsCursor3d.YDIR,yMax)));
		//				label4 = new String("Z=" + labelFormat.format(glPanel3d.cursor3d.currentSeismicVolume.getLabelFromCoor(StsCursor3d.ZDIR,zMax)));
					}
				  */
				 vec[0] = xMin;
				 vec[1] = yMin;
				 vec[2] = zMin;
				 StsGLDraw.fontOutput(gl, vec, label1, horizontalFont);
				 vec[0] = xMax;
				 StsGLDraw.fontOutput(gl, vec, label2, horizontalFont);
				 vec[0] = xMin;
				 vec[1] = yMax;
				 StsGLDraw.fontOutput(gl, vec, label3, horizontalFont);
			 }
			 else
			 {
				 GLBitmapFont horizontalFont = GLHelvetica18BitmapFont.getInstance(gl);
				 vec[2] = zMin;

				 vec[0] = xMax;
				 vec[1] = yMin;
				 StsGLDraw.fontOutput(gl, vec, "E", horizontalFont);
				 vec[0] = xMin;
				 vec[1] = yMax;
				 StsGLDraw.fontOutput(gl, vec, "N", horizontalFont);
			 }
		}
		catch(Exception e)
		{
			StsException.outputException("StsDisplayBoundingBox.display() failed.", e, StsException.WARNING);
		}
		finally
		{
            if(angle != 0.0f) gl.glPopMatrix();
            gl.glEnable(GL.GL_LIGHTING);
		}
	}
/*
	private void drawCursorCubeEdges(StsGLPanel3d glPanel3d, StsModel model, boolean isPicking)
	{
		GL gl = glPanel3d.getGL();
		if(gl == null) return;

		StsProject project = model.getProject();

		float vec[] = new float[3];

		String label1, label2, label3, label4;
		GLBitmapFont horizontalFont = new GLHelvetica12BitmapFont();
		DecimalFormat labelFormat = new DecimalFormat("###0");

		gl.glDisable(GL.GL_LIGHTING);

		gl.glLineWidth(1.0f);

		// draw lines with gray color
		StsColor.setGLColor(gl, project.getStsGridColor());

		gl.glDisable(GL.GL_LINE_STIPPLE);
	//	set_line_smoothing(current_window, TRUE);

		// Draw the bottom outer boundaries

		vec[2] = zMin;

		float xExt = (xMax - xMin) * 0.05f;
		float yExt = (yMax - yMin) * 0.05f;
		float zExt = (zMax - zMin) * 0.05f;

		gl.glBegin(GL.GL_LINE_LOOP);
		{
			vec[0] = xMin; vec[1] = yMin;
			gl.glVertex3fv(vec);
			vec[0] = xMin; vec[1] = yMax;
			gl.glVertex3fv(vec);
			vec[0] = xMax; vec[1] = yMax;
			gl.glVertex3fv(vec);
			vec[0] = xMax; vec[1] = yMin;
			gl.glVertex3fv(vec);
		}
		gl.glEnd();

		// Draw the top outer boundaries

		vec[2] = zMax;

		gl.glBegin(GL.GL_LINE_LOOP);
		{
			vec[0] = xMin; vec[1] = yMin;
			gl.glVertex3fv(vec);
			vec[0] = xMin; vec[1] = yMax;
			gl.glVertex3fv(vec);
			vec[0] = xMax; vec[1] = yMax;
			gl.glVertex3fv(vec);
			vec[0] = xMax; vec[1] = yMin;
			gl.glVertex3fv(vec);
		}
		gl.glEnd();

		if(project.getShowLabels())
		{
			StsSeismicVolumeClass lineSetClass = (StsSeismicVolumeClass)model.getCreateStsClass(StsSeismicVolume.class);
			if(lineSetClass != null)
				labelFormat = new DecimalFormat(lineSetClass.getLabelFormat());
			else
				labelFormat = new DecimalFormat("###0");

			// Draw Edge Labels
			if(!glPanel3d.getCursor3d().isGridCoordinates())  // Actual Coordinates
			{
				label1 = new String("X=" + labelFormat.format(xMin) + " Y=" +
								   labelFormat.format(yMin) + " Z=" +
								   labelFormat.format(zMin));
				label2 = new String("X=" + labelFormat.format(xMax));
				label3 = new String("Y=" + labelFormat.format(yMax));
				label4 = new String("Z=" + labelFormat.format(zMax));
			}
			else  // Relative Coordinates
			{
				label1 = new String("XL=" + labelFormat.format(glPanel3d.cursor3d.currentSeismicVolume.getLabelFromCoor(StsCursor3d.XDIR,xMin)) + " IL=" + labelFormat.format(yMin) +
									" Z=" + labelFormat.format(glPanel3d.cursor3d.currentSeismicVolume.getLabelFromCoor(StsCursor3d.ZDIR,zMin)));
				label2 = new String("XL=" + labelFormat.format(glPanel3d.cursor3d.currentSeismicVolume.getLabelFromCoor(StsCursor3d.XDIR,xMax)));
				label3 = new String("YL=" + labelFormat.format(glPanel3d.cursor3d.currentSeismicVolume.getLabelFromCoor(StsCursor3d.YDIR,yMax)));
				label4 = new String("Z=" + labelFormat.format(glPanel3d.cursor3d.currentSeismicVolume.getLabelFromCoor(StsCursor3d.ZDIR,zMax)));
			}
			vec[0] = xMin;
			vec[1] = yMin;
			vec[2] = zMin - zExt;
			StsGLDraw.fontOutput(gl, vec, label1, horizontalFont);
			vec[0] = xMax + xExt;
			vec[2] = zMin;
			StsGLDraw.fontOutput(gl, vec, label2, horizontalFont);
			vec[0] = xMin;
			vec[1] = yMax + yExt;
			StsGLDraw.fontOutput(gl, vec, label3, horizontalFont);
			vec[1] = yMin;
			vec[2] = zMax;
			StsGLDraw.fontOutput(gl, vec, label4, horizontalFont);
		}

		// Draw verticals

		gl.glBegin(GL.GL_LINES);
		{
			vec[0] = xMin; vec[1] = yMin;
			vec[2] = zMin - zExt;
			gl.glVertex3fv(vec);
			vec[2] = zMax;
			gl.glVertex3fv(vec);

			vec[0] = xMin; vec[1] = yMax;
			vec[2] = zMin;
			gl.glVertex3fv(vec);
			vec[2] = zMax;
			gl.glVertex3fv(vec);

			vec[0] = xMax; vec[1] = yMax;
			vec[2] = zMin;
			gl.glVertex3fv(vec);
			vec[2] = zMax;
			gl.glVertex3fv(vec);

			vec[0] = xMax; vec[1] = yMin;
			vec[2] = zMin;
			gl.glVertex3fv(vec);
			vec[2] = zMax;
			gl.glVertex3fv(vec);

			// Draw Extensions at X and Y Origin
			vec[2] = zMin;
			vec[0] = xMin - xExt; vec[1] = yMin;
			gl.glVertex3fv(vec);
			vec[0] = xMin;
			gl.glVertex3fv(vec);

			vec[0] = xMin; vec[1] = yMin - yExt;
			gl.glVertex3fv(vec);
			vec[1] = yMin;
			gl.glVertex3fv(vec);
		}
		gl.glEnd();

		gl.glEnable(GL.GL_LIGHTING);

	//	set_line_smoothing(current_window, FALSE);
	}
*/
}
