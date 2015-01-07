package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */
public class StsCultureLine2D extends StsCultureObject2D
{

    public static final byte SOLID = 0;
    public static final byte DASHED = 1;

    double[][] points = null;
    int stroke = SOLID;

    public StsCultureLine2D()
    {
        setStsColor(StsColor.RED);
    }

    public StsCultureLine2D(double[][] points)
    {
        this.points = points;
        setStsColor(StsColor.RED);
    }

    public StsPoint getPointAt(int index)
    {
    	if(points == null)
    		return null;
    	if(points.length-1 < index)
    		return null;
    	
    	float[] xyz = new float[3];
        xyz = currentModel.getProject().getRelativeXY(points[index][0], points[index][1]);
        return new StsPoint(xyz);
    }
    public int getNumPoints()
    {
    	if(points == null)
    		return 0;
    	
    	return points.length;
    }
    public void addPoint(double x, double y, double z, double t)
    {
        if(points == null)
        {
            points = new double[][] { {x, y, z, t} };
        }
        else
        {
            int nPoints = points.length;
            double[][] newPoints = new double[nPoints+1][];
            System.arraycopy(points, 0, newPoints, 0, nPoints);
            newPoints[nPoints] = new double[] { x, y, z, t };
            points = newPoints;
        }
    }
    /**
     * draw
     *
     * @param cultureDisplayable StsCultureDisplayable
     * @param glPanel3d
     * @todo Implement this com.Sts.Types.StsCultureObject2D method
     */
    public void draw(StsCultureDisplayable cultureDisplayable, boolean mapToSurface, byte zDomainOriginal, StsGLPanel3d glPanel3d)
    {
        draw(cultureDisplayable, glPanel3d, mapToSurface, zDomainOriginal, D3, true);
    }
    public void draw2d(StsCultureDisplayable cultureDisplayable, boolean mapToSurface, byte zDomainOriginal, StsGLPanel glPanel)
    {
        draw(cultureDisplayable, glPanel, mapToSurface, zDomainOriginal, D2, false);
    }

    public void draw(StsCultureDisplayable cultureDisplayable, StsGLPanel glPanel, boolean mapToSurface, byte zDomainOriginal, byte dimensions, boolean is3d)
    {
//        gl.glLineWidth(StsGraphicParameters.edgeLineWidth);
    	try
    	{
    		int nPoints = points.length;
            GL gl = glPanel.getGL();
    		boolean isPlanar = cultureDisplayable.isPlanar();
    		glPanel.setViewShift(gl, StsGraphicParameters.gridShift);
            // Want to plot as flat plane if not original and no velocity model.
    		if(zDomainOriginal != currentModel.getProject().getZDomain() && currentModel.getProject().velocityModel == null)
                mapToSurface = true;       
    		int zIndex = 3;
    		if(currentModel.getProject().getZDomain() == currentModel.getProject().TD_DEPTH)
    			zIndex = 2;
    		if(dimensions == D3)
    		{
    			if((mapToSurface) || (!mapToSurface && isPlanar()))
    			{
    				if(isPlanar)
    				{
    					StsPoint[] pts = new StsPoint[nPoints];
    					float z = currentModel.getCursor3d().getCurrentDirCoordinate(StsCursor3d.ZDIR);
    					for (int n = 0; n < nPoints; n++)
    					{
    						pts[n] = new StsPoint(points[n]);
    						float[] xy = currentModel.getProject().getRotatedRelativeXYFromUnrotatedAbsoluteXY(points[n][0], points[n][1]);
    						pts[n].v[0] = xy[0];
    						pts[n].v[1] = xy[1];
    						pts[n].v[2] = z;                        
    					}

    					if(stroke == SOLID)
    						StsGLDraw.drawLine(gl, stsColor, false, pts, 2);
    					else
    						StsGLDraw.drawDottedLine(gl, stsColor, false, pts, 2);
    				}
    				else
    				{
    					StsPoint[] pts = new StsPoint[nPoints];
    					for (int n = 0; n < nPoints; n++)
    					{
    						pts[n] = new StsPoint(points[n]);
    						float[] xy = currentModel.getProject().getRotatedRelativeXYFromUnrotatedAbsoluteXY(points[n][0], points[n][1]);
    						pts[n].v[0] = xy[0];
    						pts[n].v[1] = xy[1];
    						pts[n].v[2] = cultureDisplayable.getCultureZ(xy[0], xy[1]);
    					}
    					StsGLDraw.drawLine(gl, stsColor, false, pts, 2);
    				}
    			}
    			else
    			{
    				StsPoint[] pts = new StsPoint[nPoints];
    				float[] xyz = new float[4];                
    				for (int n = 0; n < nPoints; n++)
    				{
    					float[] xy = currentModel.getProject().getRotatedRelativeXYFromUnrotatedAbsoluteXY(points[n][0], points[n][1]);
    					xyz[0] = xy[0];
    					xyz[1] = xy[1];
    					if(zDomainOriginal != currentModel.getProject().getZDomain())
    					{
    						if(zDomainOriginal == StsProject.TD_DEPTH)
    							xyz[2] = (float) currentModel.getProject().velocityModel.getT(xyz[0], xyz[1], (float)points[n][2], 0.0f);
    						else
    							xyz[2] = (float) currentModel.getProject().velocityModel.getZ(xyz[0], xyz[1], (float)points[n][3]);
    					}
    					else
    						xyz[2] = (float)points[n][zIndex];
    					
    					pts[n] = new StsPoint(xyz);				
    				}
    				if(stroke == SOLID)
    					StsGLDraw.drawLine(gl, stsColor, false, pts, 2);
    				else
    					StsGLDraw.drawDottedLine(gl, stsColor, false, pts, 2);
    			}
    		}
    		else
    		{
				StsPoint[] pts = new StsPoint[nPoints];
				float z = currentModel.getCursor3d().getCurrentDirCoordinate(StsCursor3d.ZDIR);
				for (int n = 0; n < nPoints; n++)
				{
					pts[n] = new StsPoint(points[n]);
					float[] xy = currentModel.getProject().getRotatedRelativeXYFromUnrotatedAbsoluteXY(points[n][0], points[n][1]);
					pts[n].v[0] = xy[0];
					pts[n].v[1] = xy[1];
					pts[n].v[2] = z;
                    gl.glLineWidth(StsGraphicParameters.well3dLineWidth);
					if(n > 0)
					{
						if(stroke == SOLID)
						    StsGLDraw.drawLine2d(gl, stsColor, pts[n-1].getPointXYZ(), pts[n].getXYZ());
						else
							StsGLDraw.drawDottedLine2d(gl, stsColor, StsColor.BLACK, width, pts[n-1], pts[n]);
					}
				}
    		}
    		glPanel.resetViewShift(gl);
    	}
    	catch(Exception ex)
    	{
    		StsException.outputException("Error drawing cultural lines", ex, StsException.WARNING);
    	}
    }
    public int getStroke() { return stroke; }
    public void setStroke(int stroke) { this.stroke = stroke; }

    public float getMaxDepth(byte zDomainOriginal)
    {
        float maxDepth = currentModel.getProject().getDepthMax();
        float depth;
        for(int i=0; i<points.length; i++)
        {
			if(zDomainOriginal != currentModel.getProject().getZDomain())
			{
				if(currentModel.getProject().velocityModel != null)
					depth = (float) currentModel.getProject().velocityModel.getZ((float)points[i][0], (float)points[i][1], (float)points[i][3]);
				else
					depth = currentModel.getProject().getDepthMax();
			}
			else
				depth = (float)points[i][2];
			
            if(depth > maxDepth)
                maxDepth = depth;
        }
        return maxDepth;
    }
    
    public float getMaxTime(byte zDomainOriginal)
    {
    	try
    	{
    		float maxTime = currentModel.getProject().getTimeMax();
    		float time;
    		for(int i=0; i<points.length; i++)
    		{
    			if(zDomainOriginal != currentModel.getProject().getZDomain())
    			{
    				if(currentModel.getProject().velocityModel != null)
    					time = (float) currentModel.getProject().velocityModel.getT((float)points[i][0], (float)points[i][1], (float)points[i][2], 0.0f);
    				else
    					time = currentModel.getProject().getTimeMax();
    			}
    			else
    				time = (float)points[i][3];
			
    			if(time > maxTime)
    				maxTime = time;
    		}
    		return maxTime;
    	}
    	catch(Exception ex)
    	{
    		StsException.outputException("Error calculating maxTime", ex, StsException.WARNING);
    		return currentModel.getProject().getTimeMin();
    	}        
    }
    
     public float getMinDepth(byte zDomainOriginal)
    {
         float minDepth = currentModel.getProject().getDepthMin();
         float depth;
         for(int i=0; i<points.length; i++)
         {
 			if(zDomainOriginal != currentModel.getProject().getZDomain())
 			{
 				if(currentModel.getProject().velocityModel != null)
 					depth = (float) currentModel.getProject().velocityModel.getZ((float)points[i][0], (float)points[i][1], (float)points[i][3]);
 				else
 					depth = currentModel.getProject().getDepthMin();
 			}
 			else
 				depth = (float)points[i][2];
 			
             if(depth < minDepth)
                 minDepth = depth;
         }
         return minDepth;    	 
    }
     
    public float getMinTime(byte zDomainOriginal)
    {
    	try
    	{
    		float minTime = currentModel.getProject().getTimeMin();
    		float time;
    		for(int i=0; i<points.length; i++)
    		{
    			if(zDomainOriginal != currentModel.getProject().getZDomain())
    			{
    				if(currentModel.getProject().velocityModel != null)
    					time = (float) currentModel.getProject().velocityModel.getT((float)points[i][0], (float)points[i][1], (float)points[i][2], 0.0f);
    				else
    					time = currentModel.getProject().getTimeMin();
    			}
    			else
    				time = (float)points[i][3];
			
    			if(time < minTime)
    				minTime = time;
        	}
        	return minTime;
    	}
    	catch(Exception ex)
    	{
    		StsException.outputException("Error calculating minTime", ex, StsException.WARNING);
    		return currentModel.getProject().getTimeMin();
    	}
    }

    public double getXMax()
    {
        double xMax = -StsParameters.largeFloat;
        for(int i=0; i<points.length; i++)
        {
            if(points[i][0] > xMax)
                xMax = points[i][0];
        }
        return xMax;
    }

    public double getYMax()
    {
        double yMax = -StsParameters.largeFloat;
        for(int i=0; i<points.length; i++)
        {
            if(points[i][1] > yMax)
                yMax = points[i][1];
        }
        return yMax;
    }

    public double getXMin()
    {
        double xMin = -StsParameters.largeFloat;
        for(int i=0; i<points.length; i++)
        {
            if(points[i][0] < xMin)
                xMin = points[i][0];
        }
        return xMin;
    }

    public double getYMin()
    {
        double yMin = -StsParameters.largeFloat;
        for(int i=0; i<points.length; i++)
        {
            if(points[i][1] < yMin)
                yMin = points[i][1];
        }
        return yMin;
    }
}
