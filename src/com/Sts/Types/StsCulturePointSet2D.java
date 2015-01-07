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
public class StsCulturePointSet2D extends StsCultureObject2D
{
    double[][] points = null;

    public StsCulturePointSet2D()
    {
        setStsColor(StsColor.RED);
        symbolType = PT2D;
    }

    public StsCulturePointSet2D(double[][] points)
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
        xyz = currentModel.getProject().getRotatedRelativeXYFromUnrotatedAbsoluteXY(points[index][0], points[index][1]);
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
        draw(cultureDisplayable, glPanel3d, mapToSurface, zDomainOriginal, true);
    }
    public void draw2d(StsCultureDisplayable cultureDisplayable, boolean mapToSurface, byte zDomainOriginal, StsGLPanel glPanel)
    {
        draw(cultureDisplayable, glPanel, mapToSurface, zDomainOriginal, false);
    }

    public void draw(StsCultureDisplayable cultureDisplayable, StsGLPanel glPanel, boolean mapToSurface, byte zDomainOriginal, boolean is3d)
    {
    	try
    	{
    		GL gl = glPanel.getGL();
    		int nPoints = points.length;
    		boolean isPlanar = cultureDisplayable.isPlanar();
    		glPanel.setViewShift(gl, StsGraphicParameters.gridShift);
            // Want to plot as flat plane if not original and no velocity model.
    		if(zDomainOriginal != currentModel.getProject().getZDomain() && currentModel.getProject().velocityModel == null)
                mapToSurface = true;    		
    		float[] fPt = new float[3];
    		int zIndex = 3;
    		if(currentModel.getProject().getZDomain() == currentModel.getProject().TD_DEPTH)
    			zIndex = 2;    		
    		for(int n = 0; n < nPoints; n++)
    		{
    			float[] xy = currentModel.getProject().getRotatedRelativeXYFromUnrotatedAbsoluteXY(points[n][0], points[n][1]);
    			fPt[0] = xy[0];
    			fPt[1] = xy[1];
    			if((mapToSurface) || (!mapToSurface && isPlanar()))
    			{
    				if(isPlanar)
    					fPt[2] = currentModel.getCursor3d().getCurrentDirCoordinate(StsCursor3d.ZDIR);
    				else
    					fPt[2] = cultureDisplayable.getCultureZ(fPt[0],fPt[1]);
    			}
    			else
    			{
        			if(zDomainOriginal != currentModel.getProject().getZDomain())
        			{
        				if(zDomainOriginal == StsProject.TD_DEPTH)
        					fPt[2] = (float) currentModel.getProject().velocityModel.getT(fPt[0], fPt[1], (float)points[n][2], 0.0f);
        				else
        					fPt[2] = (float) currentModel.getProject().velocityModel.getZ(fPt[0], fPt[1], (float)points[n][3]);
        			}
        			else
        				fPt[2] = (float)points[n][zIndex];
    			}

                if(!is3d)
                    StsGLDraw.drawPoint(gl, fPt, stsColor, size);
                else
                {
                    switch(symbolType)
                    {
                        case PT2D:
                            StsGLDraw.drawPoint(gl, fPt, stsColor, size);
                            break;
                        case SPHERE:
                            StsGLDraw.drawSphere((StsGLPanel3d)glPanel, fPt, stsColor, size);
                            break;
                        case CUBE:
                            StsGLDraw.drawCube((StsGLPanel3d)glPanel, fPt, stsColor, size);
                            break;
                        case CYLINDER:
                            StsGLDraw.drawCylinder((StsGLPanel3d)glPanel, fPt, stsColor, size, size);
                            break;
                    }
                }
            }
    		glPanel.resetViewShift(gl);
    	}
    	catch(Exception ex)
    	{
    		StsException.outputException("Error drawing cultural point set", ex, StsException.WARNING);
    	}   	
    }

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
