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
public class StsCultureText extends StsCultureObject2D
{
    double[] point = new double[4];
    String text = null;

    public StsCultureText()
    {
        setStsColor(StsColor.WHITE);
    }

    public StsCultureText(double[] point, String text)
    {
        this.point = point;
        this.text = text;
        setStsColor(StsColor.WHITE);
    }

    public void setPoint(double x, double y, double z, double t)
    {
        point[0] = x;
        point[1] = y;
        point[2] = z;
        point[3] = t;
    }
    
    public StsPoint getPointAt(int index)
    {
    	float[] xyz = new float[3];
        xyz = currentModel.getProject().getRotatedRelativeXYFromUnrotatedAbsoluteXY(point[0], point[1]);
        return new StsPoint(xyz);    	
    }
    
    public int getNumPoints()
    {
    	return 1;
    }
    
    public void setText(String text)
    {
        this.text = text;
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
    		boolean isPlanar = cultureDisplayable.isPlanar();
    		glPanel.setViewShift(gl, StsGraphicParameters.gridShift);
            // Want to plot as flat plane if not original and no velocity model.
    		if(zDomainOriginal != currentModel.getProject().getZDomain() && currentModel.getProject().velocityModel == null)
                mapToSurface = true;    		
    		float[] fPt = currentModel.getProject().getRotatedRelativeXYFromUnrotatedAbsoluteXY(point[0], point[1]);
    		float[] xyz = new float[3];
    		xyz[0] = fPt[0];
    		xyz[1] = fPt[1];
    		int zIndex = 3;
    		if(currentModel.getProject().getZDomain() == currentModel.getProject().TD_DEPTH)
    			zIndex = 2;      	
    		if((mapToSurface) || (!mapToSurface && isPlanar()))
    		{
    			if(isPlanar)
    				xyz[2] = currentModel.getCursor3d().getCurrentDirCoordinate(StsCursor3d.ZDIR);
    			else
    				xyz[2] = cultureDisplayable.getCultureZ(xyz[0], xyz[1]);
    		}
    		else
    		{
    			if(zDomainOriginal != currentModel.getProject().getZDomain())
    			{
    				if(zDomainOriginal == StsProject.TD_DEPTH)
    					xyz[2] = (float) currentModel.getProject().velocityModel.getT(xyz[0], xyz[1], (float)point[2], 0.0f);
    				else
    					xyz[2] = (float) currentModel.getProject().velocityModel.getZ(xyz[0], xyz[1], (float)point[3]);
    			} 
    			else
    				xyz[2] = (float)point[zIndex];
    		}
    		StsGLDraw.drawPoint(gl, xyz, getStsColor(), 2);
    		if(size < 16)
    			StsGLDraw.fontHelvetica12(gl, xyz, text);
    		else
    			StsGLDraw.fontHelvetica18(gl, xyz, text);
    		glPanel.resetViewShift(gl);
    	}
    	catch(Exception ex)
    	{
    		StsException.outputException("Error drawing cultural text", ex, StsException.WARNING);
    	}		
    }
    
    public float getMaxDepth(byte zDomainOriginal)
    {
    	float depth = currentModel.getProject().getDepthMax();
		if(zDomainOriginal != currentModel.getProject().getZDomain())
		{
			if(currentModel.getProject().velocityModel != null)
				depth = (float) currentModel.getProject().velocityModel.getZ((float)point[0], (float)point[1], (float)point[3]);
			else
				depth = currentModel.getProject().getDepthMax();
		}
		else
			depth = (float)point[2];

        return depth;    	
    }
    
    public float getMaxTime(byte zDomainOriginal)
    {
    	float time = currentModel.getProject().getTimeMax();
    	try
    	{
    		if(zDomainOriginal != currentModel.getProject().getZDomain())
    		{
    			if(currentModel.getProject().velocityModel != null)
    				time = (float) currentModel.getProject().velocityModel.getT((float)point[0], (float)point[1], (float)point[2], 0.0f);
    			else
    				time = currentModel.getProject().getTimeMax();
    		}
    		else
    			time = (float)point[3];

        return time;
    	}
    	catch(Exception ex)
    	{
    		StsException.outputException("Error calculating maxTime", ex, StsException.WARNING);
    		return currentModel.getProject().getTimeMax();
    	}
    }
    public float getMinDepth(byte zDomainOriginal)
    {
    	float depth = currentModel.getProject().getDepthMin();
		if(zDomainOriginal != currentModel.getProject().getZDomain())
		{
			if(currentModel.getProject().velocityModel != null)
				depth = (float) currentModel.getProject().velocityModel.getZ((float)point[0], (float)point[1], (float)point[3]);
			else
				depth = currentModel.getProject().getDepthMax();
		}
		else
			depth = (float)point[2];

        return depth;
    }
    
    public float getMinTime(byte zDomainOriginal)
    {
    	float time = currentModel.getProject().getTimeMax();
    	try
    	{
    		if(zDomainOriginal != currentModel.getProject().getZDomain())
    		{
    			if(currentModel.getProject().velocityModel != null)
    				time = (float) currentModel.getProject().velocityModel.getT((float)point[0], (float)point[1], (float)point[2], 0.0f);
    			else
    				time = currentModel.getProject().getTimeMin();
    		}
    		else
    			time = (float)point[3];

        return time;
    	}
    	catch(Exception ex)
    	{
    		StsException.outputException("Error calculating maxTime", ex, StsException.WARNING);
    		return currentModel.getProject().getTimeMin();
    	}    	
    }

    public double getXMax()
    {
        return point[0];
    }
    public double getYMax()
    {
        return point[1];
    }
     public double getXMin()
    {
        return point[0];
    }
    public double getYMin()
    {
        return point[1];
    }
}
