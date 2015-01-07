package com.Sts.DBTypes;

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;

import javax.media.opengl.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2003</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

public class StsHorpickClass extends StsClass implements StsSerializable, StsClassDisplayable, StsClassCursorDisplayable
{
    public StsHorpickClass()
    {
        userName = "Autotracked Surfaces";
    }

    public void initializeDisplayFields()
    {
//        initColors(StsHorpick.displayFields);

        displayFields = new StsFieldBean[]
        {
            new StsBooleanFieldBean(this, "isVisible", "Enable")
        };
        super.setIsVisible(false);
    }

    public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate,
                               boolean axesFlipped, boolean xAxisReversed, boolean axisReversed)
    {
        if(currentObject == null) return;
        if(!isVisible) return;
		GL gl = glPanel3d.getGL();
//		glPanel3d.setViewShift(gl, StsGraphicParameters.vertexOnEdgeShift);
        ((StsHorpick)currentObject).display2d(gl, dirNo, dirCoordinate, axesFlipped);
//		glPanel3d.resetViewShift(gl);
        /*
        int nElements = this.getSize();
        GL gl = glPanel3d.getGL();
        for(int n = 0; n < nElements; n++)
        {
            StsHorpick horpick = (StsHorpick)getElement(n);
            horpick.display2d(gl, dirNo, dirCoordinate, axesFlipped);
        }
        */
    }

    public void drawOnCursor3d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, StsPoint[] planePoints, boolean isDragging)
    {
    }

    public void drawOn3dCurtain(StsGLPanel3d glPanel3d, StsSeismicCurtain seismicCurtain)
    {
    }

    public void displayClass(StsGLPanel3d glPanel3d)
    {
        if(currentObject == null) return;
        if(!isVisible) return;
        ((StsHorpick)currentObject).display(glPanel3d);
    }

    public StsHorpick getHorpickWithSurface(StsSurface surface)
    {
        for(int n = 0; n < getSize(); n++)
        {
            StsHorpick horpick = (StsHorpick) getElement(n);
            if (horpick.getSurface() == surface)return horpick;
        }
        return null;
    }



    public void setIsVisible(boolean isVisible)
    {
        super.setIsVisible(isVisible);
        currentModel.win3dDisplay();
    }
}
