package com.Sts.UI;

import com.Sts.MVC.View3d.*;

import javax.swing.event.*;
import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jul 3, 2008
 * Time: 11:52:20 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsNewCursor3dPanel extends StsCursor3dPanel
{
    StsCursor3dTest cursor3dTest;
    public StsNewCursor3dPanel()
    {

    }

    public StsNewCursor3dPanel(Dimension size, StsCursor3dTest cursor3d)
    {
        try
        {
            this.cursor3dTest = cursor3d;
            constructPanel(size);
            setSliderValues();

            for (int n = 0; n < 3; n++)
                sliders[n].addToGroup(sliderGroup);
        }
        catch (Exception ex)
        {
            ex.printStackTrace();
        }
    }

	public void setSliderValues()
    {
        doSetValues();
//        if(rotatedBoundingBox.zMin == StsParameters.largeFloat || rotatedBoundingBox.zMax == -StsParameters.largeFloat) return;
//        StsToolkit.runLaterOnEventThread(new Runnable() { public void run() { doSetValues(); }});
    }

    public void doSetValues()
    {
        float cursorX = cursor3dTest.getCurrentDirCoordinate(StsCursor3d.XDIR);
        float cursorY = cursor3dTest.getCurrentDirCoordinate(StsCursor3d.YDIR);
        float cursorZ = cursor3dTest.getCurrentDirCoordinate(StsCursor3d.ZDIR);
        setSliderValues(cursorX, cursorY, cursorZ);
    }

    public void setSliderValues(float cursorX, float cursorY, float cursorZ)
    {
        sliderX.initSliderValues(cursor3dTest.xMin, cursor3dTest.xMax, cursor3dTest.xInc, cursorX);
        sliderY.initSliderValues(cursor3dTest.yMin, cursor3dTest.yMax, cursor3dTest.yInc, cursorY);
        sliderZ.initSliderValues(cursor3dTest.zMin, cursor3dTest.zMax, cursor3dTest.zInc, cursorZ);
        sliderX.setValueLabel("X");
        sliderY.setValueLabel("Y");
    }

    public void stateChanged(ChangeEvent e)
    {
    }
    public void setSliderSelected(StsSliderBean slider, boolean enable)
    {
    }

    public void buttonSelected(StsSliderBean slider, boolean enable)
    {
        slider.setButtonSelected(enable);
        if(!enable) return;
        setSliderSelected(slider, enable);
        if(!slider.isSelected()) slider.setCheckBoxModelSelected(true);
    }
}
