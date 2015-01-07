package com.Sts.DBTypes;

import com.Sts.Interfaces.StsClassTimeDisplayable;
import com.Sts.Interfaces.StsClassTimeSeriesDisplayable;
import com.Sts.Interfaces.StsClassViewSelectable;
import com.Sts.Interfaces.StsDialogFace;
import com.Sts.MVC.StsModel;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.StsOkApplyCancelDialog;
import com.Sts.UI.StsSensorProperties;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;

import java.util.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsDynamicSensorClass extends StsSensorClass implements StsClassTimeDisplayable, StsClassTimeSeriesDisplayable, StsClassViewSelectable
{
   boolean displayGather = false;
   int zLimit = 5;  // Default is 5% of Z range
   float accentScale = 1.1f;
   boolean enableZLimit = false;
   boolean enable = true;

   public StsDynamicSensorClass()
   {
       userName = "Variable Location Sensors";
   }
   public void initializeFields()
   {
       displayFields = new StsFieldBean[]
       {
           new StsBooleanFieldBean(this, "enable", "Enable"),
           new StsBooleanFieldBean(this, "showDeleted", "Show Deleted"),
           new StsBooleanFieldBean(this, "displayGather", "Display Gather"),
           new StsBooleanFieldBean(this, "displayNames", "Names"),
           new StsBooleanFieldBean(this, "enableTime", "Enable Time"),
           new StsBooleanFieldBean(this, "enableSound", "Sound Accent"),
           new StsBooleanFieldBean(this, "enableAccent", "Graphic Accent"),
           new StsFloatFieldBean(this, "accentScale", 0.05f, 10.0f, "Accent Scale:", true),
           new StsFloatFieldBean(this, "scale2D", 0.1f, 1.0f, "2D Scale Factor:",true),
           new StsBooleanFieldBean(this, "enableGoTo", "Go To New Data"),
           new StsIntFieldBean(this, "goToOffset", 1, 100, "Go To Offset (*scale):", true),
           new StsBooleanFieldBean(this,"enableZLimit", "Enable T/Z Limit"),
           new StsIntFieldBean(this, "zLimit", 1, 100, "T/Z Limit(%):", true),
       };

       defaultFields = new StsFieldBean[]
       {
           new StsColorComboBoxFieldBean(this, "accentColor", "Accent Color:", StsColor.colors32),
           new StsComboBoxFieldBean(this, "defaultSpectrumName", "Spectrum:", StsSpectrumClass.cannedSpectrums)
       };
   }

    public StsObject[] getSensors()
    {

        return getObjectList();
    }

    public void setDisplayGather(boolean val)
    {
        displayGather = val;
        setDisplayField("displayGather", displayGather);
    }

    public boolean getDisplayGather()
    {
        return displayGather;
    }
    public void setZLimit(int z)
    {
        zLimit = z;
        currentModel.win3dDisplayAll();
    }

    public float computeZTCriteria()
    {
        Iterator iter = getVisibleObjectIterator();
        float zMin = StsParameters.largeFloat;
        float zMax = -StsParameters.largeFloat;
        while (iter.hasNext())
        {
            StsSensor sensor = (StsSensor) iter.next();
            float[] zRange = sensor.getZTRange();
            if (zRange == null)
                continue;
            if (zRange[0] < zMin)
                zMin = zRange[0];
            if (zRange[1] > zMax)
                zMax = zRange[1];
        }
        if (zLimit > 0.0)
            return (zMax - zMin) * ((zLimit / 100.0f) / 2.0f);
        else
            return (zMax - zMin) / 2.0f;
    }

    public float getZCriteria() { return zCriteria; }

    public int getZLimit() { return zLimit; }

    public float getAccentScale() { return accentScale; }
    public void setAccentScale(float scale) { accentScale = scale; }

    public void setEnableZLimit(boolean val)
    {
        enableZLimit = val;
        currentModel.win3dDisplayAll();
    }

    public boolean getEnableZLimit() { return enableZLimit; }

    /** Display seismic data related to the selected sensor event. */
    // TODO this needs to be reworked for multipanel window since we don't know if/where we have a gatherView
    public void displayGather(StsDynamicSensor sensor, int eventIdx, StsWin3dBase win3d)
    {
        /*
        if (displayGather)
        {
            StsViewGatherMicroseismic view = null;
            if ((gatherWindow == null) || (!gatherWindow.isVisible()))
            {
                // Create new window with gather view
   				//gatherWindow = currentModel.copyViewPreferred(currentModel.win3d);
   				//gatherWindow.glPanel3d.checkAddView(StsViewGatherMicroseismic.class);
   				StsToolbar[] toolbars = gatherWindow.getToolbars();
   				for(int i=0; i<toolbars.length; i++)
   				{
   					if(!toolbars[i].getName().equalsIgnoreCase(StsViewSelectToolbar.NAME))
   						gatherWindow.removeToolbar(toolbars[i].getName());
   				}
   			}
   			//StsViewGatherMicroseismic view = (StsViewGatherMicroseismic)gatherWindow.glPanel3d.getCurrentView();
            sensor.displayGather(eventIdx, view);

            // HACK to demonstrate the capability, really want to jump to a location in a single very long gather.
            //gatherWindow.adjustCursorXY(xy[0], xy[1]);
            //((StsViewGather3d)gatherWindow.glPanel3d.getCurrentView()).viewObjectRepaint(sensor);
            gatherWindow.repaint();
        }
        */
        if(!((StsMouseActionToolbar)currentModel.win3d.getToolbarNamed(StsMouseActionToolbar.NAME)).isReadoutSelected())
            win3d.adjustCursorXY(sensor.getXYZ(eventIdx)[0], sensor.getXYZ(eventIdx)[1]);
    }

   public void setEnable(boolean val)
   {
       if(enable == val) return;
       enable = val;
       Iterator iter = getObjectIterator();
       while(iter.hasNext())
       {
           StsDynamicSensor sensor = (StsDynamicSensor)iter.next();
           sensor.setIsVisible(val);
       }
   }


   public boolean getEnable() { return enable; }
}