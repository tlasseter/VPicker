package com.Sts.UI.Toolbars;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Mar 8, 2008
 * Time: 2:59:08 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsCenterObject extends StsAction
{
    public StsCenterObject(StsActionManager actionManager)
    {
        super(actionManager, true);
    }

    public boolean start()
    {
        return true;
    }

    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        StsGLPanel3d glPanel3d = (StsGLPanel3d)glPanel;
        StsView view = glPanel3d.getView();
        if(!(view instanceof StsView3d)) return true;

        int mouseButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
        if(mouseButtonState != StsMouse.RELEASED) return true;
        StsView3d view3d = (StsView3d)view;
        float[] parms = view3d.getViewParameters();
        StsPoint point = getPointOnObject(glPanel3d, mouse);
        if(point == null) return true;
        view3d.setViewCenter(point);
        view3d.setViewParameters(parms);
	    view3d.repaint();
        view3d.moveLockedWindows();
        glPanel3d.fireViewChangeEvent();
        actionManager.endCurrentAction();
        return true;
    }

    private StsPoint getPointOnObject(StsGLPanel3d glPanel3d, StsMouse mouse)
    {
        StsWell well = (StsWell) StsJOGLPick.pickVisibleClass3d(glPanel3d, StsWell.class, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_CLOSEST);
        if (well != null)
        {
            StsPickItem item = StsJOGLPick.pickItems[0];
		    int nSegment = item.names[1];
            StsPoint xyzmt = well.getPointOnLineNearestMouse(nSegment, mouse, glPanel3d);
		    return xyzmt.getXYZorTPoint();
        }
        StsLiveWell lwell = (StsLiveWell) StsJOGLPick.pickVisibleClass3d(glPanel3d, StsLiveWell.class, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_CLOSEST);
        if (lwell != null)
        {
            StsPickItem item = StsJOGLPick.pickItems[0];
		    int nSegment = item.names[1];
            StsPoint xyzmt = lwell.getPointOnLineNearestMouse(nSegment, mouse, glPanel3d);
		    return xyzmt.getXYZorTPoint();
        }
        StsSensor sensor = (StsDynamicSensor) StsJOGLPick.pickVisibleClass3d(glPanel3d, StsDynamicSensor.class, StsJOGLPick.PICKSIZE_LARGE, StsJOGLPick.PICK_FIRST);
        if (sensor != null)
        {
                 StsPickItem items = StsJOGLPick.pickItems[0];
                 if(items.names.length < 2) return null;
                 int sensorIdx = items.names[1];
                 if(sensor instanceof StsDynamicSensor)
				    return new StsPoint(((StsDynamicSensor)sensor).getXYZ(sensorIdx));
                 else
                    return null;
        }
        return null;
    }

    public boolean end()
    {
        return true;
    }
}
