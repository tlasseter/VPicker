package com.Sts.Actions.Wizards.WellPlan;

import com.Sts.Actions.Wizards.PlatformPlan.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsDefinePlatform extends StsWizardStep
{
    StsConfigurePlatformPanel panel;
    StsHeaderPanel header;

    StsPlatform platform = null;
    StsCursor3d cursor3d;
    StsDefineWellPlanWizard wizard = null;

    public StsDefinePlatform(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsDefineWellPlanWizard) wizard;
        panel = new StsConfigurePlatformPanel(wizard);
        header = new StsHeaderPanel();
        setPanels(panel, header);
 //       panel.setPreferredSize(new Dimension(500, 700));
        header.setTitle("Define New Platform Plan");
        header.setSubtitle("Set parameters and location interactively");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#WellPlan");
        header.setInfoText(wizardDialog, "Specify the name of the platform and set configuration parameters.\n");
    }

    public boolean start()
    {
        platform = wizard.getPlatform();
        panel.initialize(platform);
        wizard.rebuild();
        model.win3d.cursorPickSetup();
        cursor3d = model.win3d.getCursor3d();
        return true;
    }

    public void cancel()
    {
        wizard.setPlatform(null);
    }

    public boolean end()
    {
        if(!wizard.isCanceled())
        {
            if (!platform.isOriginSet())
            {
                new StsMessage(model.win3d, StsMessage.WARNING,
                    "Platform x and y  must be defined in text boxes or interactively picked.");
                return false;
            }
            platform.configurePlatform();
            if(!platform.isPersistent())platform.addToModel();
            return true;
        }
        return false;
    }

    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        StsSurfaceVertex vertex;
        StsPoint pickPoint = null;
        StsCursorPoint cursorPt;

        int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
        if (leftButtonState == StsMouse.NONE)
            return true;

        int currentDirNo = cursor3d.getCurrentDirNo();
        if (currentDirNo == StsParameters.NO_MATCH) // do we have a 3-d cursor?
        {
            logMessage("Use the 3-D Cursor Tool to define an active CURSOR plane");
            return true;
        }
        // Determine the closest intersecting cursor plane and make it current
        StsView currentView = ((StsGLPanel3d)glPanel).getView();
        if (currentView instanceof StsView3d)
        {
            StsGLPanel3d glPanel3d = (StsGLPanel3d)glPanel;
            cursorPt = cursor3d.getCursorPoint(glPanel3d, mouse);
            if(cursorPt == null)
             {
                 new StsMessage(model.win3d, StsMessage.WARNING, "Couldn't get pick on cursor.");
                 return true;
            }
            cursor3d.setCurrentDirNo(cursorPt.dirNo);
        }
        else if (currentView instanceof StsViewCursor)
            cursorPt = ( (StsViewCursor) currentView).getCursorPoint(mouse);
        else
            return true;

        pickPoint = cursorPt.point;

        // display x-y-z location
        if (pickPoint != null) logMessage("X: " + pickPoint.v[0] + " Y: " + pickPoint.v[1] + " Z: " + pickPoint.v[2]);

        if (pickPoint == null)return true;

        if (StsMouse.isButtonStateReleasedOrClicked(leftButtonState)) // add a point
        {
            double[] xy = model.getProject().getAbsoluteXYCoordinates(pickPoint);
            platform.setXOrigin(xy[0]);
            platform.setYOrigin(xy[1]);
            panel.setXOrigin(xy[0]);
            panel.setYOrigin(xy[1]);
        }
        return true;
    }
}
