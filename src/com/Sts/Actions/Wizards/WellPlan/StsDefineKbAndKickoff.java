package com.Sts.Actions.Wizards.WellPlan;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsDefineKbAndKickoff extends StsWizardStep
{
    StsDefineKbAndKickoffPanel panel;
    StsHeaderPanel header;

    StsDefineWellPlanWizard wizard = null;

    StsWellPlan plan;
    StsWellPlan wellPlan = null;
    StsCursor3d cursor3d;
    StsSeismicVelocityModel velocityModel;

    StsColor stsColor;
    String name = new String("WellPlan");
    Color color = Color.RED;
    boolean visible = false;

    static int suffix = 1;

    public StsDefineKbAndKickoff(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsDefineWellPlanWizard)wizard;
        panel = new StsDefineKbAndKickoffPanel(wizard);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 500));
        header.setTitle("Planned Well Trajectory Definition");
        header.setSubtitle("Define/pick trajectory and kickoff location");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#WellPlan");
        header.setInfoText(wizardDialog,
                           "(1) Specify sealevel depth of the of KB.\n" +
                           "(2) Specify sealevel depth of kickoff.\n" +
                           "(3) Press the Next> Button to move to target specification.\n");
    }

    public boolean start()
    {
        // Need to get paths only related to current Plan
        wellPlan = wizard.getWellPlan();
        panel.initialize(wellPlan);
        wizard.rebuild();
        cursor3d = model.win3d.getCursor3d();
        return true;
    }

    public void cancel()
    {
        wellPlan.deleteDrawPoints();
    }

    public boolean end()
    {
 //       planSet = wizard.getWellPlanSet();
        return true;
    }

    /** mouse action for 3d window */
    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        StsPoint pickPoint;
        StsCursorPoint cursorPt;

        int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
        if(leftButtonState == StsMouse.NONE) return true;
        if(StsMouse.isButtonStateReleasedOrClicked(leftButtonState)) // Give them an info message
        {
            if(wellPlan.getPlatform() != null)
            {
                new StsMessage(model.win3d, StsMessage.INFO, "Must specify KB and kickoff depths in wizard text boxes.");
            }
            else // single well, so we can pick the location interactively
            {
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
                    cursor3d.setCurrentDirNo(cursorPt.dirNo);
                }
                else if (currentView instanceof StsViewCursor)
                    cursorPt = ( (StsViewCursor) currentView).getCursorPoint(mouse);
                else
                    return true;

                pickPoint = cursorPt.point;

                // display x-y-z location
                if(pickPoint != null) logMessage("X: " + pickPoint.v[0] + " Y: " + pickPoint.v[1] + " Z: " + pickPoint.v[2]);

                if (pickPoint == null) return true;

                double[] xy = model.getProject().getAbsoluteXYCoordinates(pickPoint);
                double xOrigin = xy[0];
                double yOrigin = xy[1];
                wellPlan.setXOrigin(xOrigin);
                panel.xOriginBean.setValue(xOrigin);
                wellPlan.setYOrigin(yOrigin);
                panel.yOriginBean.setValue(yOrigin);
            }
        }
        model.win3d.win3dDisplay();
        return true;
    }

    public StsWellPlan getWellPlanSet() { return wellPlan; }
}
