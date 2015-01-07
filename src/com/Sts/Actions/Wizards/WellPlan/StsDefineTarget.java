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

public class StsDefineTarget extends StsWizardStep
{
    StsDefineTargetPanel panel;
    StsHeaderPanel header;

    StsDefineWellPlanWizard wizard;
    StsColor stsColor;
    StsCursor3d cursor3d;

    String exportName = null;
	StsWellPlanSet wellPlanSet = null;
	StsWellPlan wellPlan = null;

    public StsDefineTarget(StsWizard wizard)
    {
        super(wizard);
        this.wizard = (StsDefineWellPlanWizard)wizard;
        panel = new StsDefineTargetPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 500));
        header.setTitle("Planned Well Definition");
        header.setSubtitle("Pick the target");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#WellPlan");
        header.setInfoText(wizardDialog,
                           "(1) Select the operation mode, conduct graphical operations.\n" +
                           " **** Pick - graphically add points from the bottom of the well.\n" +
                           " **** Insert - graphically insert points between existing vertices.\n" +
                           " **** Delete - remove graphically selected vertices.\n" +
                           " **** Move - graphically pick-up and move existing vertices.\n\n" +
                           " **** As a vertex is defined it will be put in the list and it's \n" +
                           "      x,y,z and rate of change will be isVisible and editable.\n" +
                           "(1a) Alternatively, press the Import Button to read in a UT pln file.\n" +
                           "(2) Press the Finish> Button to complete the well definition.\n");
    }

    public boolean start()
    {
        wellPlanSet = wizard.getWellPlanSet();
        if(wellPlanSet == null) return false;

	    if(model.getProject().isDepth())
		{
			StsMessage.printMessage("Model is currently in depth; you must switch to time for target and mid-point picking.");
		}

		wellPlan = wellPlanSet.getCurrentWellPlan();
        panel.initialize(wellPlan);
        wizard.rebuild();
        wellPlan.adjustKBPoint();
        try
        {
            model.win3d.cursorPickSetup();
            cursor3d = model.win3d.getCursor3d();
            statusArea.setTitle("Build " + wellPlanSet.getName() + ":");
            logMessage("Pick target point on 3-D cursor plane, then press Done button.");
        }
        catch(Exception e)
        {
//            if (plannedWell != null) plannedWell.setDrawVertices(false);
            StsException.outputException(e, StsException.WARNING);
            return false;
        }
        return true;
    }

    /** mouse action for 3d window */
    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
        StsSurfaceVertex vertex;
        StsPoint pickPoint = null;
        StsCursorPoint cursorPt;

		if(model.getProject().isDepth())
		{
			StsMessage.printMessage("Model is currently in depth; you must switch to time for target and mid-point picking.");
			return true;
		}
        int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
        if(leftButtonState == StsMouse.NONE)
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
        if(pickPoint != null) logMessage("X: " + pickPoint.v[0] + " Y: " + pickPoint.v[1] + " Z: " + pickPoint.v[2]);

        if (pickPoint == null) return true;

        if(StsMouse.isButtonStateReleasedOrClicked(leftButtonState)) // add a point
        {
            StsPoint timeDepthPoint = wellPlan.computeTimeDepthPoint(pickPoint); // computed in rotated coordinates
            double[] xy = model.getProject().getAbsoluteXYCoordinates(timeDepthPoint);
            panel.setTargetX(xy[0]);
            panel.setTargetY(xy[1]);
            panel.setTargetT(timeDepthPoint.getT());
            panel.setTargetZ(timeDepthPoint.getZ());
            wellPlan.addTargetPoint(timeDepthPoint);
            wizard.enableFinish();
        }
/*
        else if (mode == DELETE)
        {
            if (leftButtonState != StsMouse.RELEASED)
                return true;

            mouse.clearButtonState(StsMouse.LEFT);

            StsMethod pickMethod = new StsMethod(plannedWell, "pickVertices", model.glPanel3d);
            StsMethodPick picker = new StsMethodPick(model.glPanel3d, pickMethod, StsMethodPick.PICKSIZE_LARGE, StsMethodPick.PICK_ALL);
            try
            {
                if(model.glPanel3d.glc.pick(picker))
                {
                    StsPickItem pickItem = picker.pickItems[0];
                    int vertexIndex = pickItem.names[0];
                    StsSurfaceVertex deleteVertex = (StsSurfaceVertex)plannedWell.getLineVertices().getElement(vertexIndex);
                    wellPlanSet.deleteLineVertex(deleteVertex);
                }
            }
            catch(Exception e)
            {
                StsException.systemError("StsDefineVertices.performMouseAction() vertex delete failed.");
            }
        }
*/
//        panel.adjustVertices(plannedWell);
        model.win3d.win3dDisplay();

        return true;
    }
/*
    public void exportPlannedWell()
    {
        plannedWell.export(StsParameters.TD_TIME_DEPTH_STRING);
    }
*/
    public void setExportName(String exportName) { this.exportName = exportName; }
    public String getExportName() { return exportName; }

    public boolean end()
    {
        statusArea.textOnly();
 //       plannedWell.setDrawVertices(false);
        wellPlanSet.setDisplayPrevious(false);
/*
        if(plannedWell.getPoints() != null)
            logMessage("Added " + plannedWell.getPoints().length + " vertices.");
        try { Thread.sleep(1000); }
        catch(Exception e) { }
*/
        logMessage(wellPlanSet.getName() + " built successfully.");
        model.win3d.win3dDisplay();
        return true;
    }

}
