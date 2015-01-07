package com.Sts.Actions.Wizards.WellPlan;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2005</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version c62e
 */

public class StsDefineActualWellStart extends StsWizardStep
{
	StsDefineActualWellStartPanel panel;
	StsHeaderPanel header;

	StsDefineWellPlanWizard wizard = null;

	StsWellPlan wellPlan = null;
	boolean isNewPlanSet;
	StsCursor3d cursor3d;
	StsSeismicVelocityModel velocityModel;

	StsColor stsColor;
	String name = new String("WellPlan");
	Color color = Color.RED;
	boolean visible = false;

	static int suffix = 1;

	public StsDefineActualWellStart(StsWizard wizard)
	{
		super(wizard);
		this.wizard = (StsDefineWellPlanWizard)wizard;
		StsWell drillingWell = this.wizard.getDrillingWell();
		panel = new StsDefineActualWellStartPanel(wizard, drillingWell);
		header = new StsHeaderPanel();
		setPanels(panel, header);
		panel.setPreferredSize(new Dimension(500, 500));
		header.setTitle("Select Drilling Well End or Sidetrack Point");
		header.setSubtitle("Accept or interactively select point.");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#WellPlan");
		header.setInfoText(wizardDialog,
						   "(1) Accept, define, or select point.\n" +
						   "(2) Accept define turn rate.\n" +
						   "(3) Press the Next> Button to move to target specification.\n");
	}

	public boolean start()
	{
		// Need to get paths only related to current Plan
		wellPlan = wizard.getWellPlan();
		isNewPlanSet = wizard.getIsNewPlanSet();
		panel.initialize(wellPlan);
		wizard.rebuild();
		cursor3d = model.win3d.getCursor3d();
		wellPlan.constructDrillingPlan(isNewPlanSet);
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
		StsSurfaceVertex vertex;
		StsPoint pickPoint = null;
		StsCursorPoint cursorPt;

		int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
		if(leftButtonState == StsMouse.NONE) return true;
		if(StsMouse.isButtonStateReleasedOrClicked(leftButtonState)) // Give them an info message
		{
			new StsMessage(model.win3d, StsMessage.INFO, "Must specify KB and kickoff in wizard text boxes.");
		}
/*
		int currentDirNo = cursor3d.getCurrentDirNo();
		if (currentDirNo == StsParameters.NO_MATCH) // do we have a 3-d cursor?
		{
			logMessage("Use the 3-D Cursor Tool to define an active CURSOR plane");
			return true;
		}
		// Determine the closest intersecting cursor plane and make it current
		StsView currentView = model.win3d.glPanel3d.getCurrentView();
		if (currentView instanceof StsView3d)
		{
			cursorPt = cursor3d.getCursorPoint(mouse);
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
			mouse.clearButtonState(StsMouse.LEFT);
			float zKickoff = pickPoint.getZ();
			planSet.setZKickoff(zKickoff);
			panel.zKickoffBean.setValue(zKickoff);
			planSet.initializeKickoffLeg();
 */
		/*
			double[] xy = model.getProject().getAbsoluteXYCoordinates(pickPoint);
			double xOrigin = xy[0];
			double yOrigin = xy[1];
			planSet.setXOrigin(xOrigin);
			panel.xOriginBean.setValue(xOrigin);
			planSet.setYOrigin(yOrigin);
			panel.yOriginBean.setValue(yOrigin);
		*/
 //       }

		model.win3d.win3dDisplay();
		return true;
	}

	public StsWellPlan getWellPlan() { return wellPlan; }
}
