
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Boundary;

import com.Sts.Actions.*;
import com.Sts.Actions.Wizards.Horizons.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;

// Here we just want to install a toolbar and then terminate the action
public class StsBuildBoundary extends StsAction
{
	StsBuildBoundaryToolbar buildBoundaryToolbar;

 	public StsBuildBoundary(StsActionManager actionManager)
    {
        super(actionManager);
        model.setActionStatus(StsHorizonsWizard.class.getName(),StsModel.ENDED);
        buildBoundaryToolbar = (StsBuildBoundaryToolbar)model.win3d.getToolbarNamed(StsBuildBoundaryToolbar.NAME);
        if(buildBoundaryToolbar != null)
        {
            buildBoundaryToolbar.setVisible(true);
            buildBoundaryToolbar.initialize(actionManager);
            buildBoundaryToolbar.revalidate();
            return;
        }
        buildBoundaryToolbar = new StsBuildBoundaryToolbar(actionManager);
        model.win3d.addToolbar(buildBoundaryToolbar);
    }

	public boolean start() { return true; }

	public void checkAddToolbar()
	{
		model.win3d.addToolbar(buildBoundaryToolbar);
    }

    public boolean checkStartAction()
    {

        if(!model.hasObjectsOfType(StsSection.class, StsParameters.BOUNDARY)) return true;

        boolean yes = StsYesNoDialog.questionValue(model.win3d,"Boundary already built.\n  Do you wish to delete boundary as well as fault framework and model (if they exist)?.");
        if(yes)
        {
            model.deleteObjectsOfType(StsSection.class, StsParameters.FAULT);
            model.deleteObjectsOfType(StsSection.class, StsParameters.BOUNDARY);
            StsBuiltModel builtModel = (StsBuiltModel)model.getCurrentObject(StsBuiltModel.class);
            if(builtModel != null) builtModel.delete();
            model.setActionStatus(StsBuildBoundary.class.getName(), StsModel.CAN_START);
        }
        return yes;
    }

    public void clearToolbar()
    {
        model.win3d.closeToolbar(buildBoundaryToolbar);
    }

    public boolean end()
    {
        model.win3d.removeToolbar(StsBuildBoundaryToolbar.NAME);
        return true;
    }
}
