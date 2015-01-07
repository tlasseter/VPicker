
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Boundary;

import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Toolbars.*;

public class StsAutoBoundary extends StsPolygonBoundary
{

	public StsAutoBoundary(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean start()
    {
		if(!actionOK())  return false;
        statusArea.setTitle("Automatically builds boundary around complete grid.");
        actionManager.endCurrentAction();
		return true;
    }

    private boolean constructCornerGridPoints()
    {
    	if( surface == null ) return false;
		gridPoints = new StsGridPoint[4];
        int ix = surface.getNCols() - 1;
        int iy = surface.getNRows() - 1;
		gridPoints[0] = surface.getGridPoint(0,0);
		gridPoints[1] = surface.getGridPoint(iy,0);
		gridPoints[2] = surface.getGridPoint(iy,ix);
		gridPoints[3] = surface.getGridPoint(0,ix);
		return true;
    }

    public boolean end()
    {
        boolean success = constructBoundary();
        if (success)
		{
			logMessage("Auto boundary successfully created.");
            model.win3d.removeToolbar(StsBuildBoundaryToolbar.NAME);
            model.win3d.win3dDisplay();
            model.setActionStatus(StsBuildBoundary.class.getName(), StsModel.ENDED);
        }
		else
			logMessage("Failed to construct rectangular boundary around model.");
        return success;
    }

	private boolean constructBoundary()
	{
		isLoop = true;
        if(!constructCornerGridPoints()) return false;
		if(!constructBoundarySections()) return false;
		return true;
	}
}
