
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Edit;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

public class StsDeleteLine extends StsAction
{
    public StsDeleteLine(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean start()
    {
        if (model.getCreateStsClass(StsLine.class).getSize() == 0)
        {
            logMessage("No wells, pseudos, or fault lines to delete. Terminating action.");
            return false;
        }
        logMessage("Select a well, pseudo, or fault to delete.");
        addAbortButton();
        return true;
    }

    public boolean end()
    {
        return true;
    }

    /** mouse action for 3d window */
   	public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
	{
        StsLine line = null;

    	int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

       	if(leftButtonState != StsMouse.RELEASED) return true;

		try
        {
            StsMainObject[] visibleLines = model.getVisibleObjectList(StsLine.class);
			if(checkDeleteLine(visibleLines, (StsGLPanel3d)glPanel)) glPanel.actionManager.endCurrentAction();

			StsMainObject[] visibleFaults = model.getVisibleObjectList(StsFaultLine.class);
			if(checkDeleteLine(visibleFaults, (StsGLPanel3d)glPanel)) glPanel.actionManager.endCurrentAction();
			return true;
        }
        catch (Exception e)
        {
			StsException.outputException("Failed during pick.", e, StsException.WARNING);
			return false;
        }
 	}

	private boolean checkDeleteLine(StsMainObject[] visibleLines, StsGLPanel3d glPanel3d)
	{
		StsLine line = null;
		try
		{
			line = (StsLine)StsJOGLPick.pickClass3d(glPanel3d, visibleLines, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_CLOSEST);
			if (line == null) return false;

			if (line.hasConnectedSection())
			{
				statusArea.textOnly();
				logMessage("Cannot delete " + line.getLabel() + ": has connected Section.  Delete section first.");
				return false;
			}
			else
			{
				line.delete();
				model.win3d.win3dDisplay();
				statusArea.textOnly();
				logMessage("Deleted " + line.getLabel() + ".");
				return true;
			}

		}
		catch (Exception e)
		{
			StsException.outputException("StsDeleteLine.checkDeleteLine() failed.", e, StsException.WARNING);
			return false;
		}
	}
}
