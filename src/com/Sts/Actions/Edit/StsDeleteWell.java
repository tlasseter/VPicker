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
import com.Sts.Types.*;
import com.Sts.Utilities.*;

public class StsDeleteWell extends StsAction
{
    public StsDeleteWell(StsActionManager actionManager)
    {
        super(actionManager);
    }

    public boolean start()
    {
        if (model.getCreateStsClass(StsWell.class).getSize() == 0)
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
        StsWell well = null;
        int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

        if (StsMouse.isButtonStateReleasedOrClicked(leftButtonState))
        {
            try
            {
                StsMainObject[] visibleWells = model.getVisibleObjectList(StsWell.class);
                well = (StsWell)StsJOGLPick.pickClass3d(glPanel, visibleWells, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_CLOSEST);
                if (well != null)
                {
                    if (well.hasConnectedSection())
                    {
                        statusArea.textOnly();
                        logMessage("Cannot delete " + well.getLabel() +
                                   ": has connected Section.  Delete section first.");
                    }
                    else
                    {
                        well.delete();
                        model.win3d.win3dDisplay();
                        statusArea.textOnly();
                        logMessage("Deleted " + well.getLabel() + ".");
                    }
                    glPanel.actionManager.endCurrentAction();
                }
                else
                {
                    statusArea.textOnly();
                    logMessage("No well/fault picked.  Try again.");
                }
            }
            catch (Exception e)
            {
                StsException.outputException("Failed during pick.", e, StsException.WARNING);
                return false;
            }
        }
        return true;
    }
}
