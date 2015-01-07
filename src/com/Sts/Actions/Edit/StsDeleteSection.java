
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

public class StsDeleteSection extends StsAction
{
    StsSectionClass sectionClass = null;

    public StsDeleteSection(StsActionManager actionManager)
    {
        super(actionManager);
    }

    public boolean start()
    {
        if (model.getNObjects(StsSection.class) < 1)
        {
            logMessage("No faults sections to delete."
                    + "  Terminating action.");
            return false;
        }

        statusArea.setTitle("Delete fault section:");
        logMessage("Select a fault section in the 3D window.");
        addAbortButton();

        sectionClass = (StsSectionClass)model.getCreateStsClass(StsSection.class);
        return true;
    }

    public boolean end()
    {
    	logMessage(" ");
		statusArea.textOnly();
        return true;
    }

    /** mouse action for 3d window */
   	public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
	{
        StsSectionEdge edge = null;
    	int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

       	if(StsMouse.isButtonStateReleasedOrClicked(leftButtonState))
       	{
        	try
            {
				edge = (StsSectionEdge)StsJOGLPick.pickClass3d( glPanel, sectionClass.getVisibleSectionEdges(), StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_CLOSEST);
                if(edge != null)
                {
                    StsSection edgeSection = edge.getSection();
                    if(edgeSection != null)
                        edgeSection.delete();
                    else
                        edge.delete();

                    model.win3d.win3dDisplay();

					glPanel.actionManager.endCurrentAction();
                }
                else
                     logMessage("No fault section edge selected: try again.");
            }
            catch (Exception e)
            {
            	logMessage("Failed to selected a fault section: try again.");
            	return false;
            }
        }
        return true;
	}

}
