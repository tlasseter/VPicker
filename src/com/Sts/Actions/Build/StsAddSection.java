
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Build;

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;
import com.Sts.UI.*;

public class StsAddSection extends StsAction
{
    private StsSection section;
    private StsFaultLine[] pickedLines = new StsFaultLine[2];
    private int nPickedLines = 0;

 	public StsAddSection(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public boolean start()
    {
        try
        {
            statusArea.setTitle("Build Fault Section: ");

            int nFaultLines = model.getStsClassSize(StsLine.class);
            if (nFaultLines < 2)
            {
                statusArea.textOnly();
                logMessage("At least two fault lines are needed.");
                return false;
            }

//            model.viewProperties.save();
//            model.viewProperties.set("Sections", true);
//            model.viewProperties.set("Patches", false);

//            model.viewProperties.set("EditingSections", true);

//            win3d.glPanel3d.set3dOverlay(true);
            model.win3d.win3dDisplay();

            logMessage("Select two lines in 3D window.");
            addAbortButton();
        }
        catch(Exception e)
        {
            StsException.outputException("StsAddSection.start() failed.",
                e, StsException.WARNING);
            return false;
        }
        return true;
    }

    public boolean end()
    {
        boolean ok = false;
		try
		{
			if(section != null)
			{
				ok = section.checkSection(model.win3d);
				if(ok)
				{
                    section.completeSection();
                    // section.constructSection();
					//section.constructRibbonSectionEdges();
				}
			}
		}
		catch(Exception e)
		{
			StsException.outputException("StsAddSection.end() failed.",
				e, StsException.WARNING);
			return false;
		}
		finally
		{
//			model.viewProperties.restore();

			StsFaultLine.clearHighlightedLines();
//			win3d.glPanel3d.set3dOverlay(false);
			model.win3d.win3dDisplay();
			statusArea.textOnly();
			return ok;
		}
    }

    /** mouse action for 3d window */
   	public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
	{
        StsFaultLine line = null;
		try
		{
	    	int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

	      	if(leftButtonState != StsMouse.RELEASED) return true;

			// leftButton released: add a point
	        StsObject[] visibleLines = model.getVisibleObjectList(StsFaultLine.class);

			line = (StsFaultLine)StsJOGLPick.pickClass3d( glPanel, visibleLines, StsJOGLPick.PICKSIZE_MEDIUM, StsJOGLPick.PICK_CLOSEST );
	        if(line == null) return true;

            logMessage("Picked " + line.getLabel() + "[" + line.getIndex() + "]");
            pickedLines[nPickedLines++] = line;
            StsFaultLine.highlightedList.add(line);

	        if(nPickedLines == 1)
            {
	            model.win3d.win3dDisplay();
                return true;
            }

            StsFaultLine[] faultSticks = getFaultSticks(pickedLines);

            section = StsSection.constructor(StsParameters.FAULT, faultSticks);
			if(section == null)
			{
				logMessage("Add Section constructor failed: see error log.");
				return false;
			}
            
	        model.win3d.win3dDisplay();
            glPanel.actionManager.endCurrentAction();
			return true;
	    }
		catch(Exception e)
		{
			StsException.outputException("StsAddSection.performMouseAction() failed.",
				e, StsException.WARNING);
			return false;
		}
	}

    StsFaultLine[] getFaultSticks(StsFaultLine[] pickedLines)
    {
        if(pickedLines == null || pickedLines.length < 2) return null;
        StsFaultStickSet stickSet0 = pickedLines[0].faultStickSet;
        StsFaultStickSet stickSet1 = pickedLines[1].faultStickSet;
        if(stickSet0 == null || stickSet1 == null) return pickedLines;
        if(stickSet0 != stickSet1)
        {
            boolean connect = StsYesNoDialog.questionValue(model.win3d,  "Picked sticks belong to two different sets.  Connect anyways?");
            if(connect) return pickedLines;
            else        return null;
        }
        int index0 = stickSet0.getFaultStickIndex(pickedLines[0]);
        int index1 = stickSet0.getFaultStickIndex(pickedLines[1]);
        if(index0 > index1)
        {
            int temp = index0;
            index0 = index1;
            index1 = temp;
        }
        int nSticks = index1 - index0 + 1;
        if(nSticks == 2)
            return pickedLines;

        Object[] stickObjects = stickSet0.faultSticks.getList().getList();
        StsFaultLine[] sticks = new StsFaultLine[nSticks];
        for(int n = 0; n < nSticks; n++)
            sticks[n] = (StsFaultLine)stickObjects[n+index0];
        return sticks;
    }
}

