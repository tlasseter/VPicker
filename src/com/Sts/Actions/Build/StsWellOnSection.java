
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
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

public class StsWellOnSection extends StsAction
{
    StsFaultLine line;
    StsCursor3d cursor3d;
    StsSection currentSection;

    /** construct a action handler for building a line on a section surface. */
	public StsWellOnSection(StsActionManager actionManager)
    {
        super(actionManager);
		addEndAllButton();
		addAbortButton();
		addEndButton();
	}

    public boolean start()
    {
		try
		{
		    currentSection = StsSection.getCurrentSection();
			if(currentSection  == null)
			{
				logMessage("No Current Section: pick an edge on a section to define one.");
				return false;
			}

			if(!currentSection.checkConstructSection()) return false;

            line = StsFaultLine.buildFault();
            line.setDrawVertices(true);
            line.setOnSection(currentSection);
	        statusArea.setTitle("Build " + line.getLabel() + ":");
            logMessage("Select a fault section and pick vertices, then press Done button.");
        }
        catch(Exception e)
        {
         	StsException.outputException(e, StsException.WARNING);
            if (line != null) line.setDrawVertices(false);
            return false;
        }
		return true;
    }

    /** mouse action for 3d window */
   	public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
	{
 		StsSurfaceVertex vertex;
        StsPoint pickPoint;

    	int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

        StsMousePoint mousePoint = mouse.getMousePoint();
        pickPoint = currentSection.getPointOnSection((StsGLPanel3d)glPanel, mousePoint);

        if(pickPoint == null) return true;

       	if(leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
        {
            // display x-y-z location
			logMessage("X: " + pickPoint.v[0] + " Y: "
            						+ pickPoint.v[1] + " Z: " + pickPoint.v[2]);
        }

       	if(StsMouse.isButtonStateReleasedOrClicked(leftButtonState))  // add a point
       	{

            vertex = line.addLineVertex(pickPoint, true, false);
            if(vertex != null)
            {
                logMessage("Added vertex: " + vertex.getIndex() +
                        " to current " + line.getLabel() + ".");
                model.win3d.win3dDisplay();
            }
        }
        return true;
	}

    /** add a line when the mouse action finishes */
  	public boolean end()
    {
        if (line != null) line.setDrawVertices(false);
        statusArea.textOnly(); // removes end/abort/endAll buttons
        logMessage(line.getLabel() + " built successfully.");
        model.win3d.win3dDisplay();
        return true;
    }
}
