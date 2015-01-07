
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
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.event.*;

public class StsFaultOnSurfaces extends StsAction
{
	private StsFaultLine line;
	private StsModelSurface[] surfaces;
	private int nSurfaces;
	private int surfaceIndex = 0;
    private StsModelSurface surface;
    private String surfaceName = null;
	private boolean[] surfacePicked;

   	/** create a new section and attach a mouse listener to it */
 	public StsFaultOnSurfaces(StsActionManager actionManager)
    {
        super(actionManager);
    }

	public void initializeRepeatAction(StsAction lastAction)
	{
		surface = ((StsFaultOnSurfaces)lastAction).getSurface();
	}

	public StsModelSurface getSurface() { return surface; }

	private boolean actionOK()
    {
    	try
        {
            StsSurfaceClass surfaceClass = (StsSurfaceClass)model.getCreateStsClass(StsModelSurface.class);
            if (surfaceClass.getSize() <= 0)
            {
                logMessage("No surfaces have been loaded: terminating action.");
                return false;
            }
	        return true;
        }
        catch(Exception e)
		{
			StsException.outputException("StsFaultOnSurfaces.actionOK() failed.",
				e, StsException.WARNING);
			return false;
		}
    }

	public boolean start()
    {
		if( !actionOK() )  return false;

        try
        {
            statusArea.setTitle("Pick fault line on surfaces: ");

			surfaces = (StsModelSurface[])model.getCastObjectList(StsModelSurface.class);
            if (surfaces == null)
            {
                statusArea.textOnly();
                logMessage("No model surfaces available.");
                return false;
            }

//            model.viewProperties.save();

			nSurfaces = surfaces.length;
			surfacePicked = new boolean[nSurfaces];
			for(int n = 0; n < nSurfaces; n++)
			{
                StsModelSurface surface = surfaces[n];
				surface.setDisplayFill(true);
				surface.setDisplayGrid(true);
				surface.setIsVisible(false);
				surfacePicked[n] = false;
            }

			surfaceIndex = 0;
			setSurface();

            line = StsFaultLine.buildFault();
            line.setDrawVertices(true);

			logMessage("Select a fault line point on surface or 'Next' button to skip.");

		    addAbortButton();
			addEndAllButton();
		    addNextButton();

            model.win3d.win3dDisplay();
            return true;
        }
        catch(Exception e)
		{
			StsException.outputException("StsFaultOnSurfaces.start() failed.",
				e, StsException.WARNING);
			return false;
		}
    }

	public void addNextButton()
    {
        statusArea.addActionButtonListener(new ButtonActionListener(), "Next", true);
	}

    protected class ButtonActionListener implements ActionListener
    {
        public void actionPerformed(ActionEvent e)
        {
			incrementSurface();
        }
    }

	private void incrementSurface()
	{
		for(int n = 0; n < nSurfaces; n++)
		{
			surfaceIndex = (surfaceIndex+1)%nSurfaces;
			if(!surfacePicked[surfaceIndex])
			{
				setSurface();
				return;
			}
		}
		// all surfaces have been picked: end action
		actionManager.endCurrentAction();
		return;
	}

	private void setSurface()
	{
		if(surface != null) surface.setIsVisible(false);
		surface = surfaces[surfaceIndex];
		surface.setIsVisible(true);
		surfaceName = surface.getName();

		if(!surface.checkIsLoaded())
        {
            statusArea.textOnly();
            logMessage("Unable to load " + surfaceName + ".");
		    actionManager.endCurrentAction();
        }
	}

    /** mouse action for 3d window */
   	public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
	{
    	int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

		StsGridPoint gridPoint = surface.getSurfacePosition(mouse, true, (StsGLPanel3d)glPanel);

       	if(StsMouse.isButtonStateReleasedOrClicked(leftButtonState))  // add a point
       	{
			StsPoint point = gridPoint.getPoint();
            StsSurfaceVertex vertex = line.addLineVertex(point, true, false);
			if(vertex == null)
				new StsMessage(model.win3d, StsMessage.ERROR, "Failed to add point to fault on surface: " + surfaceName);
            else
            {
                logMessage("Added vertex: " + vertex.getIndex() + " to current " + line.getLabel() + ".");
				surfacePicked[surfaceIndex] = true;
				incrementSurface();
                model.win3d.win3dDisplay();
            }
        }
        return true;
	}

    /** add a line when the mouse action finishes */
  	public boolean end()
    {
        try
        {
            statusArea.textOnly(); // removes end/abort/endAll buttons

            line.setDrawVertices(false);
            line.extendEnds();
            logMessage(line.getLabel() + " built successfully.");

//            model.viewProperties.restore();
            model.win3d.win3dDisplay();
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsFaultOnSurface.end() failed.",
                e, StsException.WARNING);
            return false;
        }
    }
}


