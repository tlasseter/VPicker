package com.Sts.Actions.Wizards.Color;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.Actions.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;

public class StsColorscaleAction extends StsAction implements Runnable
{
	StsColorscale originalColorscale;
	StsColorscaleDialog colorscaleDialog = null;
//    StsColorscale colorscale = null;
	StsColorscalePanel observerColorscalePanel = null;
	Frame frame = null;
//    float[] dataDist = null;
	boolean success = false;

	public StsColorscaleAction(StsActionManager actionManager, StsColorscalePanel observerColorscalePanel)
	{
		super(actionManager, true);
		// frame = actionManager.getModel().win3d;
		this.observerColorscalePanel = observerColorscalePanel;
		originalColorscale = new StsColorscale(false);
		originalColorscale.copySettingsFrom(observerColorscalePanel.getColorscale());
//		this.colorscale = observerColorscalePanel.getColorscale();
//		dataDist = observerColorscalePanel.getHistogram();
	}

	public StsColorscaleAction(StsActionManager actionManager, StsColorscalePanel observerColorscalePanel, StsWin3d win3d)
	{
		super(actionManager, true);
		this.observerColorscalePanel = observerColorscalePanel;
		this.frame = win3d;
//		this.colorscale = observerColorscalePanel.getColorscale();
//		dataDist = observerColorscalePanel.getHistogram();
	}

	public boolean start()
	{
		run();
		return true;
	}

	public void run()
	{
		try
		{
            StsToolkit.runLaterOnEventThread
            (
                new Runnable()
                {
                    public void run()
                    {
                        if(observerColorscalePanel == null)
                        {
                            StsMessageFiles.infoMessage("Need to define object first.");
                            return;
                        }
                        if(model == null)
                            colorscaleDialog = new StsColorscaleDialog("Color Spectrum Editor", observerColorscalePanel, frame);
                        else
                            colorscaleDialog = new StsColorscaleDialog("Color Spectrum Editor", model, observerColorscalePanel, frame);

                        int yPos = colorscaleDialog.getY() - 200;
                        if(yPos < 0)
                            yPos = 0;
                        colorscaleDialog.setLocation(colorscaleDialog.getX(), yPos);
                        colorscaleDialog.setVisible(true);
            //            colorscaleDialog.initData(dataDist);
                    }
                }
            );
		}
		catch(Exception e)
		{
			StsException.outputException("StsColorscaleAction.run() failed.", e, StsException.WARNING);
			return;
		}

	}

	// not referenced, not needed TJL 1/31/06
	/*
	 public void setColorscale(StsColorscale cs)
	 {
	  this.colorscale = cs;
	 }
	 */
	public boolean end()
	{
		StsColorscale editedColorscale = observerColorscalePanel.getColorscale();
		success = colorscaleDialog.getSuccess();
		if(success)
		{
			// Issue command for all selected objects.

			StsObject[] objs = model.win3d.objectTreePanel.getSelectedObjects();
			if(objs.length <= 1)
			{
				editedColorscale.commitChanges(originalColorscale);
				editedColorscale.colorsChanged();
//				model.instanceChange(editedColorscale, "colorscale changed to: " + editedColorscale.getSpectrum().getName());
			}
			else
			{
				// If we have multiselected on the objectPanel, we want to update colorscales associated with all objects selected.
				// The only issue is the scaling: do we want to set the ranges to the same.  If not, how do we indicate that?
				// As it is, each of the colorscales other than the first in the selected list will be copied from the first in the list.
				for(int n = 0; n < objs.length; n++)
				{
					updateColorscale(objs[n], editedColorscale);
				}
			}
		}
		else // action has been canceled or otherwise aborted: restore original colorscale
		{
			editedColorscale.copySettingsFrom(originalColorscale);
			editedColorscale.colorsChanged();
		}
		observerColorscalePanel.repaint();
		return success;
	}

	private void updateColorscale(StsObject obj, StsColorscale editedColorscale)
	{
		StsColorscale colorscale = null;
		colorscale = obj.getColorscaleWithName(editedColorscale.getName());
		if(colorscale == null) return;
		colorscale.commitChanges(editedColorscale);
		colorscale.colorsChanged();
	}
}