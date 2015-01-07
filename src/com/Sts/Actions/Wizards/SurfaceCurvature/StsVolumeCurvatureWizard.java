
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.SurfaceCurvature;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.event.*;

public class StsVolumeCurvatureWizard extends StsWizard implements ActionListener
{
    public StsSeismicVolume seismicVolume = null;
    //public StsPatchVolume patchVolume = null;
    public StsSelectVolume selectVolume;
    public StsPatchPick patchPick;
    public StsProcessVolume processVolume;
    public StsPatchVolume patchVolume = null;

    public StsVolumeCurvatureWizard(StsActionManager actionManager)
    {
        super(actionManager, 500, 650);
        addSteps
        (
            new StsWizardStep[]
            {
            	selectVolume = new StsSelectVolume(this),
            	patchPick = new StsPatchPick(this),
            	processVolume = new StsProcessVolume(this)
            	// patchInterp = new StsInterpolatePatches(this)
            }
        );
    }
    
    public boolean start()
    {
        disableFinish();
        dialog.addActionListener(this);
        System.runFinalization();
        System.gc();
        dialog.setTitle("Volume Curvature Analysis");
        return super.start();
    }

    public void actionPerformed(ActionEvent e)
    {
         if(e.getID() == WindowEvent.WINDOW_CLOSING)
             cancel();
    }

    public void setSelectedVolume(StsSeismicVolume volume)
    {
        if(volume == null) return;
        if(volume.fileMapRowFloatBlocks == null)
        {
            String rowFloatFilename = volume.createFloatRowVolumeFilename(volume.stemname);
            new StsMessage(model.win3d, StsMessage.ERROR, "float volume " + rowFloatFilename + " unavailable.  Can't build volume.");
            return;
        }
        seismicVolume = volume;
    }

    public StsSeismicVolume getSelectedVolume()
    {
          return seismicVolume;
    }
    
    public StsPatchVolume getPatchVolume()
    {
    	if (patchVolume == null)
    	{
    		patchVolume = new StsPatchVolume(seismicVolume);
    	}
        return patchVolume;
    }
    
    public void setPatchVolume(StsPatchVolume volume)
    {
    	patchVolume = volume;
    }

    public void setCurvatureMin(float min)
    {
        if(patchVolume == null) return;
        patchVolume.setDataMin(min);
        processVolume.panel.updateColorscale();
    }

    public float getCurvatureMin()
    {
        if(patchVolume == null) return -1.0f;
        else return patchVolume.dataMin;
    }

    public void setCurvatureMax(float max)
    {
        if(patchVolume == null) return;
        patchVolume.setDataMax(max);
        processVolume.panel.updateColorscale();
    }

    public float getCurvatureMax()
    {
        if(patchVolume == null) return 1.0f;
        else return patchVolume.dataMax;
    }

    public boolean end()
    {
         return super.end();
    }

    public void previous()
    {
         gotoPreviousStep();
    }

    public void next()
    {

        if(currentStep == selectVolume && selectVolume.panel.patchVolSelected)
            gotoStep(processVolume);
        else
            gotoNextStep();
    }

    public void finish()
    {
         super.finish();
    }

    public void cancel()
    {
        super.cancel();
    }

    public boolean exportView()
    {
          return true;
    }

    public boolean checkPrerequisites()
    {
        StsObject[] volumeObjects = model.getObjectList(StsSeismicVolume.class);
        if(volumeObjects.length == 0)
        {
             reasonForFailure = "No volumes available. Must have at least one seismic volume loaded.";
             return false;
        }
        return true;
    }

    public boolean patchPick(StsMouse mouse, StsGLPanel glPanel)
    {
        try
        {
            StsGLPanel3d glPanel3d = (StsGLPanel3d)glPanel;
            int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
            if(leftButtonState != StsMouse.RELEASED) return true;
            StsCursorPoint cursorPoint = model.win3d.getCursor3d().getCursorPoint(glPanel3d, mouse);
            if(cursorPoint == null)
            {
                StsMessageFiles.errorMessage("Failed to get pick point on cursor section.");
                return false;
            }
            patchVolume.addRemoveSelectedPatch(cursorPoint);
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsEdgeOnSurface.performMouseAction() failed.",
                    e, StsException.WARNING);
            return false;
        }
    }

    static public void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsVolumeCurvatureWizard curvatureWizard = new StsVolumeCurvatureWizard(actionManager);
        curvatureWizard.start();
    }
}