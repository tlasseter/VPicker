
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

import java.awt.event.*;

public class StsSurfaceCurvatureWizard extends StsWizard implements ActionListener
{
	public StsSurface selectedSurface = null;
    StsSelectSurface selectSurface;
    StsProcessSurface processSurface;

    public StsSurfaceCurvatureWizard(StsActionManager actionManager)
    {
        super(actionManager, 500, 600);
        addSteps
        (
            new StsWizardStep[]
            {
            	selectSurface = new StsSelectSurface(this),
            	processSurface = new StsProcessSurface(this)
            }
        );
    }
    
    public boolean start()
    {
        disableFinish();
        dialog.addActionListener(this);
        System.runFinalization();
        System.gc();
        dialog.setTitle("Surface Curvature Analysis");
        return super.start();
    }
    
    public void actionPerformed(ActionEvent e)
    {
         if(e.getID() == WindowEvent.WINDOW_CLOSING)
             cancel();
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
         gotoNextStep();
    }

    public void finish()
    {
         super.finish();
         processSurface.panel.saveButton.setEnabled(false);
    }

    public void cancel()
    {
        super.cancel();
        processSurface.panel.saveButton.setEnabled(false);
    }

    public boolean exportView()
    {
          return true;
    }
    
    public boolean checkPrerequisites()
    {
        StsSurfaceClass surfaceClass = (StsSurfaceClass)model.getStsClass("com.Sts.DBTypes.StsSurface");
        StsSurface[] surfaces = surfaceClass.getSurfaces();
        if(surfaces.length == 0)
        {
             reasonForFailure = "No surfaces available. Must have at least one surface that has been loaded in the 3D view.";
             return false;
        }
        return true;
    }
    
    static public void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        StsSurfaceCurvatureWizard curvatureWizard = new StsSurfaceCurvatureWizard(actionManager);
        curvatureWizard.start();
    }
    
    public void setSelectedSurface(StsSurface surface)
    {
          if(surface == null) return;
          selectedSurface = surface;
    }

    public StsSurface getSelectedSurface()
    {
          return selectedSurface;
    }
}