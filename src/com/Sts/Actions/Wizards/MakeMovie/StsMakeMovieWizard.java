
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.MakeMovie;

import com.Sts.Actions.Wizards.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;

import javax.imageio.*;
import java.awt.image.*;
import java.io.*;

public class StsMakeMovieWizard extends StsWizard
{
    public StsMakeMovieDefine movieDefine = null;
    public StsMakeMovieCreate movieCreate = null;

    private StsFile[] jpgFiles = new StsFile[0];

    private StsWizardStep[] mySteps =
    {
        movieDefine = new StsMakeMovieDefine(this),
        movieCreate = new StsMakeMovieCreate(this)
    };

    public StsMakeMovieWizard(StsActionManager actionManager)
    {
        super(actionManager, 400, 800);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        System.setSecurityManager(null);
        dialog.setTitle("Capture and Create a Quicktime Movie");
        
        // Verify that capture mode is off.
        if(this.getModel().win3d.captureMovie())
        {
        	this.getModel().win3d.outputMovie();
        	((StsMediaToolbar)this.getModel().win3d.getToolbarNamed(StsMediaToolbar.NAME)).resetMovieBtn();
        }
    	
        if(!super.start())
            return false;

        return true;
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
        if(currentStep == movieDefine)
        {
        	// If capturing, turn it off.
        	if(movieDefine.panel.isCapturing())
        		movieDefine.panel.startStop();
            movieCreate.constructPanel();
        }
        gotoNextStep();
    }

    public void finish()
    {
        super.finish();
    }

    public void addFile(StsFile file, String outputDirectory)
    {
        jpgFiles = (StsFile[]) StsMath.arrayAddElement(jpgFiles, file);
    }

    public void removeFile(StsFile file)
    {
        jpgFiles = (StsFile[]) StsMath.arrayDeleteElement(jpgFiles, file);
    }

    public void removeFiles()
    {
        if(jpgFiles == null) return;
        for(int i = 0; i<jpgFiles.length; i++)
            jpgFiles = (StsFile[]) StsMath.arrayDeleteElement(jpgFiles, jpgFiles[i]);

        jpgFiles = null;
    }

    public String[] getFilenames()
    {
        String[] filenames = new String[jpgFiles.length];
        for(int i=0; i<jpgFiles.length; i++)
        {
            StsFile file = (StsFile)(jpgFiles[i]);
            filenames[i] = file.getPathname();
        }
        return filenames;
    }

    public String getOutputMovieName()
    {
        return "file:" + model.getProject().getRootDirString() + model.getProject().getMediaDirString() + movieDefine.panel.getMovieName();
    }
    /* Assumes all selected images are the same size. If not any off size frames will be wrapped. */
    public boolean createMovie(StsProgressPanel panel)
    {
        try
        {
            File tFile = new File(jpgFiles[0].getPathname());
            BufferedImage buf = ImageIO.read(tFile);
            int width = buf.getWidth();
            int height = buf.getHeight();
            if(StsImagesToMovie.imagesToMovie(width, height, movieDefine.panel.getFrameRate(), getOutputMovieName(), getFilenames(), panel))
            {
                if(StsYesNoDialog.questionValue(model.win3d,"Do you want to delete all the frame images used in this movie from disk?"))
                {
                    for (int i = 0; i < jpgFiles.length; i++)
                       jpgFiles[i].delete();
                }
                return true;
            }
            else
            {
                new StsMessage(model.win3d,StsMessage.ERROR,"Unable to produce movie...Error occurred during processing.");
                return false;
            }
        }
        catch(Exception e)
        {
            return false;
        }
    }
}
