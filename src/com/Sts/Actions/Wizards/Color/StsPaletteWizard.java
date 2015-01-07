
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Color;

//import java.awt.Cursor;
//import java.io.*;
//import java.util.jar.*;
//import java.net.*;
//import javax.jnlp.*;
//import java.awt.*;
//import java.awt.event.*;
//import javax.swing.*;
//import java.util.*;
//import java.io.File;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;

public class StsPaletteWizard extends StsWizard
{
    public StsPaletteSelect selectPalette = null;
    public StsPaletteLoad loadPalette = null;

    private StsWizardStep[] mySteps =
    {
        selectPalette = new StsPaletteSelect(this),
        loadPalette = new StsPaletteLoad(this)
    };

    public StsPaletteWizard(StsActionManager actionManager)
    {
        super(actionManager,500,500);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();
        dialog.setTitle("Color Setup");
        disableFinish();
        return super.start();
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
        if(currentStep == selectPalette)
            loadPalette.constructPanel();
        gotoNextStep();
    }

    // Persist any palettes that are not already in database
    public boolean createSpectrums(StsProgressTextPanel panel)
    {
        StsAbstractFile[] files = selectPalette.getSelectedFiles();
        String[] names = selectPalette.getSelectedPaletteNames();
        String dirname = selectPalette.getSelectedDirectory();

        try
        {
            disablePrevious();
            StsSpectrumClass spectrumClass = model.getSpectrumClass();
            for (int i = 0; i < names.length; i++)
            {
                if (spectrumClass.getSpectrum(names[i]) == null)
                {
                    // Add the palette
                   StsColor[] colors = selectPalette.panel.readPalette(files[i]);
                   StsSpectrum spectrum = null;
                   if(colors.length < 255)
                	   spectrum = StsSpectrum.constructor(names[i], colors, 255);
                   else
                	   spectrum = StsSpectrum.constructor(names[i], colors);
				   if(spectrum != null)
                       panel.appendLine("Loaded palette: " + names[i]);
                   else
                       panel.appendLine("Failed to load palette: " + names[i]);
                }
                else
                {
                    panel.appendLine("Palette: " + names[i] + " already in database.");
                }

                disableCancel();
            }
            enableFinish();
        }
        catch (Exception e)
        {
            StsException.outputException("StsPaletteWizard.createSpectrums() failed.", e, StsException.WARNING);
        }
        return true;
    }

    public void finish()
    {
        super.finish();
    }
}
