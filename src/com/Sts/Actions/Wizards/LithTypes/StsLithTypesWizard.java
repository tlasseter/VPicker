
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.LithTypes;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.*;

import java.awt.*;

public class StsLithTypesWizard extends StsWizard
{
    public StsLibrarySelect librarySelect = null;
    public StsLibraryDefine libraryDefine = null;
    public StsTypeDefine typeDefine = null;
    public StsTypeCreate typeCreate = null;

    private StsTypeLibrary selectedLibrary = null;

    private StsWizardStep[] mySteps =
    {
        librarySelect = new StsLibrarySelect(this),
        libraryDefine = new StsLibraryDefine(this),
        typeDefine = new StsTypeDefine(this),
        typeCreate = new StsTypeCreate(this)
    };

    public StsLithTypesWizard(StsActionManager actionManager)
    {
        super(actionManager);
        addSteps(mySteps);
    }

    public boolean start()
    {
        System.runFinalization();
        System.gc();

        dialog.setTitle("Lithologic Library Setup");
        dialog.getContentPane().setSize(300, 300);

        if(!super.start()) return false;
        return true;
    }

    public void createNewType()
    {
        gotoStep(typeDefine);
        enableFinish();
    }

    public void createNewLibrary()
    {
        gotoStep(libraryDefine);
        enableFinish();
    }

    public boolean end()
    {
        return super.end();
    }

    public void previous()
    {
        if(currentStep == typeDefine)
        {
            gotoStep(librarySelect);
        }
        else
            gotoPreviousStep();
    }

    public void next()
    {
        if(currentStep == librarySelect)
        {
            if(selectedLibrary == null)
                new StsMessage(this.frame,StsMessage.WARNING, "No library selected, either select one or create a new one.");
            else
                gotoStep(typeDefine);
        }
        else if(currentStep == libraryDefine)
        {
            gotoStep(librarySelect);
            librarySelect.panel.reinitialize();
        }
        else
            gotoNextStep();
    }

    public void finish()
    {
        super.finish();
    }

    public StsTypeLibrary getSelectedLibrary()
    {
        if(selectedLibrary == null)
            selectedLibrary = StsTypeLibrary.getCreateDefaultLibrary();
        return selectedLibrary;
    }
    public void setSelectedLibrary(StsTypeLibrary lib)
    {
        selectedLibrary = lib;
    }

    public String getLibraryName()
    {
        return libraryDefine.getLibraryName();
    }

    public String getTypeName()
    {
        return typeDefine.getTypeName();
    }

    public Color getTypeColor()
    {
        return typeDefine.getTypeColor();
    }

    public byte[] getTypeTexture()
    {
        return typeDefine.getTypeTexture();
    }

}
