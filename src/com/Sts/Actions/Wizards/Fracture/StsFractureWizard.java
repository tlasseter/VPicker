//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Fracture;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.*;
import com.Sts.Utilities.StsMath;

import java.awt.event.*;

public class StsFractureWizard extends StsWizard
{
    public StsFractureSet selectedFractureSet = null;
    //   public StsFractureSet templateFractureSet = null;
    public StsColor stsColor = new StsColor(StsColor.RED);
    private StsFractureSetClass fractureSetClass;
    public StsSelectFractureSet selectFractureSet = null;
    public StsDefineFractureSet defineFractureSet = null;
    public StsBuildFractureSet buildFractureSet = null;

    private byte currentDisplayType;

    static String NONE = "none";

    public StsFractureWizard(StsActionManager actionManager)
    {
        super(actionManager, 400, 400);
        fractureSetClass = (StsFractureSetClass) model.getCreateStsClass(StsFractureSet.class);
        selectedFractureSet = (StsFractureSet) fractureSetClass.getCurrentObject();
        currentDisplayType = fractureSetClass.getDisplayType();
        if (currentDisplayType == StsFractureSetClass.DISPLAY_VOLUME)
            fractureSetClass.setDisplayTypeRadius();
        addSteps
                (
                        new StsWizardStep[]
                                {
                                        selectFractureSet = new StsSelectFractureSet(this),
                                        defineFractureSet = new StsDefineFractureSet(this),
                                        buildFractureSet = new StsBuildFractureSet(this)
                                }
                );
    }

    public boolean keyReleased(KeyEvent e, StsMouse mouse, StsGLPanel glPanel)
    {
        if (currentStep == buildFractureSet)
            return buildFractureSet.keyReleased(e, mouse, glPanel);
        return true;
    }

    public boolean start()
    {
        disableFinish();
        System.runFinalization();
        System.gc();
        dialog.setTitle("Fracture Definition");
        return super.start();
    }

    public void setSelectedFractureSet(Object fractureSet)
    {
        if (fractureSet == NONE) return;
        selectedFractureSet = (StsFractureSet) fractureSet;
    }
/*
    public Object getTemplateFractureSet()
        {
            return templateFractureSet;
        }
        public void setTemplateFractureSet(Object fractureSet)
        {
            if(fractureSet == NONE)
                templateFractureSet = null;
            else
            {
                templateFractureSet = (StsFractureSet)fractureSet;
                if(templateFractureSet != null)
                     selectedFractureSet.setCursorZSurfaces(templateFractureSet.getCursorZSurfaces());
            }
        }
    */

    public Object getSelectedFractureSet()
    {
        return selectedFractureSet;
    }

    public void createNewFractureSet()
    {
        selectedFractureSet = StsFractureSet.constructor();
        gotoStep(defineFractureSet);
    }

    /**
     * add the current fractureSet to the model and create a new scratch one for the next fractureSet
     */
    public void addToModel()
    {
        model.add(selectedFractureSet);
        selectedFractureSet = StsFractureSet.constructor();
    }

    public boolean end()
    {
        if (currentDisplayType == StsFractureSetClass.DISPLAY_VOLUME)
            fractureSetClass.setDisplayTypeVolume();
        return super.end();
    }

    public void previous()
    {
        gotoPreviousStep();
    }

    public void next()
    {
        if (currentStep == selectFractureSet)
        {
            if (selectedFractureSet == null)
            {
                disableFinish();
                disableNext();
                new StsMessage(frame, StsMessage.INFO, "No fracture set selected or defined.");
                return;
            } else
            {
                enableNext();
                model.setCurrentObject(selectedFractureSet);
                model.win3dDisplay();
                gotoStep(buildFractureSet);
            }
        } else if (currentStep == defineFractureSet)
        {
            if (selectedFractureSet == null)
            {
                disableFinish();
                disableNext();
                new StsMessage(frame, StsMessage.INFO, "No horizon defined.");
                return;
            }
            if (!selectedFractureSet.isPersistent())
                model.add(selectedFractureSet);
            model.setCurrentObject(selectedFractureSet);
            model.win3dDisplay();
            enableNext();
            gotoStep(buildFractureSet);
        } else if (currentStep == buildFractureSet)
        {
            model.setCurrentObject(selectedFractureSet);
            disableNext();
            enableFinish();
            gotoNextStep();
        }
    }

    public void setFractureSetColor(Object color)
    {
        stsColor = (StsColor) color;
    }

    public Object getFractureSetColor()
    {
        return stsColor;
    }

    public void setFractureSetName(String name)
    {
        selectedFractureSet.setName(name);
    }

    public String getFractureSetName()
    {
        if (selectedFractureSet == null) return NONE;
        return selectedFractureSet.getName();
    }

    public void finish()
    {
        super.finish();
    }

    public Object[] getFractureSetList()
    {
        if (fractureSetClass.getSize() == 0)
            return new Object[]{NONE};
        else
        {
            byte projectZDomain = model.getProject().zDomain;
            Object[] fractureSetObjects = fractureSetClass.getElements();
            int nFractureSetObjects = fractureSetObjects.length;
            StsFractureSet[] domainFractureSets = new StsFractureSet[nFractureSetObjects];
            int nDomainFractureSets = 0;
            for(int n = 0; n < nFractureSetObjects; n++)
            {
                StsFractureSet fractureSet = (StsFractureSet)fractureSetObjects[n];
                byte zDomain = fractureSet.getZDomain();
                if(zDomain != projectZDomain) continue;
                domainFractureSets[nDomainFractureSets++] = fractureSet;
            }
            domainFractureSets = (StsFractureSet[])StsMath.trimArray(domainFractureSets, nDomainFractureSets);
            return domainFractureSets;
        }
    }

    public Object[] getFractureSetListWithNone()
    {
        if (fractureSetClass.getSize() == 0)
            return new Object[]{NONE};
        else
        {
            Object[] newList = new Object[fractureSetClass.getSize() + 1];
            newList[0] = NONE;
            for (int i = 0; i < fractureSetClass.getSize(); i++)
                newList[i + 1] = fractureSetClass.getElement(i);
            return newList;
        }
    }

    static public void main(String[] args)
    {
        StsModel model = StsModel.constructor();
        StsActionManager actionManager = new StsActionManager(model);
        StsFractureWizard fractureWizard = new StsFractureWizard(actionManager);
        fractureWizard.start();
    }
}
