
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards;

import com.Sts.Actions.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;

import javax.swing.*;

/** defines the panels for this wizard step in a wizard sequence;
 *  subclass for concrete step and implement start and end methods
 *  for the step.
 */
abstract public class StsWizardStep extends StsAction
{
    /** wizard containing this wizardStep */
    protected StsWizard wizard;
    /** wizard dialog holding these wizardStep panels */
    public StsWizardDialog wizardDialog;
    /** panel on which user makes selections */
    protected JComponent dialogPanel;
    /** panel containing info on object(s) selected/defined */
    protected JComponent infoPanel;
    /** contains wizard description and vendor logo */
    protected JComponent headerPanel;
    /** index number for this wizardStep in sequence */
    protected int stepNumber;
    /** if not null, this is the current runnable for this step which can be canceled */
    protected StsProgressRunnable progressRunnable = null;
    /** indicates wizard step was successful */
    public boolean success = false;

    abstract public boolean start();
    abstract public boolean end();

    public StsWizardStep()
    {
    }

    public StsWizardStep(StsWizard wizard)
    {
    	this(wizard, null, null, null);
    }

    public StsWizardStep(StsWizard wizard, JComponent dialog)
    {
    	this(wizard, dialog, null, null);
    }

    public StsWizardStep(StsWizard wizard, JComponent dialog, JComponent info)
    {
        this(wizard, dialog, info, null);
    }

    public StsWizardStep(StsWizard wizard, JComponent panel, JComponent info, JComponent header)
    {
        initialize(wizard, panel, info, header);
    }

    public void initialize(StsWizard wizard, JComponent panel, JComponent info, JComponent header)
    {
        setModel(wizard.getModel());
        this.wizard = wizard;
        this.wizardDialog = wizard.dialog;
        this.dialogPanel = panel;
        this.infoPanel = info;
        this.headerPanel = header;
        stepNumber = wizard.getNextStepNumber();
    }

    public void run() { start(); }

    public void enableFinish()
    {
        wizard.enableFinish();
    }

    public void disableFinish()
    {
        wizard.disableFinish();
    }

    public void enablePrevious()
    {
        wizard.enablePrevious();
    }

    public void disablePrevious()
    {
        wizard.disablePrevious();
    }

    public void enableNext()
    {
        wizard.enableNext();
    }

    public void disableNext()
    {
        wizard.disableNext();
    }

    public void enableCancel()
    {
        wizard.enableCancel();
    }

    public void disableCancel()
    {
        wizard.disableCancel();
    }

    public int getStepNumber() { return stepNumber; }

    /** for now we need to override this method in concrete wizardSteps, but we should make it an abstract method. */
    public void cancel()
    {
        if(progressRunnable != null)
        {
            progressRunnable.cancel();
            progressRunnable = null;
        }
        else if(this instanceof StsProgressRunnable)
            doCancel();
    }

    public void doCancel()
    {

    }

    public void cancelRunnable()
    {
        if(progressRunnable == null) return;
        progressRunnable.cancel();
        progressRunnable = null;
    }

    public void setPanels(JComponent dialogPanel, JComponent headerPanel)
    {
        this.dialogPanel = dialogPanel;
        this.headerPanel = headerPanel;
    }

    public void setPanels(JComponent dialogPanel, JComponent headerPanel, JComponent infoPanel)
    {
        this.dialogPanel = dialogPanel;
        this.headerPanel = headerPanel;
        this.infoPanel = infoPanel;
    }
    /**
     * Show and Remove the instructions, if available
     */
    public void showInstructions() { if(headerPanel != null) ((StsHeaderPanel)headerPanel).showInfoText(); }
    public void removeInstructions() { if(headerPanel != null) ((StsHeaderPanel)headerPanel).hideInfoText(); }

    public JComponent getContainer()  { return dialogPanel; }
    public JComponent getInfoContainer()  { return infoPanel; }
    public JComponent getHdrContainer()  { return headerPanel; }
    public StsWizard getWizard() { return wizard; }
    public void updatePanel() { };

//    public void setDialogPanel(JComponent c)  { this.dialogPanel = c; }
    public void setInfoPanel(JComponent c)  { this.infoPanel = c; }
//    public void setHeaderPanel(JComponent c) {this.headerPanel = c;}

    public void setVisible(boolean isVisible)
    {
        wizardDialog.setWizardStep(this);
    }

    public void hide()
    {
        wizardDialog.clearMainPanel();
    }

    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
//    	System.out.println("Mouse at (" + mouse.getMousePoint().x +"," + mouse.getMousePoint().y + ")");
    	return true;
    }

    // if we are running wizard GUI standalone for testing, then the frame is null
    protected boolean isTest()
    {
        return wizard.isTest();
    }
}
