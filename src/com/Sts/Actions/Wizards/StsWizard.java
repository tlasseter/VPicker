
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards;

import com.Sts.Actions.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;
import com.Sts.Help.HelpManager;

import java.awt.*;
import java.util.*;

/** Abstract class which must be subclassed for wizards.
 *  Wizard has a fixed dialog box in which successive panels
 *  are placed for each wizardStep.  The wizard dialog has
 *  Next and Back buttons to handle stepping between wizard steps.
 */

abstract public class StsWizard extends StsAction
{
	/** execution pipe tp ensure serial execution of processes and panel updates **/
//	protected StsExecutionPipe executionPipe = new StsExecutionPipe();
    /** the dialog box for this wizard */
    public StsWizardDialog dialog = null;
    /** the list of wizardSteps */
    protected Vector steps = new Vector(5);
    /** the current wizardStep */
    protected StsWizardStep currentStep = null;
    /** indicates wizard completion status */
    protected boolean success = false;
    /** frame on which the wizard is centered */
    public Frame frame;
    /** indicates this is a test: model and/or win3d are null */
    protected boolean isTest = false;
    /** indicates the next step */
    protected int nextStepNumber = 0;
    protected boolean complete = false;

    /** sets the current step to the previous step */
    abstract public void previous();
    /** sets the current step to the next step */
    abstract public void next();

    public StsWizard(StsActionManager actionManager)
    {
        super(actionManager);

        if(model != null && model.win3d != null)
            frame = model.win3d;
        else
        {
            frame = null;
            isTest = true;
        }
        constructDialog();
        enableFinish();
    }

    public StsWizard(StsActionManager actionManager, int width, int height)
    {
        super(actionManager);

        if(model != null && model.win3d != null)
            frame = model.win3d;
        else
        {
            frame = null;
            isTest = true;
        }
        constructDialog(width, height);
        disableFinish();
    }

    private boolean constructDialog(int width, int height)
    {
        dialog = StsWizardDialog.constructor(frame, this, false);
        dialog.setPreferredSize(width, height);
        if(dialog == null) return false;
        dialog.setDefaultCloseOperation(dialog.DO_NOTHING_ON_CLOSE);
        return true;
    }

    private boolean constructDialog()
    {
        dialog = StsWizardDialog.constructor(frame, this, false);
        if(dialog == null) return false;
        dialog.setDefaultCloseOperation(dialog.DO_NOTHING_ON_CLOSE);
        return true;
    }

    public int getNextStepNumber()
    {
        return nextStepNumber++;
    }

    public void addStep(StsWizardStep step)
    {
    	steps.addElement(step);
    }

    public void addSteps(StsWizardStep[] array)
    {
    	if( array != null )
        	for( int i=0; i<array.length; i++ )
            	addStep(array[i]);
    }

    public void enableInstructions()
    {
        for(int i=0; i<steps.size(); i++)
            ((StsWizardStep)steps.elementAt(i)).showInstructions();
    }

    public void setHelpSet()
    {
        dialog.setHelpSet("Wizards/" + getName() + ".hs");
    }

    public void setHelpSet(String hsName)
    {
        dialog.setHelpSet(hsName);
    }

    public void run()
	{
        StsToolkit.runWaitOnEventThread(new Runnable() { public void run() { start(); }});
    }

    public boolean start()
    {

    	if(!initialize()) return false;
		dialog.pack();
        dialog.setLocationRelativeTo(frame);
    	dialog.setVisible(true);
        complete = false;
        return true;
    }

    public boolean end()
    {
        if(currentStep != null)
            success = currentStep.end();
    	dialog.setVisible(false);
        dialog.dispose();
        model.refreshObjectPanel();
        if(success) model.setActionStatus(getClass().getName(), StsModel.STARTED);
		return success;
    }

    public void help()
    {
        // No Help Available
    }

    public boolean gotoStep(StsWizardStep wizardStep)
    {
        if(wizardStep == null)
        {
            enableFinish();
            return false;
        }

        if(wizardStep == currentStep)
            return true;
        
        if(currentStep != null)
        {
            if(!currentStep.end())
                return false;
            currentStep.hide();
        }

        currentStep = wizardStep;
        enablePrevious();
        enableNext();
        currentStep.setVisible(true);
        
        if(currentStep instanceof Runnable)
        {
            Thread thread = new Thread((Runnable)currentStep);
            thread.start();
            return true;
        }
        boolean ok = currentStep.start();
        dialog.repaint();
        return ok;
    }

    public boolean gotoNextStep()
    {
        StsWizardStep nextStep = getNextStep();
        return gotoStep(nextStep);
    }

    public boolean gotoPreviousStep()
    {
        StsWizardStep prevStep = getPreviousStep();
        return gotoStep(prevStep);
    }

    public boolean hasPrevious()
    {
        int index = currentStep.getStepNumber();
        return index > 0;
    }

    public boolean hasNext()
    {
        int index = currentStep.getStepNumber();
        return index < steps.size()-1;
    }

    public StsWizardStep getPreviousStep()
    {
        int index = currentStep.getStepNumber()-1;
        if(index < 0) return null;
        return (StsWizardStep)steps.elementAt(index);
    }

    public StsWizardStep getNextStep()
    {
        int index = currentStep.getStepNumber()+1;
        if(index >= steps.size()) return null;
        return (StsWizardStep)steps.elementAt(index);
    }

    public boolean gotoFirstStep()
    {
        return gotoStep((StsWizardStep)steps.firstElement());
    }

    private boolean initialize()
    {
        if(steps == null || steps.size() == 0)
        {
            StsException.systemError("StsWizard.initialize() failed. No wizard steps.");
            return false;
        }
        return gotoFirstStep();
//        currentStep = (StsWizardStep) steps.elementAt(0);
//        gotoStep(currentStep);
//        return true;
    }

    public void checkBlockButtons()
    {
        disablePrevious();
        disableNext();
    }

    public void checkUnblockButtons(boolean success)
    {
        enablePrevious();
        if(success)
            enableNext();
        else
            disableNext();
    }

    public void enableNext(boolean enable)
    {
        if(enable)
            enableNext();
        else
            disableNext();
    }

    public void enableNext()
    {
        dialog.enableNext(hasNext());
    }

    public void disableNext()
    {
        dialog.disableNext();
    }

    public void enablePrevious()
    {
        dialog.enablePrevious(hasPrevious());
    }

    public void disablePrevious()
    {
        dialog.disablePrevious();
    }

    public void enableCancel()
    {
        dialog.enableCancel();
    }

    public void disableCancel()
    {
        dialog.disableCancel();
    }

    public void enableFinish()
    {
        dialog.enableFinish();
    }

    public void disableFinish()
    {
        dialog.disableFinish();
    }

    public void finish()
    {
//    	success = currentStep.end();
//		executionPipe.stop();
        actionManager.endCurrentAction();
	    currentStep = null;
    }

    public boolean isComplete() { return complete; }

    /** call this method when a wizard has finished loading objects */
    public void completeLoading(boolean success)
    {
        this.success = success;
        disableCancel();
        model.enableDisplay();
        if(success)
        {
            model.setActionStatus(getClass().getName(), StsModel.STARTED);
            enableFinish();
        }
        else
            enableCancel();
    }

    public void cancel()
	{
		if(currentStep != null)
			currentStep.cancel();
		actionManager.cancelCurrentAction();
		dialog.setVisible(false);
		dialog.dispose();
		model.refreshObjectPanel();
        model.enableDisplay();
		try
		{
//            System.out.println("StsWizard.cancel().threadComplete() disabled.");
//            actionManager.threadComplete();
		}
		catch(Exception e)
        {
            StsException.outputWarningException(this, "cancel", e);
        }
	}
/*
	public void stopProcess()
	{
//		executionPipe.killAllProcesses(this);
		actionManager.cancelCurrentAction();
	}
*/
    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel)
    {
//    	System.out.println("Mouse at (" + mouse.getMousePoint().x +"," + mouse.getMousePoint().y + ")");
        if( currentStep != null) currentStep.performMouseAction(mouse, glPanel);
    	return true;
    }

    public void rebuild()
    {
        dialog.pack();
    }

    // If this is a GUI test, then frame is null; otherwise it's model.win3d
    public boolean isTest() { return isTest; }

    public void itemSelected() {}
}
