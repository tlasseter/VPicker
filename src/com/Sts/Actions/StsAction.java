
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions;

import com.Sts.MVC.*;
import com.Sts.MVC.View3d.StsView;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import java.awt.event.*;

/** Abstract class which must be subclassed.  Actions using the glPanel3d should
 *  subclass from StsAction3d @see com.Sts.Actions.StsAction3d.
 *  Responds to start, performMouseAction, and end calls from the action Manager.
 *  Start may in turn call run to run a separate thread for the action.
 */
abstract public class StsAction
{
    /** controller for this action */
	public StsActionManager actionManager;
    /** convenience copy of model */
	public StsModel model;
    /** panel for this action */
 //   public StsGLPanel glPanel;
    /** status area where temporary buttons and progress bar are placed */
	protected StsStatusArea statusArea;
    /** Failure reason */
	protected String reasonForFailure = null;
    /** indicates action has started */
    protected boolean started = false;
    /** Indicates this is a passive action, so it can interrupt another action without consequences.
     *  If false, then the user must approve the interrupt.
     */
    protected boolean canInterrupt = false;
    /** can be repeated */
	protected boolean repeatable = false;

	public StsAction()
	{
	}

    public StsAction(StsActionManager actionManager)
    {
	    this.actionManager = actionManager;
    	setModel(actionManager.getModel());
//		this.glPanel = actionManager.getGLPanel();
    }

    public StsAction(StsActionManager actionManager, boolean canInterrupt)
    {
		this(actionManager);
		this.canInterrupt = canInterrupt;
    }

    public StsAction(StsModel model)
    {
		setModel(model);
    }

    public StsAction(StsModel model, boolean canInterrupt)
    {
		setModel(model);
        this.canInterrupt = canInterrupt;
    }

	public void setModel(StsModel model)
	{
		if(model == null)
		{
			StsException.systemError("StsAction.setModel() failed. Model is null.");
			return;
		}
		this.model = model;
//        this.windowActionManager = model.mainWindowActionManager;
//        this.glPanel = model.glPanel3d;
        setStatusArea();
		if(statusArea != null) statusArea.textOnly();
	}

    private void setStatusArea()
    {
        if(model == null || model.win3d == null || model.win3d.statusArea == null) return;
        this.statusArea = model.win3d.statusArea;
    }

    public void setActionButtonText()
    {
        String val = getName();
        if (statusArea != null && val != null)
            statusArea.setFunctionButtonText(val);
        else
            System.out.println("setActionButtonText null found "+ (val ==null ? "val null":val) + (statusArea == null ? "status null":statusArea));
    }

    public void setActionButtonTextToEnd()
    {
        if(statusArea != null)
            statusArea.setFunctionButtonText("End");
    }
    
    public StsModel getModel() { return model; }
    public StsActionManager getActionManager() { return actionManager; }

    public boolean checkPrerequisites() { return true; }
    public String getReasonForFailure() { return reasonForFailure; }
    public boolean isStarted() { return started; }
    public void setStarted(boolean b) { this.started = b; }
    public boolean canInterrupt() { return canInterrupt; }
    public void setCanInterrupt(boolean b) { canInterrupt = b; }
	public boolean isRepeatable() { return repeatable; }
    public void setIsRepeatable(boolean b) { repeatable = b; }
    public boolean start() { return true;}
    public boolean performMouseAction(StsMouse mouse, StsGLPanel glPanel) { return true; }
    public boolean performMouseAction(StsMouse mouse, StsView view) { return true; }
    public boolean end() { return true;}
	public void addAbortButton() { statusArea.addAbortButton(); }
	public void addUndoButton() { statusArea.addUndoButton(); }
	public void addEndButton() { statusArea.addEndButton(); }
    public boolean checkStartAction() { return true; }
 //   public boolean hasOpenTransaction() { return openTransaction; }
 //   public void setOpenTransaction(boolean open) { openTransaction = open; }
    public boolean keyPressed(KeyEvent e, StsMouse mouse) { return true; }
    public boolean keyReleased(KeyEvent e, StsMouse mouse, StsGLPanel glPanel) { return true; }
	public boolean keyReleased(KeyEvent e, StsMouse mouse, StsGLJPanel glJPanel) { return true; }

    public void viewChanged(Object viewObject) { }

    /** implement as needed in subclasses to handle adding/removing toolbar
     *  these will be called by action manager when an action is started/stopped.
     */
	public void checkAddToolbar() { }

    /** implement as needed in subclasses to handle adding/removing toolbar
     *  these will be called by action manager when an action is started/stopped.
     */
    public void clearToolbar()
    {
    }

	public void addEndAllButton()
	{
		// repeatable = false;
		statusArea.addEndAllButton(this);
	}

    public void endAll()
    {
        System.out.println("end all called");
    }

    public void setActionButtonCommand(String command)
	{
		statusArea.setActionButtonCommand(command);
	}

	public void initializeRepeatAction(StsAction lastAction) { }
/*
	public void fireChangeEvent()
	{
		glPanel.getActionManager().fireChangeEvent();
	}
*/
	public void repeatAction()
	{
        model.repeatAction();
	}

	public void logMessage(String msg)
	{
		StsMessageFiles.logMessage(msg);
	}

	public void infoMessage(String msg)
	{
		StsMessageFiles.infoMessage(msg);
	}

	public void errorMessage(String msg)
	{
		StsMessageFiles.errorMessage(msg);
	}
    public String getName()
    {
        return StsToolkit.getSimpleClassname(this);
    }

    public boolean canBeInterrupted()
    {
        model.commit();
        return true;
    }

/*
	public void addAbortButton()
    {
        JButton abortButton = status.addAbortButton();
        abortButton.addActionListener(new AbortButtonActionListener());
	}

    protected class AbortButtonActionListener implements ActionListener
    {
        public void actionPerformed(ActionEvent e)
        {
			glPanel.getActionManager().endCurrentAction();
        }
    }

	public void addEndButton()
    {
        JButton endButton = status.addEndButton();
        endButton.addActionListener(new EndButtonActionListener());
	}

    protected class EndButtonActionListener implements ActionListener
    {
        public void actionPerformed(ActionEvent e)
        {
			if(e.getActionCommand().equals("End"))
			    glPanel.getActionManager().endCurrentAction();
        }
    }

	public void addEndAllButton()
    {
        JButton endAllButton = status.addEndAllButton();
        endAllButton.addActionListener(new EndAllButtonActionListener());
	}

    protected class EndAllButtonActionListener implements ActionListener
    {
        public void actionPerformed(ActionEvent e)
        {
			glPanel.getActionManager().endAllAction();
        }
    }
*/
/*
    public StsGLPanel getGlPanel()
    {
        return glPanel;
    }

    public StsGLPanel3d getGlPanel3d()
    {
         return (StsGLPanel3d)glPanel;
    }

    public void setGlPanel(StsGLPanel glPanel)
    {
        this.glPanel = glPanel;
    }
*/
}
