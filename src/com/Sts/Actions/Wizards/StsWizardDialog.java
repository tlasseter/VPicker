//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards;

import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;
import com.Sts.Help.HelpManager;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

/** dialog window for wizard in which all wizard step panels are shown.
 *  Has buttons for controlling operations: next, back, finish, cancel.
 *  Has headerPanel with vendor logo, info
 */
public class StsWizardDialog extends JDialog implements ActionListener
{
    StsWizard wizard;
    StsWizardStep wizardStep = null;
    /** contentPanel holds mainPanel, progressBar (optional), and buttonPanel */
    StsJPanel contentPanel = StsJPanel.addInsets();
    /** mainPanel holds headerPanel, infoPanel, and dialogPanel */
    StsJPanel mainPanel = StsJPanel.addInsets();
    /** buttonPanel has prev, next, finish, and cancel buttons.
     *  Each of these may be enabled/disabled.
     */
    StsJPanel buttonPanel = StsJPanel.addInsets();

    StsButton prevButton;
	StsButton nextButton;
    StsButton finishButton;
	StsButton cancelButton;
    StsButton helpButton;

    Border border = BorderFactory.createLoweredBevelBorder();

    static public final String PREVIOUS = "<< Back";
	static public final String NEXT = "Next >>";
	static public final String FINISH = "Finish";
	static public final String CANCEL = "Cancel";
    static public final String HELP = "Help";

    public StsButton[] buttons;

    Vector actionListeners = null;

    public StsWizardDialog()
	{
		this(null, null, false);
	}

	private StsWizardDialog(Frame frame, StsWizard wizard, boolean modal)
	{
		super(frame, modal);
        try
		{
            this.wizard = wizard;
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}

    static public StsWizardDialog constructor(Frame frame, StsWizard wizard, boolean modal)
    {
        StsWizardDialog wizardDialog;
        try
        {
            wizardDialog = new StsWizardDialog(frame, wizard, modal);
            wizardDialog.constructDialog();
            return wizardDialog;
        }
        catch(Exception e)
        {
            return null;
        }
    }

    static public StsWizardDialog constructor(Frame frame, boolean modal)
    {
        return constructor(frame, null, modal);
    }

    public void setPreferredSize(int width, int height)
    {
        setPreferredSize(new Dimension(width, height));
    }

    public void constructDialog()
	{
        createButtons();
        getContentPane().add(contentPanel);
        contentPanel.gbc.fill = GridBagConstraints.BOTH;
        contentPanel.gbc.weighty = 1.0;
        contentPanel.add(mainPanel);
        contentPanel.gbc.fill = GridBagConstraints.NONE;
        contentPanel.gbc.weighty = 0.0;
        contentPanel.add(buttonPanel);
    }

    public void setHelpSet(String hsName)
    {
        try
        {
            if(!HelpManager.getHelpManager(this).addHelpSet(hsName))
            {
                helpButton.setToolTipText("No wizard specific help available");
                helpButton.setEnabled(false);
                return;
            }
            helpButton.setEnabled(true);            
            HelpManager.setSourceHelpActionListener(helpButton, HelpManager.WIZARD, this);  // JavaHelp
        }
        catch(Exception ex)
        {
            System.err.println("Unable to set helpset: " + hsName);
            ex.printStackTrace();
        }
    }

    private void createButtons()
    {
        prevButton = new StsButton(PREVIOUS, wizard, "previous");
        nextButton = new StsButton(NEXT, wizard, "next");
        finishButton = new StsButton(FINISH, wizard, "finish");
        cancelButton = new StsButton(CANCEL, wizard, "cancel");
        helpButton = new StsButton(HELP, wizard, "help");
        buttons = new StsButton[] { prevButton, nextButton, finishButton, cancelButton, helpButton };
        for(int n = 0; n < buttons.length; n++)
            buttons[n].setMargin(new Insets(2, 5, 2, 5));
        helpButton.setEnabled(false);

		buttonPanel.setBorder(BorderFactory.createEtchedBorder());
        addButtons();
    }

    private void addButtons()
    {
        buttonPanel.addToRow(prevButton);
        buttonPanel.addToRow(nextButton);
        buttonPanel.addToRow(finishButton);
        buttonPanel.addToRow(cancelButton);
        buttonPanel.addEndRow(helpButton);
    }

    public void clearMainPanel()
    {
        mainPanel.removeAll();
    }

    public void setWizardStep(StsWizardStep wizardStep)
    {
        this.wizardStep = wizardStep;
        clearMainPanel();
        mainPanel.gbc.fill = GridBagConstraints.BOTH;
        mainPanel.gbc.weighty = 0;
        if(wizardStep.headerPanel != null) mainPanel.add(wizardStep.headerPanel);
        if(wizardStep.infoPanel != null) mainPanel.add(wizardStep.infoPanel);
        mainPanel.gbc.weighty = 1;
        if(wizardStep.dialogPanel != null) mainPanel.add(wizardStep.dialogPanel);
        wizard.rebuild();
    }

    public void addToolbar(StsToolbar tb)
	{
        buttonPanel.removeAll();
        buttonPanel.addToRow(tb);
        addButtons();
        pack();
    }

    public void removeToolbar(StsToolbar tb)
	{
		buttonPanel.remove(tb);
	}
    
    private class runEnableButton implements Runnable
    {
		String _label;
		boolean _b;

        runEnableButton(String label, boolean b)
        {
			_label = label;
			_b = b;
		}
		public void run()
        {
           for (int i = 0; i < buttons.length; i++)
           {
               if (buttons[i].getText().equalsIgnoreCase(_label))
               {
                   buttons[i].setEnabled(_b);
                   break;
               }
           }
		}
	}

    public void enableNext(boolean hasNext)
    {
        enableButton(nextButton, hasNext);
    }

    public void disableNext()
    {
        enableButton(nextButton, false);
    }

    public void enablePrevious(boolean hasPrevious)
    {
        enableButton(prevButton, hasPrevious);
    }

    public void disablePrevious()
    {
        enableButton(prevButton, false);
    }

    public void enableCancel()
    {
        enableButton(cancelButton, true);
    }

    public void disableCancel()
    {
        enableButton(cancelButton, false);
    }

    public void enableFinish()
    {
        enableButton(finishButton, true);
    }

    public void disableFinish()
    {
        enableButton(finishButton, false);
    }

    static public void enableButton(StsButton button_, boolean b_)
    {
        final StsButton button = button_;
        final boolean enable = b_;
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                button.setEnabled(enable);
            }
        };
        StsToolkit.runLaterOnEventThread(runnable);
    }

    public void enableButton(String label, boolean b)
	{

		runEnableButton runnable = new runEnableButton(label,b);
        StsToolkit.runWaitOnEventThread(runnable);
	}

	public JComponent getHeaderPane()
	{
		return wizardStep.headerPanel;
    }

	public JComponent getDialogPane()
	{
		return wizardStep.dialogPanel;
	}

	public JComponent getInfoPane()
	{
		return wizardStep.infoPanel;
	}

	public synchronized void removeActionListener(ActionListener l)
	{
		if (actionListeners != null && actionListeners.contains(l))
		{
			Vector v = (Vector) actionListeners.clone();
			v.removeElement(l);
			actionListeners = v;
		}
	}

	public synchronized void addActionListener(ActionListener l)
	{
		Vector v = actionListeners == null ? new Vector(2) :
			(Vector) actionListeners.clone();
		if (!v.contains(l))
		{
			v.addElement(l);
			actionListeners = v;
		}
	}

	protected void fireActionPerformed(ActionEvent e)
	{
		if (actionListeners != null)
		{
			Vector listeners = actionListeners;
			int count = listeners.size();
			for (int i = 0; i < count; i++)
			{
				( (ActionListener) listeners.elementAt(i)).actionPerformed(e);
			}
		}
	}

	public void actionPerformed(ActionEvent e)
	{
		fireActionPerformed(e);
	}

    protected void processWindowEvent(WindowEvent e)
	{
		if (e.getID() == WindowEvent.WINDOW_CLOSING)
		{
			fireActionPerformed(new ActionEvent(e.getSource(), e.getID(), StsWizardDialog.CANCEL));
            setVisible(false);
		    dispose();
        }
		else
		{
			super.processWindowEvent(e);
		}
	}
	/*
	   static public void main(String[] args)
	   {
		StsWizardDialog dialog = new StsWizardDialog();
		int width = 250;
		int height = 500;
//        Dimension size = new Dimension(250, 500);
//        dialog.setSize(width, height);
//        JDialog dialog = new JDialog();
//        JPanel panel = (JPanel)dialog.getDialogPane();
		JPanel panel = new JPanel(new GridBagLayout());
		dialog.setDialogPane(panel);
		StsGLPanel2d leftPanel = new StsGLPanel2d(width, height);
//        leftPanel.setMaximumSize(size);
//        leftPanel.setPreferredSize(size);
//        leftPanel.setBackground(Color.RED);
		panel.add(leftPanel, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));
		System.out.println("First panel added");
		dialog.pack();
		Dimension size = dialog.getSize();
	 System.out.println("Wizard dialog packed: " + size.width + " " + size.height);
		dialog.setVisible(true);
		try { Thread.currentThread().sleep(1000); } catch(Exception e) { }

		StsGLPanel2d ritePanel = new StsGLPanel2d(width, height);
//        JPanel ritePanel = new JPanel();
//        ritePanel.setMaximumSize(size);
//        ritePanel.setPreferredSize(size);
	  //       ritePanel.setBackground(Color.GREEN);
		panel.add(ritePanel, new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));
//        Container dialogContainer = dialog.getContentPane();
//        dialogContainer.add(panel);
		System.out.println("Second panel added.");
		dialog.pack();
		size = dialog.getSize();
	 System.out.println("Wizard dialog packed: " + size.width + " " + size.height);
//        dialog.setVisible(true);
	   }
	 }
	 */
	/*
	  static public void main(String[] args)
	  {
	 StsWizardDialog dialog = new StsWizardDialog();
	 int width = 250;
	 int height = 500;
	  JPanel panel = new JPanel(new GridBagLayout());
	 JPanel dummyPanel = new JPanel();
	 dummyPanel.setPreferredSize(new Dimension(width, height));
	 panel.add(dummyPanel, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));
	 dialog.setDialogPane(panel);

	 dialog.pack();
	 Dimension size = dialog.getSize();
//        System.out.println("Wizard dialog packed: " + size.width + " " + size.height);
	 dialog.setVisible(true);

	 try
	 {
	  Thread.currentThread().sleep(1000);
	 }
	 catch (Exception e)
	 {}

	 panel.remove(dummyPanel);
	 StsGLPanel2d leftPanel = new StsGLPanel2d(width, height);
	 panel.add(leftPanel,
		 new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH,
			  new Insets(0, 0, 0, 0), 0, 0));
	 System.out.println("First panel added");
	 dialog.pack();
	 size = dialog.getSize();
//        System.out.println("Wizard dialog packed: " + size.width + " " + size.height);

	 try
	 {
	  Thread.currentThread().sleep(1000);
	 }
	 catch (Exception e)
	 {}

	 StsGLPanel2d ritePanel = new StsGLPanel2d(width, height);
	 panel.add(ritePanel,
		 new GridBagConstraints(1, 0, 1, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.BOTH,
			  new Insets(0, 0, 0, 0), 0, 0));
	 System.out.println("Second panel added.");
	 dialog.pack();
	 size = dialog.getSize();
//        System.out.println("Wizard dialog packed: " + size.width + " " + size.height);
	  }
	 */
}
