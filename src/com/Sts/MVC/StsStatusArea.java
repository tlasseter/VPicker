package com.Sts.MVC;
/**
 * <p>Title:        StsMessageFile</p>
 * <p>Description:  Class used to configure and maintain the status area. The status
 * area is used to present workflow step specific action controls and progress indicators.</p>
 * <p>Copyright:    Copyright (c) 2001</p>
 * <p>Company:      4D Systems LLC</p>
 * @author T.J.Lasseter
 * @version 1.0
 */

import com.Sts.Actions.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.border.*;
import java.awt.*;
import java.awt.event.*;

public class StsStatusArea extends StsJPanel implements StsStatusUI
{
    StsActionManager actionManager;
    /** The status area panel */
 //   static public StsStatusArea statusArea;
    /** User interface update delay = 5 */
    static public int UI_UPDATE_DELAY = 5;

	protected JProgressBar progress = new JProgressBar();
	protected JLabel status = new JLabel();
	protected JLabel titleBar = new JLabel();
    protected StsButton actionButton = new StsButton("End", this, "end");
    protected JButton abortButton = new JButton();
    protected JButton undoButton = new JButton();
    protected StsButton endAllButton = new StsButton("End all", this, "endAll");
    protected StsButton functionButton;

    private boolean statusOn = false;
    private boolean actionButtonOn = false;
    private boolean abortButtonOn = false;
	private boolean undoButtonOn = false;
    private StsAction endAllAction;
	private ActionListener actionListener = null;

    private boolean progressOn = false;
    private boolean titleOn = false;
	private int scaling = 1;

    /**
     * Default status area constructor
     * */
    public StsStatusArea(StsActionManager actionManager)
    {
        super(false);
        this.actionManager = actionManager;
        try { jbInit(); }
        catch(Exception e) { }
    }

    private void jbInit() throws Exception
    {
        Color bg = progress.getBackground();
      	Border border = BorderFactory.createLoweredBevelBorder();

    	titleBar.setBorder(border);
    	titleBar.setBackground(bg);

        status.setBackground(bg);
    	status.setBorder(border);
//        StsMessageFiles.logMessage(" ");
        functionButton = new StsButton("No Action", actionManager, "endCurrentAction");
        abortButton.setText(" ");
        undoButton.setText(" ");
        // actionButton.setText(" ");
        setBackground(bg);
        // Dimension preferredSize = actionButton.getPreferredSize();
        // Insets borderInset = border.getBorderInsets(status);
        // preferredSize.height += (borderInset.bottom + borderInset.top);
        // setPreferredSize(preferredSize);
        addAll();
		setVisible(true);
	}

    /**
     * Static add all status area components, title, progress bar, and buttons
     */
	public void staticAddAll() { addAll(); }

    /**
     * Add all status area components, title, progress bar, and buttons
     */
    public void addAll()
    {
        StsToolkit.runLaterOnEventThread
        (
            new Runnable()
            {
                public void run()
                {
                    gbc.anchor = GridBagConstraints.WEST;
                    gbc.fill = GridBagConstraints.NONE;
                    StsJPanel functionPanel = new StsJPanel();
                    functionPanel.addToRow(functionButton);
                    if( titleOn )functionPanel.addToRow(titleBar);
                    if(statusOn) functionPanel.addToRow(status);
                    if( progressOn ) functionPanel.addToRow(progress);
                    addToRow(functionPanel);

                    StsJPanel buttonPanel = new StsJPanel();
                    if( actionButtonOn) buttonPanel.addToRow(actionButton);
                    if( abortButtonOn) buttonPanel.addToRow(abortButton);
                    if( undoButtonOn) buttonPanel.addToRow(undoButton);
                    if(endAllAction != null) buttonPanel.addToRow(endAllButton);
                    addToRow(buttonPanel);

                    validate();
                }
            }
        );
    }

    public void removeAll()
    {
        StsToolkit.runLaterOnEventThread
        (
            new Runnable()
            {
                public void run()
                {
                    doRemoveAll();
                }
            }
        );
    }

    private void doRemoveAll()
    {
        super.removeAll();
    }
    /**
     * Add the progress bar
     */
    public void addProgress()
    {
        removeAll();
        progressOn = true;
        addAll();
    }

    /**
     * Remove the progress bar
     */
    public void removeProgress()
    {
    	removeAll();
        progressOn = false;
        setProgress(0);
        setMinimum(0);
        setMaximum(1);
        addAll();
    }

    /**
     * Set the title
     * @param title_ the title string
     */
    public void setTitle(String title_)
    {
        final String title = title_;
        StsToolkit.runLaterOnEventThread
        (
            new Runnable()
            {
                public void run()
                {
                    removeAll();
                    titleOn = true;
                    titleBar.setText(title);
                    addAll();
                    paintImmediately();
                }
            }
        );
    }

    private void paintImmediately()
	{
		paintImmediately(0, 0, getWidth(), getHeight());
	}

	public void setStatus(String text_)
	{
        final String text = text_;
        StsToolkit.runLaterOnEventThread
        (
            new Runnable()
            {
                public void run()
                {
                    status.setText(text);
                    statusOn = true;
                    paintImmediately();
                }
            }
        );
	}

	public void clearStatus()
	{
        StsToolkit.runLaterOnEventThread
        (
            new Runnable()
            {
                public void run()
                {
                    status.setText("");
                    statusOn = false;
                    paintImmediately();
               }
            }
        );
	}

    public void end()
    {
        if(actionManager != null) actionManager.endAction();
    }

    public void endAll()
    {
        if(actionManager != null) actionManager.endAction();
    }

    /**
     * Add the Add Button
     */

    public void setFunctionButtonText(String text_)
    {
        final String text = text_;
        StsToolkit.runLaterOnEventThread(new Runnable() { public void run() { getFunctionButton().setText(text); } } );
    }

    public JButton getFunctionButton() { return functionButton; }

    /**
     * Add the Abort Button
     */
	public void addAbortButton() { abortButtonOn = true; }
    public JButton getAbortButton() { return abortButton; }

    /**
     * Add the Undo Button
     */
	public void addUndoButton() { undoButtonOn = true; }
    /**
     * Add the EndAll Button
     */
    public void addEndAllButton(StsAction action)
    {
        endAllAction = action;
    }

    /**
     * Add the RemoveAll Button
     */
	public void removeAllButtons()
	{
		actionButton = null;
		abortButtonOn = false;
		undoButtonOn = false;
		endAllAction = null;
		addAll();
		paintImmediately();
	}
    /**
     * Add the End Button
     */
	public void addEndButton()
	{
		actionButton.setText("End");
		actionButton.setActionCommand("End");
	}

    /**
     * Add a action listener to the Add Button
     * @param newActionListener the action listener
     * @param text button text
     * @param visible button visibility
     * @return JButton
     */
	public JButton addActionButtonListener(ActionListener newActionListener, String text, boolean visible)
    {
        ActionListener[] listeners = actionButton.getActionListeners();
        if(listeners != null)
        {
            for(int n = 0; n < listeners.length; n++)
                actionButton.removeActionListener(listeners[n]);
        }
		actionButton.addActionListener(newActionListener);
		actionListener = newActionListener;
        actionButton.setText(text);
        actionButton.setActionCommand(text);
    	actionButton.setBackground(progress.getBackground());

		if(visible)
		{
			removeAll();
			actionButtonOn = true;
			addAll();
		}
        return actionButton;
    }

    /**
     * Add a action listener to the End Button
     * @param text button text
     * @param visible button visibility
     * @return JButton
     */
	public JButton addEndButtonListener(ActionListener newEndListener, String text, boolean visible)
    {
		actionButton.addActionListener(newEndListener);
        actionButton.setText(text);
        actionButton.setActionCommand(text);
    	actionButton.setBackground(progress.getBackground());

		if(visible)
		{
			removeAll();
			actionButtonOn = true;
			addAll();
		}

		return actionButton;
    }
    /**
     * Add a action listener to the Undo Button
     * @param text button text
     * @param visible button visibility
     * @return JButton
     */
	public JButton addUndoButtonListener(ActionListener newUndoListener, String text, boolean visible)
    {
		undoButton.addActionListener(newUndoListener);
        undoButton.setText(text);
        undoButton.setActionCommand(text);
    	undoButton.setBackground(progress.getBackground());

		if(visible)
		{
			removeAll();
			undoButtonOn = true;
			addAll();
		}

        return undoButton;
    }
    /**
     * Add a action listener to the Abort Button
     * @param text button text
     * @param visible button visibility
     * @return JButton
     */
	public JButton addAbortButtonListener(ActionListener newAbortListener, String text, boolean visible)
    {
		abortButton.addActionListener(newAbortListener);
        abortButton.setText("Abort");
        abortButton.setActionCommand("Abort");
    	abortButton.setBackground(progress.getBackground());

		if(visible)
		{
			removeAll();
			abortButtonOn = true;
			addAll();
		}

        return abortButton;
    }
    /**
     * Set action button command
     * @param command action command
     */
	public void setActionButtonCommand(String command)
	{
   	    removeAll();
		actionButton.setActionCommand(command);
		actionButton.setText(command);
		addAll();
	}

    /**
     * Set statusArea to text only
     */
    public void textOnly()
    {
        StsToolkit.runLaterOnEventThread
        (
            new Runnable()
            {
                public void run()
                {
                    removeAll();
                    actionButtonOn = false;
                    abortButtonOn = false;
                    undoButtonOn = false;
                    endAllAction = null;
                    titleOn = false;
                    progressOn = false;
                    addAll();
                }
            }
        );
    }
    /**
     * Set statusArea progress bar maximum
     * @params max maximum value
     */
    public void setMaximum(float max)
	{
		setMaximum((int)max);
	}

    /**
     * Set statusArea minimum progress bar value
     * @params min minimum value
     */
    public void setMinimum(float min) { setMinimum((int)min); }

    /**
     * Set progress bar value
     * @params f value of progress bar
     */
	public void setProgress(float f)
	{
		progress.setValue((int)(f*scaling));
		progress.repaint();
    }
/*
    public void setProgress(float f)
	{
        doUpdate(new updateProgress((int)(f*scaling)));
    }

    private void doUpdate(Runnable update)
    {
        if( SwingUtilities.isEventDispatchThread() ) update.run();
        else SwingUtilities.invokeLater(update);
    }

    class updateMax implements Runnable
    {
        int max;
        updateMax(int i) { max = i; }
        public void run() { progress.setMaximum(max); }
    }
*/
    /**
     * Set maximum progress bar value
     * @param max_ maximum
     */
    public void setMaximum(int max_)
    {
        scaling = 1;
        final int max = max_;
        StsToolkit.runLaterOnEventThread(new Runnable() { public void run() { progress.setMaximum(max); } } );
    }

    /**
     * Set the progress bar minimum
     * @param min_ minimum
     */
    public void setMinimum(int min_)
    {
        scaling = 1;
        final int min = min_;
        StsToolkit.runLaterOnEventThread(new Runnable() { public void run() { progress.setMinimum(min); } } );
    }

    /**
     * Set the progress bar value
     * @param n_ value
     */

	public void setProgress(int n_)
	 {
        final int n = n_;
        StsToolkit.runLaterOnEventThread
        (
            new Runnable()
            {
                public void run()
                {
                     progress.setValue(n);
                     paintImmediately();
                }
            }
        );
    }

    /**
     * Get the progress bar value
     * @return value
     */
    public int getProgress()
    {
    	return progress.getValue();
    }
/*
    class updateText implements Runnable
    {
        String string;
        updateText(String msg) { string = msg; }
        public void run() { StsMessageFiles.logMessage(string); }
    }
*/
    /**
     * Write the message to the screen
     * @param msg message
     */
	public void setText(String msg) { setText(msg, 0, true); }
    /**
     * Write the message to the screen
     * @param msg message
     * @param msec delay in milliseconds
     */
	public void setText(String msg, int msec) { setText(msg, msec, true); }
    /**
     * Write the message to the screen and log and sleep to let the user read the message
     * @param msg message
     * @param msec delay in microseconds
     * @param writeToLog true to write to log
     */
    public void setText(String msg, int msec, boolean writeToLog)
    {
		StsMessageFiles.logMessage(msg);
 //       doUpdate(new updateText(msg));
        if (msec>0) sleep(msec);
//        if (writeToLog && logFile != null) logFile.writeln(msg);
    }
    /**
     * Threaded sleep
     * @param msec length of sleep in microseconds
     */
	public void sleep(int msec)
    {
        if (msec<1) msec = 1;
        if (msec>10000) msec = 10000;
        try { Thread.currentThread().sleep(msec); }
        catch(Exception e) { }
    }


	static public void main(String[] args)
	 {
         StsModel model = StsModel.constructor();
         StsStatusArea statusArea = new StsStatusArea(model.mainWindowActionManager);
		 statusArea.addProgress();
		 com.Sts.Utilities.StsToolkit.createDialog(statusArea, false, 300, 30);
//		 statusArea.setTitle("title of panel .....");
		 statusArea.setMaximum(1000);
		 statusArea.setMinimum(0);
		 statusArea.setStatus("Status text...");

		 for( int i=0; i<1000; i++ )
		 {
			 statusArea.setProgress(i);
			 try { Thread.currentThread().sleep(1); }
			 catch(Exception e) { }
		}
    }
//    public void setLogFile(StsLogFile logFile) { this.logFile = logFile; }
}
