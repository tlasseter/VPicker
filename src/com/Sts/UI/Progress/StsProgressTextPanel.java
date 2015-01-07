package com.Sts.UI.Progress;

import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

public class StsProgressTextPanel extends StsJPanel
{
    int nRows;
    int nCols;

    protected StsJPanel textPanel;
	protected JScrollPane scrollPane;
	public JTextArea textOutput = null;

    public StsProgressTextPanel() { }

    /** create and construct panel.*/
    static public StsProgressTextPanel constructor(int nRows, int nCols)
    {
        try
        {
            StsProgressTextPanel panel = new StsProgressTextPanel(nRows, nCols);
            panel.constructPanel();
            return panel;
        }
        catch(Exception e)
        {
            StsException.outputException("StsProgressPanel.construct() failed.", e, StsException.WARNING);
            return null;
        }
    }

    /** Only create panel; call constructPanel() when other elements and attributes have been assigned. */
    static public StsProgressTextPanel create(int nRows, int nCols)
    {
        try
        {
            return new StsProgressTextPanel(nRows, nCols);
        }
        catch(Exception e)
        {
            StsException.outputException("StsProgressPanel.construct() failed.", e, StsException.WARNING);
            return null;
        }
    }

    private StsProgressTextPanel(int nRows, int nCols)
	{
        super();
        this.nRows = nRows;
        this.nCols = nCols;
    }

    public void constructPanel()
    {
        gbc.insets = new Insets(5, 5, 5, 5);
        UIManager.put("ProgressBar.selectionBackground",Color.BLACK);

		if(nRows > 0)
		{
			textOutput = new JTextArea(nRows, nCols);
			textOutput.setBackground(Color.WHITE);
            textOutput.setFont(new java.awt.Font("Serif", 1, 12));
            textOutput.setToolTipText("");
			textOutput.setLineWrap(false);
			scrollPane = new JScrollPane(textOutput);
//			scrollPane.setPreferredSize(textOutput.getPreferredScrollableViewportSize());
			scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
			scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
            scrollPane.getHorizontalScrollBar().setValue(0);
            scrollPane.getVerticalScrollBar().setValue(0); 
//			textPanel = new StsJPanel();
//			textPanel.gbc.fill = GridBagConstraints.BOTH;
//			textPanel.setBorder(BorderFactory.createEtchedBorder());
//			textPanel.add(scrollPane);
            gbc.fill = GridBagConstraints.BOTH;
            gbc.weighty = 1.0;
            setBorder(BorderFactory.createEtchedBorder());
            add(scrollPane);
//            add(textPanel);
		}
	}

    public synchronized void clearText()
	{
        Runnable clearTextRunnable = new Runnable()
        {
            synchronized public void run()
            {
                textOutput.setText("");
                textOutput.update(textOutput.getGraphics());
            }
        };
        StsToolkit.runLaterOnEventThread(clearTextRunnable);
    }
/*
    public void clearText()
    {
        textOutput.setText("");
        textOutput.update(textOutput.getGraphics());
    }
*/
/*
    public void appendLine(String line)
    {
        if(textOutput == null)return;
        textOutput.append(line + "\n");
        textOutput.update(textOutput.getGraphics());
        Point pt = new Point(0, (int)scrollPane.getViewport().getViewSize().getHeight());
        scrollPane.getViewport().setViewPosition(pt);
    }
*/
    public synchronized void appendLine(String line)
    {
        final String panelLine = line;
        Runnable runnable = new Runnable()
        {
            synchronized public void run()
            {
                if(textOutput == null) return;
                textOutput.append(panelLine + "\n");
//                textOutput.update(textOutput.getGraphics());
                Point pt = new Point(0, (int) scrollPane.getViewport().getViewSize().getHeight());
                scrollPane.getViewport().setViewPosition(pt);
            }
        };
        StsToolkit.runLaterOnEventThread(runnable);
    }

    // TODO Doesn't work yet.  Change design to JTextPane so we can support multiple fonts/colors
    public synchronized void appendErrorLine(String line)
    {
        final String panelLine = line;
        Runnable runnable = new Runnable()
        {
            synchronized public void run()
            {
                if(textOutput == null)return;
                textOutput.setForeground(Color.RED);
                textOutput.append(panelLine + "\n");
                textOutput.setForeground(Color.BLACK);
//                textOutput.update(textOutput.getGraphics());
                Point pt = new Point(0, (int) scrollPane.getViewport().getViewSize().getHeight());
                scrollPane.getViewport().setViewPosition(pt);

            }
        };
        StsToolkit.runLaterOnEventThread(runnable);
    }
    static public void main(String[] args)
	{
        Runnable winThread = new Runnable()
        {
            public void run()
            {
                StsProgressTextPanel textPanel = StsProgressTextPanel.constructor(5, 20);
                ProgressTextTest test =  new ProgressTextTest(textPanel);
                StsGroupBox groupBox = new StsGroupBox("Test Progress Panel");
                StsButton cancelButton = new StsButton("Cancel", "Cancel process", test, "cancel");
                groupBox.add(cancelButton);
                textPanel.add(groupBox);
                textPanel.constructPanel();
                StsToolkit.createDialog(textPanel, false);

		        test.run();
            }
        };
        StsToolkit.runLaterOnEventThread(winThread); 
	}
}

class ProgressTextTest
{
	StsProgressTextPanel textPanel;
    int max = 100;
    boolean cancel = false;

    ProgressTextTest(StsProgressTextPanel textPanel)
	{
        this.textPanel = textPanel;
    }

    public void cancel() { cancel = true; }

    public void run()
	{
		Runnable progressTest = new Runnable()
		{
			public void run()
			{
				runSomething();
			}
		};

		Thread progressTestThread = new Thread(progressTest);
		progressTestThread.start();
	}

    void runSomething()
	{
		try
		{
			textPanel.appendLine("Start message.");
			for(int n = 1; n <= max; n++)
			{
                if(cancel) return;
				textPanel.appendLine("Step " + n);
                StsToolkit.sleep(100);
            }
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}
}

