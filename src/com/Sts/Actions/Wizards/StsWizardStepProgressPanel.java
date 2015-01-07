package com.Sts.Actions.Wizards;

import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

public class StsWizardStepProgressPanel extends StsJPanel
{
    public StsWizard wizard;
    public StsWizardStep wizardStep;
    protected JPanel textPanel;
    protected JScrollPane jScrollPane1;
    protected JTextArea textOutput;
    public JProgressBar progressBar;
    public JLabel progressLbl;

    public int progressBarMaxValue = 0;
    public int progressBarValue = 0;
    public String panelLine = "";

    public StsWizardStepProgressPanel()
    {
//        initializePanel();
    }

    public StsWizardStepProgressPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        setLayout(new GridBagLayout());
    }

    public void addProgressPanel(int nRows, int nCols)
    {
        textPanel = new JPanel(new BorderLayout());
        textPanel.setBorder(BorderFactory.createEtchedBorder());
        textPanel.setDebugGraphicsOptions(0);

        jScrollPane1 = new JScrollPane();
//        scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
//        scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

        textOutput = new JTextArea(nRows, nCols);
        textOutput.setBackground(SystemColor.menu);
        textOutput.setToolTipText("");
        textOutput.append("Built text output.\n");
        textOutput.setLineWrap(false);
//        textOutput.setMinimumSize(new Dimension(100, 100));

        jScrollPane1.getViewport().add(textOutput);
        textPanel.add(jScrollPane1, BorderLayout.CENTER);

        progressBar = new JProgressBar();
        progressLbl = new JLabel();

        gbc.fill = GridBagConstraints.BOTH;
        addEndRow(textPanel);
        gbc.fill = GridBagConstraints.HORIZONTAL;
        addEndRow(progressBar);
        addEndRow(progressLbl);
    }

    public void clearScreen()
    {
        int nLines = 0;
        String line = null;
        initializeProgressBar(0);
        textOutput.setText("");
        textOutput.update(textOutput.getGraphics());
    }

    public void appendLine(String line)
    {
        panelLine = line;
        StsToolkit.runLaterOnEventThread(new Runnable() { public void run() { swingAppendLine(); }});
    }

    public void swingAppendLine()
    {
        textOutput.append(panelLine + "\n");
        textOutput.update(textOutput.getGraphics());
        Point pt = new Point(0, (int) jScrollPane1.getViewport().getViewSize().getHeight());
        jScrollPane1.getViewport().setViewPosition(pt);
    }

    public void initializeProgressBar(int maxValue)
    {
        progressBarMaxValue = maxValue;
        StsToolkit.runLaterOnEventThread(new Runnable() { public void run() { swingInitializeProgressBar(); }});
    }

    public void swingInitializeProgressBar()
    {
//		System.out.println("progress bar initialized to max of 100");
        progressBar.setMaximum(100); // 100 percent
        progressBar.setValue(0);
        progressLbl.setText("0%");
    }

    public void setProgressBarMax(int maxValue)
    {
        progressBarMaxValue = maxValue;
        StsToolkit.runLaterOnEventThread(new Runnable() { public void run() { swingSetProgressBarMax(); }});
    }

    public void swingSetProgressBarMax()
    {
        int percent = (int) (100.0f * (float) progressBarValue / progressBarMaxValue);
        progressBar.setValue(percent);
        progressLbl.setText(percent + "%");
    }

    public void incrementProgress()
    {
        progressBar.setValue(++progressBarValue);
    }

    public void setProgressBarValue(int value)
    {
        progressBarValue = value;
        StsToolkit.runLaterOnEventThread(new Runnable() { public void run() { swingSetProgressBarValue(); }});
    }

    public void swingSetProgressBarValue()
    {
        int percent = (int) (100.0f * (float) progressBarValue / progressBarMaxValue);
//		System.out.println("swingSetProgressBarValue to percent " + percent);
        progressBar.setValue(percent);
        progressLbl.setText(percent + "%");
    }

    public void setProgressBarFinished()
    {
        setProgressBarValue(progressBarMaxValue);
    }
    public void end()
    {
        wizardStep.end();
    }

    static public void main(String[] args)
    {
        JDialog dialog = new JDialog((Frame) null, "Panel Test", true);
        StsWizardStepProgressPanel progressPanel = new StsWizardStepProgressPanel(null, null);
        progressPanel.addProgressPanel(0, 0);
        dialog.getContentPane().add(progressPanel);
        dialog.pack();
        dialog.setVisible(true);
    }
}
