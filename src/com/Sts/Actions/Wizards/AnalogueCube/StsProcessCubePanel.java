package com.Sts.Actions.Wizards.AnalogueCube;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.IO.*;
import com.Sts.SeismicAttributes.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author T.Lasseter
 * @version 1.1
 */

public class StsProcessCubePanel extends StsJPanel implements ActionListener, ChangeListener
{
    StsAnalogueCubeWizard analogueWizard;
    StsWizardStep wizardStep;
    StsProgressPanel progressPanel;
    StsToggleButton dryRunBtn;
    StsButton fullRunBtn;
    boolean runForeground = false;
    float correlThreshold = 0.95f;
//    boolean killIt = false;
    static String PREVIEW = "Preview Mode";
    static String STANDBY = "Standby Mode";
    static String EXECUTE = "Execute Mode";

	String panelLine = "";

    String executeCmd = null;
    JButton killProcessBtn;
    StsSliderBean correlationSlider;
    JPanel executePanel;
    JLabel optMethodLabel;
    JPanel optMethodPanel;
    JRadioButton runningAvg;
    JRadioButton spiralAvg;
    JRadioButton spiralMax;
    ButtonGroup opMethodGrp;
    JLabel modeLabel;

    GridBagLayout gridBagLayout1 = new GridBagLayout();
    static float correlThresholdMin = -1.0f;

    public StsProcessCubePanel(StsAnalogueCubeWizard wizard, StsWizardStep wizardStep)
    {
        analogueWizard = wizard;
        this.wizardStep = wizardStep;
        initializePanel();
    }

    public void initializePanel()
    {
        this.setLayout(gridBagLayout1);

        correlationSlider = new StsSliderBean(false);
        correlationSlider.initSliderValues(correlThresholdMin, 1.0f, 0.01f, correlThreshold);
        correlationSlider.addChangeListener(this);
        correlationSlider.setValueLabel("Correl Min");
        correlationSlider.showCheckbox(false);

        optMethodLabel = new JLabel();
        optMethodLabel.setFont(new java.awt.Font("Dialog", 0, 11));
        optMethodLabel.setHorizontalAlignment(SwingConstants.LEFT);
        optMethodLabel.setText("Optimization Method:");

        runningAvg = new JRadioButton();
        runningAvg.setText("Running Avg.");
        runningAvg.addActionListener(this);

        spiralAvg = new JRadioButton();
        spiralAvg.setSelected(true);
        spiralAvg.addActionListener(this);
        spiralAvg.setText("Spiral Avg.");

        spiralMax = new JRadioButton();
        spiralMax.setText("Spiral Max.");
        spiralMax.addActionListener(this);

        opMethodGrp = new ButtonGroup();
        opMethodGrp.add(runningAvg);
        opMethodGrp.add(spiralAvg);
        opMethodGrp.add(spiralMax);

        optMethodPanel = new JPanel();
        optMethodPanel.setLayout(new GridBagLayout());
        optMethodPanel.add(optMethodLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(3, 1, 0, 0), 0, 0));
        optMethodPanel.add(runningAvg, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 6, 0, 0), 0, -5));
        optMethodPanel.add(spiralAvg, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
            , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 0, 0, 0), 7, -3));
        optMethodPanel.add(spiralMax, new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
            , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 0, 0, 0), 7, -2));

        modeLabel = new JLabel();
        modeLabel.setFont(new java.awt.Font("Serif", 1, 16));
        modeLabel.setForeground(Color.blue);
        modeLabel.setHorizontalAlignment(SwingConstants.CENTER);
        modeLabel.setText("Stand-By Mode");

        dryRunBtn = new StsToggleButton(" Preview ", "Runs analog cube on current cursor plane.", this, "dryRun");
        fullRunBtn = new StsButton(" Execute ", "Keep wizard open and run.", this, "fullRun");
        killProcessBtn = new StsButton(" Cancel ", "Kill the current process", this, "killRun");
        killProcessBtn.setEnabled(false);

        executePanel = new JPanel(new GridBagLayout());
        executePanel.setBorder(null);
        executePanel.setDebugGraphicsOptions(0);
		executePanel.add(modeLabel,   new GridBagConstraints(0, 0, 3, 1, 0.0, 0.0
				,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 1, 0, 1), 274, 0));
        executePanel.add(dryRunBtn,  new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 1, 3, 0), 0, 0));
        executePanel.add(fullRunBtn,  new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 5, 3, 5), 2, 0));
        executePanel.add(killProcessBtn,  new GridBagConstraints(2, 1, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 0, 3, 4), 8, 0));


        progressPanel = StsProgressPanel.constructorWithCancelButton(5, 50);
        this.add(progressPanel,  new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
                ,GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(3, 2, 0, 0), 0, 0));
        this.add(optMethodPanel,  new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
                ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(3, 2, 0, 0), 0, 0));
        this.add(executePanel,  new GridBagConstraints(0, 2, 1, 1, 1.0, 0.0
                ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(3, 2, 0, 0), 0, 0));
    }

    public void killRun()
    {
        fullRunBtn.setEnabled(true);
        dryRunBtn.setSelected(false);
        killProcessBtn.setEnabled(true);
        wizardStep.enableFinish();
        wizardStep.enableCancel();
        wizardStep.enablePrevious();

        progressPanel.appendLine("Killing Process...");
		progressPanel.initialize(0);
		analogueWizard.killRun();
//        killIt = true;
        analogueWizard.endDataPlaneCalculation();
        dryRunComplete();
    }
/*
    public boolean isKilled()
    {
        return killIt;
    }
*/
    public void dryRun()
    {
//        killIt = false;
        dryRunStart();
        boolean isDryRun = dryRunBtn.isSelected();
        analogueWizard.dryRun(this, isDryRun);
        modeLabel.setText(PREVIEW);
    }

    public void fullRun()
    {
//        killIt = false;
        executionStart();
        analogueWizard.fullRun(this);
    }

    public void executionStart()
    {
//        dryRunBtn.setSelected(false);
        dryRunBtn.setEnabled(false);
        killProcessBtn.setEnabled(true);
        fullRunBtn.setEnabled(false);
        wizardStep.disableFinish();
        wizardStep.disableCancel();
        wizardStep.disablePrevious();
        modeLabel.setText(EXECUTE);
    }

    public void executionComplete()
    {
//        dryRunBtn.setSelected(false);
        dryRunBtn.setEnabled(true);
        killProcessBtn.setEnabled(false);
        fullRunBtn.setEnabled(true);
		fullRunBtn.setSelected(false);
        wizardStep.enableFinish();
        wizardStep.disableCancel();
        wizardStep.enablePrevious();
        modeLabel.setText(STANDBY);
    }

    public void buttonsReinitialize()
    {
        killProcessBtn.setEnabled(false);
        fullRunBtn.setEnabled(true);
        dryRunBtn.setEnabled(true);
		dryRunBtn.setSelected(false);
        fullRunBtn.setSelected(false);
        wizardStep.enableFinish();
        wizardStep.enablePrevious();
        wizardStep.enableNext();
    }

    public void resetDryRunBtn()
    {
		dryRunBtn.setSelected(false);
        modeLabel.setText(STANDBY);
    }
    public void dryRunComplete()
    {
        killProcessBtn.setEnabled(false);
        fullRunBtn.setEnabled(true);
        dryRunBtn.setEnabled(true);
//        dryRunBtn.setSelected(false);
		 fullRunBtn.setSelected(false);

        wizardStep.enableFinish();
        wizardStep.enableCancel();
        wizardStep.enablePrevious();
        modeLabel.setText(STANDBY);
    }

    public void dryRunStart()
    {
        killProcessBtn.setEnabled(true);
        fullRunBtn.setEnabled(true);
        wizardStep.disableFinish();
        wizardStep.disableCancel();
        wizardStep.disablePrevious();
        modeLabel.setText(PREVIEW);
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();
    }

    public void displayHeader(String cmd)
    {
        StsAnalogueVolumeConstructor analogueVolumeConstructor = analogueWizard.getAnalogueVolume();
        if (analogueVolumeConstructor == null)
        {
            return;
        }
        String name = analogueVolumeConstructor.getVolume().getName();
        progressPanel.appendLine("Executing: " + cmd + "\n\n");
 //       progressPanel.textOutput.update(progressPanel.textOutput.getGraphics());
        progressPanel.appendLine("                             Name: " + name + "\n");
        progressPanel.appendLine("                    Target PostStack3d: " + analogueVolumeConstructor.getTargetVolume().getName() + "\n");
        progressPanel.appendLine("                    Source PostStack3d: " + analogueVolumeConstructor.getSourceVolume().getName() + "\n");
        if (analogueVolumeConstructor.getTargetSubVolume() != null)
        {
            progressPanel.appendLine("            SubVolume Constraints: " + analogueVolumeConstructor.getTargetSubVolume().getName() +
                                "\n");
        }
        else
        {
            progressPanel.appendLine("            SubVolume Constraints: None\n");

        }
        if (analogueVolumeConstructor.isDataFloat()) progressPanel.appendLine("              Analysis Resolution: 32-bit float\n");
        else progressPanel.appendLine("              Analysis Resolution: 8-bit\n");

		progressPanel.appendLine("                       is complex: " + analogueVolumeConstructor.isIsComplex() + "\n");
		progressPanel.appendLine("           Correlation Threashold: " + analogueWizard.getCorrelThreshold() + "\n\n");

		StsSeismicBoundingBox seismicBoundingBox = analogueVolumeConstructor.getVolume();
		StsSubVolume targetSubVolume = analogueVolumeConstructor.getTargetSubVolume();
		if(targetSubVolume != null)
		{
			StsGridBoundingBox targetBoundingBox = analogueVolumeConstructor.getTargetSubVolume().getGridBoundingBox();
			if (targetBoundingBox != null)
			{
				progressPanel.appendLine("                       Target Min Line: " + seismicBoundingBox.getNearestBoundedRowNumFromY(targetBoundingBox.yMin) + "\n");
				progressPanel.appendLine("                      Target Max Line: " + seismicBoundingBox.getNearestBoundedRowNumFromY(targetBoundingBox.yMax) + "\n");
				progressPanel.appendLine("                     Target Min XLine: " + seismicBoundingBox.getNearestBoundedColNumFromX(targetBoundingBox.xMin) + "\n");
				progressPanel.appendLine("                     Target Max XLine: " + seismicBoundingBox.getNearestBoundedColNumFromX(targetBoundingBox.xMax) + "\n");
				progressPanel.appendLine("                         Target Min Z: " + targetBoundingBox.getZMin() + "\n");
				progressPanel.appendLine("                         Target Max Z: " + targetBoundingBox.getZMax() + "\n");
			}
		}
		else
		{
			progressPanel.appendLine("                      Target Min Line: " + seismicBoundingBox.getRowNumMin() + "\n");
			progressPanel.appendLine("                      Target Max Line: " + seismicBoundingBox.getRowNumMax() + "\n");
			progressPanel.appendLine("                     Target Min XLine: " + seismicBoundingBox.getColNumMin() + "\n");
			progressPanel.appendLine("                     Target Max XLine: " + seismicBoundingBox.getColNumMax() + "\n");
			progressPanel.appendLine("                         Target Min Z: " + seismicBoundingBox.getZMin() + "\n");
			progressPanel.appendLine("                         Target Max Z: " + seismicBoundingBox.getZMax() + "\n");
		}
		StsBoxSubVolume sourceSubVolume = analogueVolumeConstructor.getSourceSubVolume();
		if (sourceSubVolume != null)
		{
			progressPanel.appendLine("                      Source Min Line: " + sourceSubVolume.getRowNumMin() + "\n");
			progressPanel.appendLine("                      Source Max Line: " + sourceSubVolume.getRowNumMax() + "\n");
			progressPanel.appendLine("                     Source Min XLine: " + sourceSubVolume.getColNumMin() + "\n");
			progressPanel.appendLine("                     Source Max XLine: " + sourceSubVolume.getColNumMax() + "\n");
			progressPanel.appendLine("                         Source Min Z: " + sourceSubVolume.getZMin() + "\n");
			progressPanel.appendLine("                         Source Max Z: " + sourceSubVolume.getZMax() + "\n");
		}
//        progressPanel.textOutput.update(progressPanel.textOutput.getGraphics());
    }

    private void displayFile()
    {
        StsSeismicVolume attributeVolume = analogueWizard.getAttributeVolume();

        StsFile file = StsFile.constructor(attributeVolume.stsDirectory + "AnalogueCube.log");
        StsAsciiFile asciiFile = new StsAsciiFile(file);
        if(asciiFile.openReadWithErrorMessage()) return;
        int nLines = 0;
        String line = null;
//        progressPanel.textOutput.setText("");
//        progressPanel.textOutput.update(progressPanel.textOutput.getGraphics());
        try
        {
            line = asciiFile.readLine();
            while (line != null)
            {
                line = line + "\n";
                progressPanel.appendLine(line);
                line = asciiFile.readLine();
            }
//            progressPanel.textOutput.update(progressPanel.textOutput.getGraphics());
        }
        catch (Exception g)
        {
            System.out.println("Error reading " + attributeVolume.stsDirectory + "AnalogueCube.log");
        }
        asciiFile.close();
    }

    public boolean wasRunForeground()
    {
        return runForeground;
    }

    public void stateChanged(ChangeEvent e)
    {
        Object source = e.getSource();
        //
        // If any of the search criterion is changed the preview mode must be exitted and
        // the user will need to re-start it.
        if (source instanceof StsSliderBean)
        {
            correlThreshold = ( (StsSliderBean) source).getValue();
            analogueWizard.setCorrelThreshold(correlThreshold);
        }
        else if(source instanceof JRadioButton)
        {
            correlThreshold = ((StsSliderBean) source).getValue();
            analogueWizard.setCorrelThreshold(correlThreshold);
        }
    }

    public float getCorrelationThreshold()
    {
        return correlThreshold;
    }

    public byte getOptimizationMethod()
    {
        if (runningAvg.isSelected())
        {
            return analogueWizard.RUNNING_AVERAGE;
        }
        if (spiralAvg.isSelected())
        {
            return analogueWizard.SPIRAL_AVERAGE;
        }
        if (spiralMax.isSelected())
        {
            return analogueWizard.SPIRAL_MAXIMUM;
        }
        else
        {
            return analogueWizard.SPIRAL_AVERAGE;
        }
    }

}
